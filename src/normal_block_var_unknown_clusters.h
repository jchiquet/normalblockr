#ifndef NORMALBLOCKR_NORMAL_BLOCK_VAR_UNKNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_VAR_UNKNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_var_base.h"
#include "normal_block_data.h"
#include "utils_arma.h"

// Normal-block model with unknown clusters, inferred through a variational
// EM (mean-field on W and on the cluster membership C), templated on the
// residual-noise policy (see noise_models.h). Equivalent of the R6 class
// NormalBlockVarUnknownClusters (R/NormalBlockVarUnknownClusters.R); see
// inst/normal_block_models.qmd §3/§4 for the VEM recursion.
template <typename NoisePolicy>
class NormalBlockVarUnknownClusters : public NormalBlockVarBase {
  arma::mat C_;       // p x q, variational membership probabilities (tau)
  arma::vec alpha_;   // q, prior cluster probabilities
  arma::mat M_;       // n x q, variational mean of W | Y
  arma::vec S_;       // q, variational variance of W | Y (constant across rows)
  arma::mat Gamma_;   // q x q, posterior covariance of W | Y (not part of the model state)
  bool fixed_tau_;    // if true, C/alpha are not re-estimated (stability selection)
  arma::mat R_;       // n x p, Y - X*B computed at the start of the (V)EM step,
                      // reused unchanged by both E_step() and M_step()
  arma::mat Sigma_hat_; // q x q, M-step covariance estimate, cached for objective()
                        // (recomputing it there would be identical: nothing
                        // changes M_/S_ between M_step() and objective())

  void E_step() override {
    R_ = data_.Y - XB();

    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    Gamma_ = arma::inv_sympd(Omega_ + arma::diagmat(C_.t() * dm1_));
    M_ = R_ * dm1C * Gamma_;
    S_ = Gamma_.diag();

    if (q_ > 1 && !fixed_tau_) {
      arma::mat RtM = R_.t() * M_;
      RtM.each_col() %= dm1_;
      arma::vec rowfactor = arma::vectorise(arma::sum(arma::square(M_), 0)) + data_.n * S_;
      arma::mat term2 = dm1_ * rowfactor.t();
      arma::mat term3(data_.p, q_);
      term3.each_row() = arma::log(alpha_).t();
      arma::mat eta = RtM - 0.5 * term2 + term3 - 1.0;
      C_ = nb_utils::clip_probabilities(nb_utils::softmax_rows(eta));
    }
  }

  void M_step() override {
    arma::mat MCT = M_ * C_.t();
    set_B(data_.XtXm1 * (data_.X.t() * (data_.Y - MCT)));
    R_ = data_.Y - XB();

    arma::mat M2plusS = arma::square(M_);
    M2plusS.each_row() += S_.t();
    arma::mat term = M2plusS * C_.t();
    arma::vec ddiag = arma::vectorise(arma::mean(arma::square(R_) - 2.0 * (R_ % MCT) + term, 0));
    dm1_ = NoisePolicy::update_dm1(ddiag);

    alpha_ = arma::vectorise(arma::mean(C_, 0));
    Sigma_hat_ = M_.t() * M_ / data_.n + arma::diagmat(S_);
    Omega_ = estimate_omega(Sigma_hat_);
  }

public:
  NormalBlockVarUnknownClusters(const NormalBlockData& data,
                              const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                              const arma::mat& C0, const arma::vec& alpha0,
                              const arma::mat& M0, const arma::vec& S0,
                              double sparsity, const arma::mat& sparsity_weights,
                              bool fixed_tau = false) :
    NormalBlockVarBase(data, C0.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    C_(C0), alpha_(alpha0), M_(M0), S_(S0), fixed_tau_(fixed_tau) {
    Sigma_hat_ = M_.t() * M_ / data.n + arma::diagmat(S_); // matches the very first objective() call in run_em()
  }

  // ELBO at the current parameter values (see inst/normal_block_models.qmd
  // §3). When sparsity_ > 0, Omega_ is no longer the exact inverse of
  // Sigma_hat, so the trace/penalty correction that otherwise cancels out
  // must be added back.
  double objective() const override {
    double log_det_Omega = arma::log_det_sympd(Omega_);
    double sum_log_dm1    = arma::sum(arma::log(dm1_));
    double sum_log_S      = arma::sum(arma::log(S_));
    double log2pie        = 1.0 + std::log(2.0 * arma::datum::pi);

    double J = -0.5 * data_.n * data_.p * log2pie + 0.5 * data_.n * sum_log_dm1;
    J += 0.5 * data_.n * log_det_Omega + 0.5 * data_.n * sum_log_S;
    J += arma::accu(C_ * arma::log(alpha_)) - nb_utils::sum_xlogx(C_);

    if (sparsity_ > 0.0) {
      J += 0.5 * data_.n * q_ - 0.5 * data_.n * arma::trace(Omega_ * Sigma_hat_);
      J -= sparsity_ * arma::accu(arma::abs(sparsity_weights_ % Omega_));
    }
    return J;
  }

  const arma::mat& C() const { return C_; }
  const arma::vec& alpha() const { return alpha_; }
  const arma::mat& M() const { return M_; }
  const arma::vec& S() const { return S_; }

  std::unique_ptr<NormalBlockVarBase> clone() const override {
    return std::make_unique<NormalBlockVarUnknownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockVarBase& other) override {
    copy_tracked_state_from(other);
    const auto& o = static_cast<const NormalBlockVarUnknownClusters<NoisePolicy>&>(other);
    C_ = o.C_;
    alpha_ = o.alpha_;
    M_ = o.M_;
    S_ = o.S_;
    Gamma_ = o.Gamma_;
    R_ = o.R_;
    Sigma_hat_ = o.Sigma_hat_;
  }
  // Same SQUAREM mechanism as NormalBlockVarKnownClusters, but relies on a
  // different validity argument: every objective() comparison happens right
  // after a fresh E_step()+M_step() pair, exactly the regime the profiled
  // ELBO shortcut is already valid in. See inst/normal_block_models.qmd
  // ("SQUAREM") for the speedup measured on real data.
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_UNKNOWN_CLUSTERS_H
