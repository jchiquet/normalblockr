#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_UNKNOWN_CLUSTERS_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_UNKNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_base.h"
#include "zi_normal_block_data.h"
#include "utils_arma.h"
#include "zi_closed_form_solvers.h"

// Zero-inflated normal-block model with unknown clusters, templated on the
// residual-noise policy (ZIDiagonalNoise / ZISphericalNoise, see
// zi_noise_models.h). Equivalent of the R6 class ZINormalBlockUnknownClusters
// (R/ZINormalBlockUnknownClusters.R).
//
// Unlike the non zero-inflated counterpart (NormalBlockUnknownClusters), the
// variational variance of W | Y is *row-dependent* (the zero-inflation mask
// varies row by row), hence S_ is a matrix (n x q) rather than a vector of
// length q. The zero-inflation mask reweights B's normal equations
// per-column and M's per-row, but both subproblems stay exactly quadratic,
// so both are solved by direct linear systems rather than an iterative
// optimizer (see zi_closed_form_solvers.h).
template <typename NoisePolicy>
class ZINormalBlockUnknownClusters : public NormalBlockBase {
  const ZINormalBlockData& zi_data_;
  arma::mat C_;       // p x q, variational membership probabilities (tau)
  arma::vec alpha_;   // q, prior cluster probabilities
  arma::mat M_;       // n x q, variational mean of W | Y
  arma::mat S_;       // n x q, variational variance of W | Y (row-dependent)
  bool fixed_tau_;    // if true, C/alpha are not re-estimated (stability selection)
  arma::mat R_;       // n x p, Y - X*B computed at the start of the (V)EM step,
                       // reused unchanged by both E_step() and M_step()
  arma::mat DM1_;     // n x p, repmat(dm1_, n) % zeros_bar, computed once in
                      // E_step() and reused in M_step() (dm1_ does not change
                      // in between -- it is only updated at the end of M_step())
  arma::mat Sigma_hat_; // q x q, M-step covariance estimate, cached for objective()
                        // (recomputing it there would be identical: nothing
                        // changes M_/S_ between M_step() and objective())

  void E_step() override {
    R_ = zi_data_.Y - zi_data_.X * B_;
    DM1_ = arma::repmat(dm1_.t(), zi_data_.n, 1) % zi_data_.zeros_bar;

    M_ = nb_optim::solve_tau_ridge(DM1_, R_, C_, Omegaq_);

    arma::mat DM1C = DM1_ * C_;                // n x q
    DM1C.each_row() += Omegaq_.diag().t();
    S_ = 1.0 / DM1C;

    if (q_ > 1 && !fixed_tau_) {
      arma::mat term1 = -0.5 * (DM1_.t() * (arma::square(M_) + S_));  // p x q
      arma::mat term2 = (DM1_ % R_).t() * M_;                          // p x q
      arma::mat eta = term1 + term2;
      eta.each_row() += (arma::log(alpha_) - 1.0).t();
      C_ = nb_utils::clip_probabilities(nb_utils::softmax_rows(eta));
    }
  }

  void M_step() override {
    arma::mat MCT = M_ * C_.t();
    B_ = nb_optim::solve_wls(DM1_, zi_data_.Y, zi_data_.X, MCT);
    R_ = zi_data_.Y - zi_data_.X * B_;

    arma::mat A = arma::square(R_) - 2.0 * (R_ % MCT) + (arma::square(M_) + S_) * C_.t();
    arma::vec weighted_ssq = arma::vectorise(arma::sum(zi_data_.zeros_bar % A, 0));
    dm1_ = NoisePolicy::update_dm1(weighted_ssq, zi_data_.nY, zi_data_.npY);

    alpha_ = arma::vectorise(arma::mean(C_, 0));
    Sigma_hat_ = M_.t() * M_ / zi_data_.n + arma::diagmat(arma::vectorise(arma::mean(S_, 0)));
    Omegaq_ = estimate_omega(Sigma_hat_);
  }

public:
  ZINormalBlockUnknownClusters(const ZINormalBlockData& data,
                                const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                                const arma::mat& C0, const arma::vec& alpha0,
                                const arma::mat& M0, const arma::mat& S0,
                                double sparsity, const arma::mat& sparsity_weights,
                                bool fixed_tau = false) :
    NormalBlockBase(data, C0.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    zi_data_(data), C_(C0), alpha_(alpha0), M_(M0), S_(S0), fixed_tau_(fixed_tau) {
    Sigma_hat_ = M_.t() * M_ / data.n + arma::diagmat(arma::vectorise(arma::mean(S_, 0))); // matches the very first objective() call in run_em()
  }

  // ELBO at the current parameter values, mirrors compute_loglik in
  // R/ZINormalBlockUnknownClusters.R: same shape as the non zero-inflated criterion
  // (NormalBlockUnknownClusters::objective), but normalized by nY/npY
  // instead of n/p, plus the fixed zero-inflation contribution
  // (zi_cond_mean); sum(log(S)) ranges over all n*q entries since S is now
  // row-dependent.
  double objective() const override {
    double log_det_Omegaq = arma::log_det_sympd(Omegaq_);
    double weighted_sum_log_dm1 = arma::accu(zi_data_.nY % arma::log(dm1_));
    double sum_log_S = arma::accu(arma::log(S_));
    double log2pie = 1.0 + std::log(2.0 * arma::datum::pi);

    double J = -0.5 * zi_data_.npY * log2pie + 0.5 * weighted_sum_log_dm1;
    J += 0.5 * zi_data_.n * log_det_Omegaq + 0.5 * sum_log_S;
    J += arma::accu(C_ * arma::log(alpha_)) - nb_utils::sum_xlogx(C_);
    J += zi_data_.zi_cond_mean;

    if (sparsity_ > 0.0) {
      J += 0.5 * zi_data_.n * q_ - 0.5 * zi_data_.n * arma::trace(Omegaq_ * Sigma_hat_);
      J -= sparsity_ * arma::accu(arma::abs(sparsity_weights_ % Omegaq_));
    }
    return J;
  }

  const arma::mat& C() const { return C_; }
  const arma::vec& alpha() const { return alpha_; }
  const arma::mat& M() const { return M_; }
  const arma::mat& S() const { return S_; }

  std::unique_ptr<NormalBlockBase> clone() const override {
    return std::make_unique<ZINormalBlockUnknownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockBase& other) override {
    copy_tracked_state_from(other);
    const auto& o = static_cast<const ZINormalBlockUnknownClusters<NoisePolicy>&>(other);
    C_ = o.C_;
    alpha_ = o.alpha_;
    M_ = o.M_;
    S_ = o.S_;
    R_ = o.R_;
    DM1_ = o.DM1_;
    Sigma_hat_ = o.Sigma_hat_;
  }
  // Same mechanism as NormalBlockUnknownClusters (see its
  // supports_acceleration()): no sign bug in this profiled ELBO (same "+1"
  // shortcut structure), and every objective() call inside
  // try_squarem_step() happens right after a fresh VE-step+M-step pair, the
  // exact regime that shortcut is valid in. Validated on synthetic
  // zero-inflated data (see git history): clean, growing speedup with q
  // (e.g. q=25: 2378 plain VEM iterations vs 279 accelerated), and the only
  // non-monotonicity observed (q=2, q=5) is a pre-existing artifact of this
  // model's first iteration specifically (reproduces identically with
  // acceleration off), unrelated to this gate.
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_UNKNOWN_CLUSTERS_H
