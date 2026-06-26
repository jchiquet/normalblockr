#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_base.h"
#include "zi_normal_block_data.h"
#include "zi_closed_form_solvers.h"
#include "utils_arma.h"

// Zero-inflated normal-block model with known clusters (fixed indicator
// matrix C), templated on the residual-noise policy (ZIDiagonalNoise /
// ZISphericalNoise, see zi_noise_models.h). Equivalent of the R6 class
// ZINormalBlockKnownClusters (R/ZINormalBlockKnownClusters.R).
//
// Unlike the non zero-inflated counterpart (NormalBlockKnownClusters), the
// posterior covariance of W | Y is *row-dependent* (the zero-inflation mask
// varies row by row), hence Gamma_ is a cube (q x q x n) of n distinct q x q
// matrices rather than a single shared matrix. B's normal equations no
// longer share a common cross-product across columns (the zero-inflation
// mask reweights each column independently), so it is solved one column at
// a time (see zi_closed_form_solvers.h).
template <typename NoisePolicy>
class ZINormalBlockKnownClusters : public NormalBlockBase {
  const ZINormalBlockData& zi_data_;
  arma::mat C_;       // p x q, fixed cluster-indicator matrix
  arma::cube Gamma_;  // q x q x n, posterior covariance of W | Y (one matrix per row)
  arma::mat Mu_;      // n x q, posterior mean of W | Y
  arma::mat R_;       // n x p, Y - X*B computed at the start of the (V)EM step,
                       // reused unchanged by both E_step() and M_step()
  arma::mat dm1_mat_; // n x p, repmat(dm1_, n) % zeros_bar, computed once in
                      // E_step() and reused in M_step() (dm1_ does not change
                      // in between -- it is only updated at the end of M_step())
  arma::mat sum_slices() const {
    arma::mat out(q_, q_, arma::fill::zeros);
    for (arma::uword i = 0; i < Gamma_.n_slices; ++i) out += Gamma_.slice(i);
    return out;
  }

  void E_step() override {
    R_ = zi_data_.Y - zi_data_.X * B_;
    dm1_mat_ = arma::repmat(dm1_.t(), zi_data_.n, 1) % zi_data_.zeros_bar;
    arma::mat dm1C = dm1_mat_ * C_;             // n x q
    arma::mat Rdm1C = (R_ % dm1_mat_) * C_;     // n x q

    for (arma::uword i = 0; i < zi_data_.n; ++i) {
      Gamma_.slice(i) = arma::inv_sympd(Omegaq_ + arma::diagmat(dm1C.row(i).t()));
      Mu_.row(i) = Rdm1C.row(i) * Gamma_.slice(i);
    }
  }

  void M_step() override {
    arma::mat MuCT = Mu_ * C_.t();
    B_ = nb_optim::solve_wls(dm1_mat_, zi_data_.Y, zi_data_.X, MuCT);
    R_ = zi_data_.Y - zi_data_.X * B_;

    arma::mat RmMuCT = R_ - MuCT;
    arma::mat CgC(zi_data_.n, zi_data_.p);
    for (arma::uword i = 0; i < zi_data_.n; ++i) CgC.row(i) = (C_ * Gamma_.slice(i).diag()).t();
    arma::mat A = arma::square(RmMuCT) + CgC;

    arma::vec weighted_ssq = arma::vectorise(arma::sum(zi_data_.zeros_bar % A, 0));
    dm1_ = NoisePolicy::update_dm1(weighted_ssq, zi_data_.nY, zi_data_.npY);

    arma::mat Sigma_hat = (Mu_.t() * Mu_ + sum_slices()) / zi_data_.n;
    Omegaq_ = estimate_omega(Sigma_hat);
  }

public:
  ZINormalBlockKnownClusters(const ZINormalBlockData& data, const arma::mat& C,
                              const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                              double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockBase(data, C.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    zi_data_(data), C_(C), Gamma_(C.n_cols, C.n_cols, data.n), Mu_(arma::zeros(data.n, C.n_cols)) {
    for (arma::uword i = 0; i < data.n; ++i) Gamma_.slice(i) = arma::eye(C.n_cols, C.n_cols);
  }

  // General (non-profiled) marginal log-likelihood of Y given the fixed
  // zero-inflation mask, valid at *any* (B_, dm1_, Omegaq_) -- see
  // NormalBlockKnownClusters::objective() for the non-ZI derivation this
  // mirrors. Because the zero-inflation mask varies row by row, each row i
  // has its own marginal: Y_i,obs ~ N(B^T X_i,obs, Sigma_Y,i), restricted to
  // i's *observed* (non-zero-inflated) columns, with Sigma_Y,i = D_i,obs +
  // C_obs Omegaq^{-1} C_obs^T. The same Woodbury/determinant-lemma trick
  // applies per row, through the row's own posterior precision
  // Gamma^{(i),-1} = Omegaq + diag((W C)_{i.}) (the matrix E_step() already
  // inverts to get Gamma_.slice(i)/Mu_.row(i) -- recomputed fresh here, see
  // the non-ZI objective() for why):
  //
  //   log|Sigma_Y,i|                 = -sum_{j obs} log(dm1_j) - log|Omegaq| - log|Gamma^{(i)}|
  //   (R_i,obs)^T Sigma_Y,i^{-1} R_i,obs = sum_{j obs} dm1_j R_ij^2 - mu_i^T Gamma^{(i),-1} mu_i
  //
  // summed over rows (the per-row log(dm1) sum becomes the usual nY-weighted
  // sum once accumulated across all n rows). This replaces the same kind of
  // profiled shortcut (and the same sign bug, +log|Gamma^{(i),-1}| instead
  // of +log|Gamma^{(i)}|) described in NormalBlockKnownClusters::objective().
  double objective() const override {
    arma::mat R = zi_data_.Y - zi_data_.X * B_;
    arma::mat dm1_mat = arma::repmat(dm1_.t(), zi_data_.n, 1) % zi_data_.zeros_bar;
    arma::mat dm1C = dm1_mat * C_;          // n x q
    arma::mat Rdm1C = (R % dm1_mat) * C_;   // n x q

    double sum_log_Gamma_inv = 0.0;
    double quad_mu = 0.0;
    for (arma::uword i = 0; i < zi_data_.n; ++i) {
      arma::mat Gamma_inv_i = Omegaq_ + arma::diagmat(dm1C.row(i).t());
      auto Gamma_fresh_i = nb_utils::inv_and_log_det_sympd(Gamma_inv_i);
      sum_log_Gamma_inv += Gamma_fresh_i.log_det;
      arma::vec mu_i = (Rdm1C.row(i) * Gamma_fresh_i.inv).t();
      quad_mu += arma::as_scalar(mu_i.t() * Gamma_inv_i * mu_i);
    }

    double log_det_Omegaq = arma::log_det_sympd(Omegaq_);
    double weighted_sum_log_dm1 = arma::accu(zi_data_.nY % arma::log(dm1_));
    double SSQ_w = arma::accu(dm1_mat % arma::square(R));

    double J = -0.5 * zi_data_.npY * std::log(2.0 * arma::datum::pi);
    J += 0.5 * weighted_sum_log_dm1;
    J += 0.5 * zi_data_.n * log_det_Omegaq;
    J -= 0.5 * sum_log_Gamma_inv; // -log|Gamma^{-1}| = +log|Gamma|, summed over rows
    J -= 0.5 * (SSQ_w - quad_mu);
    J += zi_data_.zi_cond_mean;
    return J;
  }

  const arma::cube& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }

  std::unique_ptr<NormalBlockBase> clone() const override {
    return std::make_unique<ZINormalBlockKnownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockBase& other) override {
    copy_tracked_state_from(other);
    const auto& o = static_cast<const ZINormalBlockKnownClusters<NoisePolicy>&>(other);
    Gamma_ = o.Gamma_;
    Mu_ = o.Mu_;
    R_ = o.R_;
    dm1_mat_ = o.dm1_mat_;
  }
  // Same general/always-valid objective() as NormalBlockKnownClusters (see
  // its supports_acceleration()), so the same objective-comparison-gated
  // SQUAREM applies here too: the per-row precision matrices
  // Omegaq + diag(dm1C.row(i)) add a *non-negative* diagonal to Omegaq (dm1
  // > 0, C is 0/1), which by Weyl's inequality can only raise -- never
  // lower -- the smallest eigenvalue relative to Omegaq's own. So whenever
  // the extrapolated Omegaq candidate itself passes state_is_feasible(),
  // every row's matrix is at least as well-conditioned; there is no
  // separate per-row numerical risk left for state_is_feasible() to miss,
  // unlike under the earlier conditioning-only safety net this gate used to
  // rely on.
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_KNOWN_CLUSTERS_H
