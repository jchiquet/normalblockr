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
  arma::mat Sigma_hat_; // q x q, M-step covariance estimate, cached for objective()
                        // (recomputing it there would be identical: nothing
                        // changes Mu_/Gamma_ between M_step() and objective())
  arma::vec log_det_Gamma_; // n, log-determinant of each Gamma_ slice, a free
                            // byproduct of the per-row Cholesky factorization
                            // already needed to invert it in E_step() --
                            // cached for objective() instead of running n
                            // more factorizations via log_det_sympd()

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
      auto Gamma_inv = nb_utils::inv_and_log_det_sympd(Omegaq_ + arma::diagmat(dm1C.row(i).t()));
      Gamma_.slice(i) = Gamma_inv.inv;
      log_det_Gamma_(i) = Gamma_inv.log_det;
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

    Sigma_hat_ = (Mu_.t() * Mu_ + sum_slices()) / zi_data_.n;
    Omegaq_ = estimate_omega(Sigma_hat_);
  }

public:
  ZINormalBlockKnownClusters(const ZINormalBlockData& data, const arma::mat& C,
                              const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                              double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockBase(data, C.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    zi_data_(data), C_(C), Gamma_(C.n_cols, C.n_cols, data.n), Mu_(arma::zeros(data.n, C.n_cols)) {
    for (arma::uword i = 0; i < data.n; ++i) Gamma_.slice(i) = arma::eye(C.n_cols, C.n_cols);
    Sigma_hat_ = (Mu_.t() * Mu_ + sum_slices()) / data.n; // matches the very first objective() call in run_em()
    log_det_Gamma_ = arma::zeros(data.n); // log_det(I) = 0 for every row
  }

  // Log-likelihood criterion, mirrors compute_loglik in
  // R/ZINormalBlockKnownClusters.R: same shape as the non zero-inflated
  // criterion (NormalBlockKnownClusters::objective), but normalized by
  // nY/npY instead of n/p (only non-zero-inflated observations contribute
  // to the Gaussian part), plus the fixed zero-inflation contribution
  // (zi_cond_mean), and one log-determinant per row (Gamma_ is per-row).
  double objective() const override {
    double log_det_Omegaq = arma::log_det_sympd(Omegaq_);
    double weighted_sum_log_dm1 = arma::accu(zi_data_.nY % arma::log(dm1_));
    double log2pie = 1.0 + std::log(2.0 * arma::datum::pi);

    double J = -0.5 * zi_data_.npY * log2pie + 0.5 * weighted_sum_log_dm1;
    J += 0.5 * zi_data_.n * log_det_Omegaq + 0.5 * arma::accu(log_det_Gamma_);
    J += zi_data_.zi_cond_mean;

    if (sparsity_ > 0.0) {
      J += 0.5 * zi_data_.n * q_ - 0.5 * zi_data_.n * arma::trace(Omegaq_ * Sigma_hat_);
      J -= sparsity_ * arma::accu(arma::abs(sparsity_weights_ % Omegaq_));
    }
    return J;
  }

  const arma::cube& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_KNOWN_CLUSTERS_H
