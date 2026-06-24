#ifndef NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_base.h"
#include "normal_block_data.h"

// Normal-block model with known clusters (fixed indicator matrix C), templated
// on the residual-noise policy (DiagonalNoise / SphericalNoise, see
// noise_models.h). Equivalent of the R6 class NormalBlockKnownClusters
// (R/NormalBlockKnownClusters.R), implementing the EM recursion of section 2
// (diagonal D) / section 3 (spherical D) of normal_block_calculations_v2.pdf,
// summarized in 6.1/6.2.
template <typename NoisePolicy>
class NormalBlockKnownClusters : public NormalBlockBase {
  arma::mat C_;      // p x q, fixed cluster-indicator matrix
  arma::mat Gamma_;  // q x q, posterior covariance of W | Y
  arma::mat Mu_;     // n x q, posterior mean of W | Y

  void E_step() override {
    // Gamma = (Omegaq + C^T diag(dm1) C)^{-1}, Mu = R diag(dm1) C Gamma
    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    Gamma_ = arma::inv_sympd(Omegaq_ + arma::diagmat(C_.t() * dm1_));
    arma::mat R = data_.Y - data_.X * B_;
    Mu_ = R * dm1C * Gamma_;
  }

  void M_step() override {
    arma::mat YmMuCT = data_.Y - Mu_ * C_.t();
    B_ = data_.XtXm1 * (data_.X.t() * YmMuCT);
    arma::mat resid = YmMuCT - data_.X * B_;
    arma::vec ddiag = arma::vectorise(arma::mean(arma::square(resid), 0));
    ddiag += C_ * Gamma_.diag();
    dm1_ = NoisePolicy::update_dm1(ddiag);
    arma::mat Sigma_hat = Mu_.t() * Mu_ / data_.n + Gamma_;
    Omegaq_ = estimate_omega(Sigma_hat);
  }

public:
  NormalBlockKnownClusters(const NormalBlockData& data, const arma::mat& C,
                            const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                            double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockBase(data, C.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    C_(C), Gamma_(arma::eye(C.n_cols, C.n_cols)), Mu_(arma::zeros(data.n, C.n_cols)) {}

  // Log-likelihood criterion, see "Criterion" in section 6.1/6.2 of
  // normal_block_calculations_v2.pdf. When sparsity_ > 0, Omegaq_ is no
  // longer the exact inverse of Sigma_hat, so the trace/penalty correction
  // that otherwise cancels out must be added back (mirrors the `if
  // (self$sparsity > 0)` branch of compute_loglik in R/NormalBlockKnownClusters.R).
  double objective() const override {
    double log_det_Omegaq = arma::log_det_sympd(Omegaq_);
    double log_det_Gamma  = arma::log_det_sympd(Gamma_);
    double sum_log_dm1    = arma::sum(arma::log(dm1_));
    double log2pie        = 1.0 + std::log(2.0 * arma::datum::pi);

    double J = -0.5 * (data_.n * data_.p * log2pie - data_.n * sum_log_dm1);
    J += 0.5 * data_.n * (log_det_Omegaq + log_det_Gamma);

    if (sparsity_ > 0.0) {
      arma::mat Sigma_hat = Gamma_ + Mu_.t() * Mu_ / data_.n;
      J += 0.5 * data_.n * q_ - 0.5 * data_.n * arma::trace(Omegaq_ * Sigma_hat);
      J -= sparsity_ * arma::accu(arma::abs(sparsity_weights_ % Omegaq_));
    }
    return J;
  }

  const arma::mat& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H
