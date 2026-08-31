#ifndef NORMALBLOCKR_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_mean_base.h"
#include "normal_block_data.h"

// Mean-block model with a known clustering: plain EM alternating B and
// Sigma/Omega, the clustering C staying fixed. Equivalent of the R6 class
// NormalBlockMeanKnownClusters (R/NormalBlockMeanKnownClusters.R).
class NormalBlockMeanKnownClusters : public NormalBlockMeanBase {
  arma::mat C_; // p x q, fixed cluster indicator

  // B = (X'X)^-1 X'Y Omega C (C' Omega C)^-1, then Sigma from the residuals.
  void M_step() override {
    arma::mat CtOC = (C_.t() * Omega_) * C_;
    B_ = (((data_.XtXm1 * data_.XtY) * Omega_) * C_) * arma::inv(CtOC);
    arma::mat R = data_.Y - (data_.X * B_) * C_.t();
    Omega_ = estimate_omega(R.t() * R / data_.n);
  }

  void E_step() override {}

public:
  NormalBlockMeanKnownClusters(const NormalBlockData& data, const arma::mat& C,
                               const arma::mat& B0, const arma::mat& Omega0,
                               double sparsity, const arma::mat& sparsity_weights,
                               bool accelerate) :
    NormalBlockMeanBase(data, C.n_cols, B0, Omega0, sparsity, sparsity_weights, accelerate),
    C_(C) {}

  std::unique_ptr<NormalBlockMeanBase> clone() const override {
    return std::make_unique<NormalBlockMeanKnownClusters>(*this);
  }
  void restore_from(const NormalBlockMeanBase& other) override {
    copy_tracked_state_from(other);
  }

  double objective() const override {
    arma::mat R = data_.Y - (data_.X * B_) * C_.t();
    return -0.5 * data_.n * data_.p * std::log(2.0 * arma::datum::pi)
           + 0.5 * data_.n * arma::log_det_sympd(Omega_)
           - 0.5 * arma::accu(R % (R * Omega_));
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H
