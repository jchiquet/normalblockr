#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "zi_normal_block_mean_base.h"

// Zero-inflated mean-block model with a known clustering: plain EM alternating
// B and Sigma, C staying fixed. Equivalent of the R6 class
// ZINormalBlockMeanKnownClusters (R/ZINormalBlockMeanKnownClusters.R).
class ZINormalBlockMeanKnownClusters : public ZINormalBlockMeanBase {
  arma::mat C_; // p x q, fixed cluster indicator

  void M_step() override {
    update_B(C_);
    refresh_stats();
    update_omega(s_from(C_));
  }

  void E_step() override {}

public:
  ZINormalBlockMeanKnownClusters(const ZINormalBlockData& data, const arma::mat& C,
                                 const arma::mat& B0, const arma::mat& Omega0,
                                 bool accelerate, const std::string& cov_structure) :
    ZINormalBlockMeanBase(data, C.n_cols, B0, Omega0, accelerate, cov_structure),
    C_(C) {}

  double objective() const override {
    refresh_stats();
    return gaussian_objective(C_);
  }

  std::unique_ptr<NormalBlockEMBase> clone() const override {
    return std::make_unique<ZINormalBlockMeanKnownClusters>(*this);
  }
  void restore_from(const NormalBlockEMBase& other) override {
    copy_tracked_state_from(static_cast<const ZINormalBlockMeanKnownClusters&>(other));
  }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_KNOWN_CLUSTERS_H
