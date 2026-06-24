#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_DATA_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_DATA_H

#include <RcppArmadillo.h>
#include "normal_block_data.h"

// Zero-inflation bookkeeping added to NormalBlockData when Y has structural
// zeros (R/NormalBlockData.R: zeros, zeros_bar, npY, nY), plus the fixed
// log-likelihood contribution of the zero-inflation component, computed once
// in NormalBlockBase$initialize (R/NormalBlockBase.R: private$ZI_cond_mean). The zero-inflation
// probabilities themselves (kappa, estimated upfront by a fixed logistic
// regression on X0) are never revisited by the (V)EM recursion, so only
// their fixed log-likelihood contribution needs to be carried over here.
class ZINormalBlockData : public NormalBlockData {
public:
  arma::mat zeros_bar;  // n x p, 1 where Y != 0, 0 where Y == 0
  arma::vec nY;          // p, number of non-zero observations per variable
  double npY;            // total number of non-zero observations
  double zi_cond_mean;   // fixed contribution of the (pre-estimated) ZI component to the log-likelihood

  ZINormalBlockData(const arma::mat& Y_, const arma::mat& X_,
                     const arma::mat& zeros_bar_, double zi_cond_mean_) :
    NormalBlockData(Y_, X_), zeros_bar(zeros_bar_),
    nY(arma::vectorise(arma::sum(zeros_bar_, 0))),
    npY(arma::accu(zeros_bar_)),
    zi_cond_mean(zi_cond_mean_) {}
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_DATA_H
