#ifndef NORMALBLOCKR_NORMAL_BLOCK_DATA_H
#define NORMALBLOCKR_NORMAL_BLOCK_DATA_H

#include <RcppArmadillo.h>

// Equivalent of the R6 class NormalBlockData (R/NormalBlockData.R), restricted to the
// fields required by the non zero-inflated models: the response matrix Y,
// the design matrix X and the cross-product (XtXm1) reused at every M-step.
class NormalBlockData {
public:
  arma::mat Y;     // n x p, response matrix
  arma::mat X;     // n x d, design matrix
  int n;            // number of samples
  int p;            // number of variables
  int d;            // number of covariates
  arma::mat XtX;    // d x d, reused by the mean-block M-step (M = B' XtX B)
  arma::mat XtXm1;  // d x d, inverse of X^T X, reused at every M-step
  arma::mat XtY;    // d x p, X^T Y -- mirrors the R6 field of the same name
                    // (R/NormalBlockData.R). Unread by the variance-block
                    // (V)EM step, which always needs X^T(Y - latent
                    // contribution) instead; used by the mean-block one.

  NormalBlockData(const arma::mat& Y_, const arma::mat& X_) :
    Y(Y_), X(X_), n(Y_.n_rows), p(Y_.n_cols), d(X_.n_cols) {
    XtX   = X.t() * X;
    XtXm1 = arma::inv_sympd(XtX);
    XtY   = X.t() * Y;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_DATA_H
