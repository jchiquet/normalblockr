#ifndef NORMALBLOCKR_NORMAL_BLOCK_DATA_H
#define NORMALBLOCKR_NORMAL_BLOCK_DATA_H

#include <RcppArmadillo.h>

// Equivalent of the R6 class NB_data (R/NB_data-Class.R), restricted to the
// fields required by the non zero-inflated models: the response matrix Y,
// the design matrix X and the cross-products reused at every M-step.
class NormalBlockData {
public:
  arma::mat Y;     // n x p, response matrix
  arma::mat X;     // n x d, design matrix
  int n;            // number of samples
  int p;            // number of variables
  int d;            // number of covariates
  arma::mat XtXm1;  // d x d, inverse of X^T X
  arma::mat XtY;    // d x p, X^T Y

  NormalBlockData(const arma::mat& Y_, const arma::mat& X_) :
    Y(Y_), X(X_), n(Y_.n_rows), p(Y_.n_cols), d(X_.n_cols) {
    XtXm1 = arma::inv_sympd(X.t() * X);
    XtY   = X.t() * Y;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_DATA_H
