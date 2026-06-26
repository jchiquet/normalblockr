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
  arma::mat XtXm1;  // d x d, inverse of X^T X, reused at every M-step
  arma::mat XtY;    // d x p, X^T Y -- mirrors the R6 field of the same name
                    // (R/NormalBlockData.R), used there for the OLS heuristic
                    // initialization of B; the (V)EM step itself always needs
                    // X^T(Y - latent contribution) instead of the plain X^T Y,
                    // so this field is otherwise unread on the C++ side

  NormalBlockData(const arma::mat& Y_, const arma::mat& X_) :
    Y(Y_), X(X_), n(Y_.n_rows), p(Y_.n_cols), d(X_.n_cols) {
    XtXm1 = arma::inv_sympd(X.t() * X);
    XtY   = X.t() * Y;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_DATA_H
