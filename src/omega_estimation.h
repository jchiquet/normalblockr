#ifndef NORMALBLOCKR_OMEGA_ESTIMATION_H
#define NORMALBLOCKR_OMEGA_ESTIMATION_H

#include <RcppArmadillo.h>
#include "graphical_lasso.h"

// Equivalent of the shared private method `NormalBlockVarBase$get_Omega` (R/NormalBlockVarBase.R),
// used identically by the M-step of both NormalBlockVarKnownClusters and NormalBlockVarUnknownClusters:
// estimate the precision matrix of the blocks from its covariance estimate
// Sigma_hat (q x q). When `sparsity <= 0`, this is a plain matrix inversion;
// otherwise it runs the in-package graphical lasso (src/graphical_lasso.h).
// That solver used to be glassoFast, reached by calling back into R from
// inside this loop -- see graphical_lasso.h for why it isn't any more.
namespace nb_omega {

// Armadillo equivalent of ensure_pd() (R/utils.R): the graphical lasso can
// return a precision matrix that is not quite positive definite (an EM
// iterate's Sigma can be badly conditioned, especially for the p x p Sigma of
// the mean-block family), which would then make log_det_sympd() throw. The
// Cholesky attempt makes the common case free; only a failing candidate pays
// for the eigendecomposition.
inline arma::mat ensure_pd(const arma::mat& M, double floor_value = 1e-6) {
  arma::mat sym = arma::symmatu(M), R;
  if (arma::chol(R, sym)) return sym;
  arma::vec eigval;
  arma::mat eigvec;
  if (!arma::eig_sym(eigval, eigvec, sym)) return sym;
  eigval = arma::clamp(eigval, floor_value, eigval.max());
  return eigvec * arma::diagmat(eigval) * eigvec.t();
}

inline arma::mat estimate(const arma::mat& Sigma_hat, double sparsity, const arma::mat& sparsity_weights) {
  if (sparsity <= 0.0) {
    return arma::inv_sympd(Sigma_hat);
  }

  nb_glasso::Result glasso_out = nb_glasso::solve(Sigma_hat, sparsity * sparsity_weights);

  if (glasso_out.X.has_nan()) {
    Rcpp::warning("GLasso fails, the penalty is probably too small and the system badly "
                  "conditionned (reciprocal condition number = %f). "
                  "Sending back the original matrix and its inverse (unpenalized).",
                  arma::rcond(Sigma_hat));
    return arma::inv_sympd(Sigma_hat);
  }
  return ensure_pd(glasso_out.X); // symmpart(), plus a guard against a non-PD glasso output
}

} // namespace nb_omega

#endif // NORMALBLOCKR_OMEGA_ESTIMATION_H
