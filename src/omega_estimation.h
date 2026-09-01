#ifndef NORMALBLOCKR_OMEGA_ESTIMATION_H
#define NORMALBLOCKR_OMEGA_ESTIMATION_H

#include <RcppArmadillo.h>

// Equivalent of the shared private method `NormalBlocVarkBase$get_Omega` (R/NormalBlockVarBase.R),
// used identically by the M-step of both NormalBlockVarKnownClusters and NormalBlockVarUnknownClusters:
// estimate the precision matrix of the blocks from its covariance estimate
// Sigma_hat (q x q). When `sparsity <= 0`, this is a plain matrix inversion;
// otherwise it calls back into R's glassoFast::glassoFast() (graphical
// lasso). There is no Armadillo reimplementation of glasso in this package,
// and since q (the number of clusters) is small, the R round-trip cost is
// negligible compared to the O(n x p) work done in the rest of the (V)EM
// step.
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

  // Looked up on every call, deliberately: caching it in a `static` keeps an
  // R object alive outside R's protection discipline, across garbage
  // collections and namespace reloads. A stale SEXP then corrupts memory
  // non-deterministically -- heap corruption or a segfault, depending on the
  // run. The lookup is negligible next to the graphical lasso itself.
  Rcpp::Function glassoFast(Rcpp::Environment::namespace_env("glassoFast")["glassoFast"]);

  Rcpp::List glasso_out = glassoFast(Rcpp::wrap(Sigma_hat), Rcpp::wrap(sparsity * sparsity_weights));
  arma::mat Wi = Rcpp::as<arma::mat>(glasso_out["wi"]);

  if (Wi.has_nan()) {
    Rcpp::warning("GLasso fails, the penalty is probably too small and the system badly "
                  "conditionned (reciprocal condition number = %f). "
                  "Sending back the original matrix and its inverse (unpenalized).",
                  arma::rcond(Sigma_hat));
    return arma::inv_sympd(Sigma_hat);
  }
  return ensure_pd(Wi); // symmpart(), plus a guard against a non-PD glasso output
}

} // namespace nb_omega

#endif // NORMALBLOCKR_OMEGA_ESTIMATION_H
