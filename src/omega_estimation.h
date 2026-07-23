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

inline arma::mat estimate(const arma::mat& Sigma_hat, double sparsity, const arma::mat& sparsity_weights) {
  if (sparsity <= 0.0) {
    return arma::inv_sympd(Sigma_hat);
  }

  static Rcpp::Function glassoFast(Rcpp::Environment::namespace_env("glassoFast")["glassoFast"]);

  Rcpp::List glasso_out = glassoFast(Rcpp::wrap(Sigma_hat), Rcpp::wrap(sparsity * sparsity_weights));
  arma::mat Wi = Rcpp::as<arma::mat>(glasso_out["wi"]);

  if (Wi.has_nan()) {
    Rcpp::warning("GLasso fails, the penalty is probably too small and the system badly "
                  "conditionned (reciprocal condition number = %f). "
                  "Sending back the original matrix and its inverse (unpenalized).",
                  arma::rcond(Sigma_hat));
    return arma::inv_sympd(Sigma_hat);
  }
  return 0.5 * (Wi + Wi.t()); // equivalent of Matrix::symmpart(glasso_out$wi)
}

} // namespace nb_omega

#endif // NORMALBLOCKR_OMEGA_ESTIMATION_H
