// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include "normal_block_data.h"
#include "normal_block_var_types.h"
#include "zi_normal_block_data.h"
#include "zi_normal_block_var_types.h"

namespace {

// arma::vec wraps to a NumericVector with a (n, 1) "dim" attribute; strip it
// so that R-side fields look like plain vectors (consistent with their R6
// counterparts, e.g. private$dm1).
Rcpp::NumericVector to_rvector(const arma::vec& v) {
  return Rcpp::NumericVector(v.begin(), v.end());
}

template <typename Model>
Rcpp::List known_clusters_result(const Model& model) {
  return Rcpp::List::create(
    Rcpp::Named("B")         = model.B(),
    Rcpp::Named("dm1")       = to_rvector(model.dm1()),
    Rcpp::Named("Omega")    = model.Omega(),
    Rcpp::Named("gamma")     = model.Gamma(),
    Rcpp::Named("mu")        = model.Mu(),
    Rcpp::Named("objective") = Rcpp::wrap(model.objective_trace()),
    Rcpp::Named("niter")     = model.niter()
  );
}

template <typename Model>
Rcpp::List unknown_clusters_result(const Model& model) {
  return Rcpp::List::create(
    Rcpp::Named("B")         = model.B(),
    Rcpp::Named("dm1")       = to_rvector(model.dm1()),
    Rcpp::Named("Omega")    = model.Omega(),
    Rcpp::Named("C")         = model.C(),
    Rcpp::Named("alpha")     = to_rvector(model.alpha()),
    Rcpp::Named("M")         = model.M(),
    Rcpp::Named("S")         = to_rvector(model.S()),
    Rcpp::Named("objective") = Rcpp::wrap(model.objective_trace()),
    Rcpp::Named("niter")     = model.niter()
  );
}

template <typename Model>
Rcpp::List ZINormalBlockVarKnownClusters_result(const Model& model) {
  return Rcpp::List::create(
    Rcpp::Named("B")         = model.B(),
    Rcpp::Named("dm1")       = to_rvector(model.dm1()),
    Rcpp::Named("Omega")    = model.Omega(),
    Rcpp::Named("gamma")     = model.Gamma(), // q x q x n array
    Rcpp::Named("mu")        = model.Mu(),
    Rcpp::Named("objective") = Rcpp::wrap(model.objective_trace()),
    Rcpp::Named("niter")     = model.niter()
  );
}

template <typename Model>
Rcpp::List ZINormalBlockVarUnknownClusters_result(const Model& model) {
  return Rcpp::List::create(
    Rcpp::Named("B")         = model.B(),
    Rcpp::Named("dm1")       = to_rvector(model.dm1()),
    Rcpp::Named("Omega")    = model.Omega(),
    Rcpp::Named("C")         = model.C(),
    Rcpp::Named("alpha")     = to_rvector(model.alpha()),
    Rcpp::Named("M")         = model.M(),
    Rcpp::Named("S")         = model.S(), // n x q matrix
    Rcpp::Named("objective") = Rcpp::wrap(model.objective_trace()),
    Rcpp::Named("niter")     = model.niter()
  );
}

} // namespace

//' Fit a normal-block model with known clusters (Rcpp/Armadillo core)
//'
//' Equivalent of the R6 class NormalBlockVarKnownClusters (R/NormalBlockVarKnownClusters.R),
//' implementing sections 2/3 (summarized in 6.1/6.2) of
//' normal_block_calculations_v2.pdf. Initialization is entirely done on the R
//' side: this function only runs the EM recursion from the parameters given.
//'
//' @param Y response matrix (n x p)
//' @param X design matrix (n x d)
//' @param C fixed cluster-indicator matrix (p x q)
//' @param B0 initial regression coefficients (d x p)
//' @param dm1_0 initial inverse variance per variable (length p)
//' @param Omega0 initial precision matrix of the blocks (q x q)
//' @param sparsity sparsity penalty applied to Omega through the graphical
//' lasso (glassoFast); 0 means an unpenalized inversion
//' @param sparsity_weights q x q matrix of per-pair penalty weights (see
//' R/NormalBlockVarBase.R, `sparsity_weights`); only used when sparsity > 0
//' @param noise_covariance either "diagonal" or "spherical"
//' @param niter maximum number of EM iterations
//' @param threshold convergence threshold on the objective increment
//' @return a list with the fitted parameters (B, dm1, Omega, gamma, mu), the
//' objective (log-likelihood) trace and the number of iterations performed
//' @noRd
// [[Rcpp::export]]
Rcpp::List NormalBlockVarKnownClusters_fit(const arma::mat& Y, const arma::mat& X, const arma::mat& C,
                                  arma::mat B0, arma::vec dm1_0, arma::mat Omegaq0,
                                  double sparsity, arma::mat sparsity_weights,
                                  std::string noise_covariance,
                                  int niter, double threshold) {
  NormalBlockData data(Y, X);

  if (noise_covariance == "diagonal") {
    norm_block_var_cov_diag_noise_known_clusters model(data, C, B0, dm1_0, Omegaq0, sparsity, sparsity_weights);
    model.run_em(niter, threshold);
    return known_clusters_result(model);
  } else if (noise_covariance == "spherical") {
    norm_block_var_cov_spherical_noise_known_clusters model(data, C, B0, dm1_0, Omegaq0, sparsity, sparsity_weights);
    model.run_em(niter, threshold);
    return known_clusters_result(model);
  }
  Rcpp::stop("noise_covariance must be either \"diagonal\" or \"spherical\"");
}

//' Fit a normal-block model with unknown clusters (Rcpp/Armadillo core, VEM)
//'
//' Equivalent of the R6 class NormalBlockVarUnknownClusters (R/NormalBlockVarUnknownClusters.R), implementing
//' sections 4/5 (summarized in 6.3) of normal_block_calculations_v2.pdf.
//' Initialization is entirely done on the R side: this function only runs the
//' VEM recursion from the parameters given.
//'
//' @inheritParams NormalBlockVarKnownClusters_fit
//' @param C0 initial variational membership probabilities (p x q)
//' @param alpha0 initial cluster prior probabilities (length q)
//' @param M0 initial variational mean of the cluster effects (n x q)
//' @param S0 initial variational variance of the cluster effects (length q)
//' @param fixed_tau if TRUE, the variational membership probabilities are not
//' re-estimated (useful for stability selection)
//' @return a list with the fitted parameters (B, dm1, Omega, C, alpha, M, S),
//' the ELBO trace and the number of iterations performed
//' @noRd
// [[Rcpp::export]]
Rcpp::List NormalBlockVarUnknownClusters_fit(const arma::mat& Y, const arma::mat& X,
                                    arma::mat B0, arma::vec dm1_0, arma::mat Omegaq0,
                                    arma::mat C0, arma::vec alpha0, arma::mat M0, arma::vec S0,
                                    double sparsity, arma::mat sparsity_weights,
                                    std::string noise_covariance, bool fixed_tau,
                                    int niter, double threshold) {
  NormalBlockData data(Y, X);

  if (noise_covariance == "diagonal") {
    norm_block_var_cov_diag_noise_unknown_clusters model(data, B0, dm1_0, Omegaq0, C0, alpha0, M0, S0,
                                                      sparsity, sparsity_weights, fixed_tau);
    model.run_em(niter, threshold);
    return unknown_clusters_result(model);
  } else if (noise_covariance == "spherical") {
    norm_block_var_cov_spherical_noise_unknown_clusters model(data, B0, dm1_0, Omegaq0, C0, alpha0, M0, S0,
                                                          sparsity, sparsity_weights, fixed_tau);
    model.run_em(niter, threshold);
    return unknown_clusters_result(model);
  }
  Rcpp::stop("noise_covariance must be either \"diagonal\" or \"spherical\"");
}

//' Fit a zero-inflated normal-block model with known clusters (Rcpp/Armadillo core)
//'
//' Equivalent of the R6 class ZINormalBlockVarKnownClusters (R/ZINormalBlockVarKnownClusters.R).
//' Initialization (including the fixed zero-inflation logistic regression)
//' is entirely done on the R side: this function only runs the EM recursion
//' from the parameters given. The zero-inflation mask makes the B-update's
//' normal equations column-specific, so it is solved one column at a time
//' (see src/zi_closed_form_solvers.h).
//'
//' @param Y response matrix (n x p)
//' @param X design matrix (n x d)
//' @param zeros_bar zero-inflation mask (n x p), 1 where Y != 0, 0 where Y == 0
//' @param zi_cond_mean fixed log-likelihood contribution of the (pre-estimated)
//' zero-inflation component (private$ZI_cond_mean in R/NormalBlockVarBase.R)
//' @param C fixed cluster-indicator matrix (p x q)
//' @param B0 initial regression coefficients (d x p)
//' @param dm1_0 initial inverse variance per variable (length p)
//' @param Omega0 initial precision matrix of the blocks (q x q)
//' @param sparsity sparsity penalty applied to Omega through the graphical
//' lasso (glassoFast); 0 means an unpenalized inversion
//' @param sparsity_weights q x q matrix of per-pair penalty weights
//' @param noise_covariance either "diagonal" or "spherical"
//' @param niter maximum number of EM iterations
//' @param threshold convergence threshold on the objective increment
//' @return a list with the fitted parameters (B, dm1, Omega, gamma, mu -- gamma
//' is a q x q x n array, one posterior covariance matrix per row), the
//' log-likelihood trace and the number of iterations performed
//' @noRd
// [[Rcpp::export]]
Rcpp::List ZINormalBlockVarKnownClusters_fit(const arma::mat& Y, const arma::mat& X,
                                  const arma::mat& zeros_bar, double zi_cond_mean, const arma::mat& C,
                                  arma::mat B0, arma::vec dm1_0, arma::mat Omega0,
                                  double sparsity, arma::mat sparsity_weights,
                                  std::string noise_covariance,
                                  int niter, double threshold) {
  ZINormalBlockData data(Y, X, zeros_bar, zi_cond_mean);

  if (noise_covariance == "diagonal") {
    zi_norm_block_var_cov_diag_noise_known_clusters model(data, C, B0, dm1_0, Omegaq0, sparsity, sparsity_weights);
    model.run_em(niter, threshold);
    return ZINormalBlockVarKnownClusters_result(model);
  } else if (noise_covariance == "spherical") {
    zi_norm_block_var_cov_spherical_noise_known_clusters model(data, C, B0, dm1_0, Omegaq0, sparsity, sparsity_weights);
    model.run_em(niter, threshold);
    return ZINormalBlockVarKnownClusters_result(model);
  }
  Rcpp::stop("noise_covariance must be either \"diagonal\" or \"spherical\"");
}

//' Fit a zero-inflated normal-block model with unknown clusters (Rcpp/Armadillo core, VEM)
//'
//' Equivalent of the R6 class ZINormalBlockVarUnknownClusters (R/ZINormalBlockVarUnknownClusters.R).
//' Initialization is entirely done on the R side: this function only runs
//' the VEM recursion from the parameters given. The zero-inflation mask
//' makes the B- and M-updates' normal equations column-/row-specific, so
//' both are solved through direct linear systems rather than an iterative
//' optimizer (see src/zi_closed_form_solvers.h).
//'
//' @inheritParams ZINormalBlockVarKnownClusters_fit
//' @param C0 initial variational membership probabilities (p x q)
//' @param alpha0 initial cluster prior probabilities (length q)
//' @param M0 initial variational mean of the cluster effects (n x q)
//' @param S0 initial variational variance of the cluster effects (n x q,
//' row-dependent because of the zero-inflation mask)
//' @param fixed_tau if TRUE, the variational membership probabilities are not
//' re-estimated (useful for stability selection)
//' @return a list with the fitted parameters (B, dm1, Omega, C, alpha, M, S),
//' the ELBO trace and the number of iterations performed
//' @noRd
// [[Rcpp::export]]
Rcpp::List ZINormalBlockVarUnknownClusters_fit(const arma::mat& Y, const arma::mat& X,
                                    const arma::mat& zeros_bar, double zi_cond_mean,
                                    arma::mat B0, arma::vec dm1_0, arma::mat Omega0,
                                    arma::mat C0, arma::vec alpha0, arma::mat M0, arma::mat S0,
                                    double sparsity, arma::mat sparsity_weights,
                                    std::string noise_covariance, bool fixed_tau,
                                    int niter, double threshold) {
  ZINormalBlockData data(Y, X, zeros_bar, zi_cond_mean);

  if (noise_covariance == "diagonal") {
    zi_norm_block_var_cov_diag_noise_unknown_clusters model(data, B0, dm1_0, Omegaq0, C0, alpha0, M0, S0,
                                                         sparsity, sparsity_weights, fixed_tau);
    model.run_em(niter, threshold);
    return ZINormalBlockVarUnknownClusters_result(model);
  } else if (noise_covariance == "spherical") {
    zi_norm_block_var_cov_spherical_noise_unknown_clusters model(data, B0, dm1_0, Omegaq0, C0, alpha0, M0, S0,
                                                              sparsity, sparsity_weights, fixed_tau);
    model.run_em(niter, threshold);
    return ZINormalBlockVarUnknownClusters_result(model);
  }
  Rcpp::stop("noise_covariance must be either \"diagonal\" or \"spherical\"");
}
