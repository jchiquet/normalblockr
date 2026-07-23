#ifndef NORMALBLOCKR_NOISE_MODELS_H
#define NORMALBLOCKR_NOISE_MODELS_H

#include <RcppArmadillo.h>

// Noise policies controlling the shape of the residual covariance D = diag(1/dm1).
// `ddiag` is the per-variable mean squared residual (already including the
// uncertainty correction term, e.g. C * diag(Gamma)); `update_dm1` returns
// the corresponding precision vector dm1 = 1 / variance, of length p.
//
// This is the only place where the diagonal/spherical models differ: the
// (V)EM criterion is written generically in terms of the vector dm1 in both
// cases (see normal_block_var_known_clusters.h / normal_block_var_unknown_clusters.h).

// D diagonal, variable-specific variance: dm1_j = 1 / ddiag_j (R/NormalBlockVarKnownClusters.R, "diagonal" branch)
struct DiagonalNoise {
  static arma::vec update_dm1(const arma::vec& ddiag) {
    return 1.0 / ddiag;
  }
};

// D = xi * I_p, common variance across variables (R/NormalBlockVarKnownClusters.R, "spherical" branch)
struct SphericalNoise {
  static arma::vec update_dm1(const arma::vec& ddiag) {
    double dm1_common = 1.0 / arma::mean(ddiag);
    return arma::vec(ddiag.n_elem, arma::fill::value(dm1_common));
  }
};

#endif // NORMALBLOCKR_NOISE_MODELS_H
