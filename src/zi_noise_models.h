#ifndef NORMALBLOCKR_ZI_NOISE_MODELS_H
#define NORMALBLOCKR_ZI_NOISE_MODELS_H

#include <RcppArmadillo.h>

// Zero-inflated counterparts of DiagonalNoise/SphericalNoise (noise_models.h).
// `weighted_ssq` is the per-variable sum of squared, zero-inflation-masked
// residuals (`colSums(zeros_bar * A)` in R/zi_known_clusters-Class.R /
// R/zi_unknown_clusters-Class.R); the denominator is the per-variable (resp. total)
// count of non-zero observations (nY/npY) rather than the sample size n,
// since dm1 is estimated only from the non-zero-inflated residuals. Both
// arguments are passed to both policies so that call sites can stay
// noise-policy-agnostic.
struct ZIDiagonalNoise {
  static arma::vec update_dm1(const arma::vec& weighted_ssq, const arma::vec& nY, double /*npY*/) {
    return nY / weighted_ssq;
  }
};

struct ZISphericalNoise {
  static arma::vec update_dm1(const arma::vec& weighted_ssq, const arma::vec& /*nY*/, double npY) {
    double dm1_common = npY / arma::accu(weighted_ssq);
    return arma::vec(weighted_ssq.n_elem, arma::fill::value(dm1_common));
  }
};

#endif // NORMALBLOCKR_ZI_NOISE_MODELS_H
