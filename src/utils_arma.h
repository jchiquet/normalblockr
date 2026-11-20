#ifndef NORMALBLOCKR_UTILS_ARMA_H
#define NORMALBLOCKR_UTILS_ARMA_H

#include <RcppArmadillo.h>
#include <limits>

// Small numerical helpers, Armadillo equivalents of R/utils.R, used by the
// variational E-step of the unknown-clusters models.

namespace nb_utils {

constexpr double kEps = std::numeric_limits<double>::epsilon();

// Row-wise softmax (equivalent of `apply(eta, 1, softmax)` then transposed back).
inline arma::mat softmax_rows(const arma::mat& eta) {
  arma::mat out(eta.n_rows, eta.n_cols);
  for (arma::uword i = 0; i < eta.n_rows; ++i) {
    arma::rowvec row = eta.row(i);
    double b = row.max();
    arma::rowvec e = arma::exp(row - b);
    out.row(i) = e / arma::sum(e);
  }
  return out;
}

// x * log(x), set to 0 when x is (numerically) 0 (equivalent of `xlogx`).
inline double xlogx(double x) {
  return (x < kEps) ? 0.0 : x * std::log(x);
}

inline double sum_xlogx(const arma::mat& x) {
  double out = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) out += xlogx(x(i));
  return out;
}

// keeps variational probabilities (tau) away from the {0, 1} boundaries,
// NaN set to `zero` (equivalent of `clip_probabilities` in R/utils.R)
inline arma::mat clip_probabilities(arma::mat x, double zero = kEps) {
  x.transform([zero](double v) {
    if (std::isnan(v)) return zero;
    return std::max(std::min(v, 1.0 - zero), zero);
  });
  return x;
}

// Inverse and log-determinant of a symmetric positive-definite matrix,
// computed from a single Cholesky factorization (A = R^T R) instead of two
// separate ones (as arma::inv_sympd(A) followed by arma::log_det_sympd(A)
// would do): log_det(A) = 2 * sum(log(diag(R))) is a free byproduct of the
// factorization already needed for the inverse.
struct SympdInverse {
  arma::mat inv;
  double log_det;
};

inline SympdInverse inv_and_log_det_sympd(const arma::mat& A) {
  arma::mat R = arma::chol(A);
  arma::mat Rinv = arma::inv(arma::trimatu(R));
  return SympdInverse{Rinv * Rinv.t(), 2.0 * arma::accu(arma::log(R.diag()))};
}

} // namespace nb_utils

#endif // NORMALBLOCKR_UTILS_ARMA_H
