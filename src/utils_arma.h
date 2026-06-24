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

// adds machine's 0 to elements equal to 0 in x (equivalent of `check_zero_boundary`)
inline arma::mat check_zero_boundary(arma::mat x, double zero = kEps) {
  x.transform([zero](double v) { return (std::isnan(v) || v < zero) ? zero : v; });
  return x;
}

// removes machine's 0 to elements equal to 1 in x (equivalent of `check_one_boundary`)
inline arma::mat check_one_boundary(arma::mat x, double zero = kEps) {
  x.transform([zero](double v) { return std::isnan(v) ? zero : std::min(v, 1.0 - zero); });
  return x;
}

// equivalent of `check_zero_boundary(check_one_boundary(x))`, used to keep
// variational probabilities (tau) away from the {0, 1} boundaries.
inline arma::mat clip_probabilities(const arma::mat& x, double zero = kEps) {
  return check_zero_boundary(check_one_boundary(x, zero), zero);
}

} // namespace nb_utils

#endif // NORMALBLOCKR_UTILS_ARMA_H
