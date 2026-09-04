#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_BASE_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_BASE_H

#include <RcppArmadillo.h>
#include <string>
#include "zi_normal_block_data.h"
#include "normal_block_em_base.h"

// Abstract base of the zero-inflated mean-block family, counterpart of
// NormalBlockMeanBase (normal_block_mean_base.h) for data with structural
// zeros. The zero-inflation probabilities are fitted once on the R side and
// only their fixed log-likelihood contribution (zi_cond_mean) reaches the
// recursion, so the mask enters the Gaussian part as a 0/1 weight:
//
//   l = -npY/2 log(2pi) + 1/2 sum_j nY_j log(w_j)
//       - 1/2 sum_j w_j s_j + zi_cond_mean,   s_j = sum_ik m_ij T_jk R_ijk^2
//
// which is why Sigma is restricted to its diagonal and spherical shapes here.
// A full Sigma would tie the variables together inside each row, and the mask
// makes each row's set of observed variables different: every row would then
// need its own submatrix inverse, i.e. a missing-data EM rather than a
// reweighting (see R/ZINormalBlockMeanBase.R for the user-facing error).
//
// Everything the recursion needs is carried by two p x q statistics, both
// linear in the mask:
//
//   A_jk = sum_i m_ij Z_ik^2,   P_jk = sum_i m_ij Y_ij Z_ik,   Z = X B
//
// so no p x p matrix is ever formed.
class ZINormalBlockMeanBase : public NormalBlockEMBase {
protected:
  const ZINormalBlockData& data_;
  int q_;
  arma::mat B_;      // d x q, one column per cluster
  arma::mat Omega_;  // p x p diagonal, kept in that shape for symmetry with
                     // the non-ZI family (R-side model_par, split(), nb_param)
  bool accelerate_;
  std::string cov_structure_; // "diagonal" or "spherical"

  arma::mat Ym_;     // n x p, Y with the masked entries zeroed out
  arma::vec ssq_;    // p, sum_i m_ij Y_ij^2

  // refreshed from B_ by refresh_stats(); mutable so that objective() can
  // resynchronize them after a SQUAREM extrapolation has moved B_
  mutable arma::mat A_, P_; // p x q

  virtual void M_step() = 0;
  virtual void E_step() = 0; // no-op for known clusterings
  void em_cycle() override { M_step(); E_step(); }

  void refresh_stats() const {
    arma::mat Z = data_.X * B_;
    A_ = data_.zeros_bar.t() * arma::square(Z);
    P_ = Ym_.t() * Z;
  }

  // s_j = sum_i m_ij E_T[(Y_ij - Z_i,c_j)^2]; the variational variance
  // correction (Lambda in the non-ZI family) is already inside the
  // T-weighted A_ term.
  arma::vec s_from(const arma::mat& T) const {
    return ssq_ - 2.0 * arma::sum(T % P_, 1) + arma::sum(T % A_, 1);
  }

  // One weighted least squares per cluster: the mask makes the weight of
  // observation i for cluster k, a_ik = sum_j m_ij w_j T_jk, depend on i, so
  // the single (X'X)^-1 solve of the non-ZI M-step splits into q d x d ones.
  void update_B(const arma::mat& T) {
    arma::rowvec w = Omega_.diag().t();
    arma::mat Mw = data_.zeros_bar; Mw.each_row() %= w;
    arma::mat Yw = Ym_;             Yw.each_row() %= w;
    arma::mat AW = Mw * T; // n x q
    arma::mat BW = Yw * T; // n x q
    for (int k = 0; k < q_; ++k) {
      arma::mat Xa = data_.X;
      Xa.each_col() %= AW.col(k);
      B_.col(k) = arma::solve(arma::symmatu(data_.X.t() * Xa), data_.X.t() * BW.col(k));
    }
  }

  void update_omega(const arma::vec& s) {
    arma::vec ss = arma::clamp(s, arma::datum::eps, arma::datum::inf);
    arma::vec w = (cov_structure_ == "spherical")
      ? arma::vec(data_.p, arma::fill::value(data_.npY / arma::accu(ss)))
      : data_.nY / ss;
    Omega_ = arma::diagmat(w);
  }

  double gaussian_objective(const arma::mat& T) const {
    arma::vec w = Omega_.diag();
    return -0.5 * data_.npY * std::log(2.0 * arma::datum::pi)
           + 0.5 * arma::accu(data_.nY % arma::log(w))
           - 0.5 * arma::accu(w % s_from(T))
           + data_.zi_cond_mean;
  }

  bool supports_acceleration() const override { return accelerate_; }

  void copy_tracked_state_from(const ZINormalBlockMeanBase& other) {
    B_ = other.B_;
    Omega_ = other.Omega_;
  }

public:
  ZINormalBlockMeanBase(const ZINormalBlockData& data, int q,
                        const arma::mat& B0, const arma::mat& Omega0,
                        bool accelerate, const std::string& cov_structure) :
    data_(data), q_(q), B_(B0), Omega_(Omega0), accelerate_(accelerate),
    cov_structure_(cov_structure) {
    Ym_  = data_.Y % data_.zeros_bar;
    ssq_ = arma::vectorise(arma::sum(Ym_ % data_.Y, 0));
    refresh_stats();
  }

  virtual ~ZINormalBlockMeanBase() = default;

  const arma::mat& B() const { return B_; }
  const arma::mat& Omega() const { return Omega_; }

protected:
  // Only Omega's diagonal is tracked: the rest of the matrix is structurally
  // zero and extrapolating it would just add p^2 - p no-ops.
  arma::vec get_state() const override {
    return arma::join_vert(arma::vectorise(B_), Omega_.diag());
  }

  void set_state(const arma::vec& p) override {
    arma::uword nB = B_.n_elem;
    B_ = arma::reshape(p.subvec(0, nB - 1), B_.n_rows, B_.n_cols);
    Omega_ = arma::diagmat(p.subvec(nB, p.n_elem - 1));
  }

  bool state_is_feasible(const arma::vec& p) const override {
    arma::uword nB = B_.n_elem;
    return p.is_finite() && arma::all(p.subvec(nB, p.n_elem - 1) > 0.0);
  }
};

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_BASE_H
