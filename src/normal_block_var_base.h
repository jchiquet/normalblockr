#ifndef NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <memory>
#include "normal_block_data.h"
#include "normal_block_em_base.h"
#include "omega_estimation.h"

// Abstract base class, equivalent of the R6 class NormalBlockVarBase
// (R/NormalBlockVarBase.R) stripped of initialization/heuristics (done in R;
// resulting parameters are passed in at construction). Concrete subclasses
// implement E_step(), M_step() and objective(). The iteration loop and the
// SQUAREM acceleration live in NormalBlockEMBase, shared with the mean-block
// family. See inst/normal_block_models.qmd.
class NormalBlockVarBase : public NormalBlockEMBase {
protected:
  const NormalBlockData& data_;
  int q_;
  arma::mat B_;       // d x p, regression coefficients
  arma::vec dm1_;     // p, inverse variance per variable (1 / diag(D))
  arma::mat Omega_;  // q x q, precision matrix of the blocks
  double sparsity_;          // sparsity penalty applied to Omega (0 = unpenalized)
  arma::mat sparsity_weights_; // q x q, per-pair penalty weights (see R/NormalBlockBase.R)

  virtual void E_step() = 0;
  virtual void M_step() = 0;
  void em_cycle() override { E_step(); M_step(); }

  // Plain inversion when sparsity_ <= 0, graphical lasso otherwise; see
  // omega_estimation.h.
  arma::mat estimate_omega(const arma::mat& Sigma_hat) const {
    return nb_omega::estimate(Sigma_hat, sparsity_, sparsity_weights_);
  }

  // Lazily-cached X*B_ (see inst/normal_block_models.qmd). Subclasses must
  // update B_ through set_B() (never assign B_ directly) for the cache to
  // stay correct.
  void set_B(const arma::mat& new_B) {
    B_ = new_B;
    XB_valid_ = false;
  }
  const arma::mat& XB() const {
    if (!XB_valid_) {
      XB_cache_ = data_.X * B_;
      XB_valid_ = true;
    }
    return XB_cache_;
  }

  // Copies back the B_/dm1_/Omega_ slice of `other`; used by every
  // subclass's restore_from() (can't just be `*this = other`, see clone()).
  void copy_tracked_state_from(const NormalBlockVarBase& other) {
    B_ = other.B_;
    dm1_ = other.dm1_;
    Omega_ = other.Omega_;
    XB_valid_ = false;
  }

public:
  NormalBlockVarBase(const NormalBlockData& data, int q,
                  const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omega0,
                  double sparsity, const arma::mat& sparsity_weights) :
    data_(data), q_(q), B_(B0), dm1_(dm1_0), Omega_(Omega0),
    sparsity_(sparsity), sparsity_weights_(sparsity_weights) {}

  virtual ~NormalBlockVarBase() = default;

  const arma::mat& B() const { return B_; }
  const arma::vec& dm1() const { return dm1_; }
  const arma::mat& Omega() const { return Omega_; }

protected:
  mutable arma::mat XB_cache_;
  mutable bool XB_valid_ = false;

  // p = vec(B_, dm1_, Omega_): the parameters SQUAREM extrapolates over.
  arma::vec get_state() const override {
    return arma::join_vert(arma::join_vert(arma::vectorise(B_), dm1_), arma::vectorise(Omega_));
  }

  void set_state(const arma::vec& p) override {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    B_      = arma::reshape(p.subvec(0, nB - 1), B_.n_rows, B_.n_cols);
    dm1_    = p.subvec(nB, nB + nd - 1);
    Omega_ = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omega_.n_rows, Omega_.n_cols);
    XB_valid_ = false;
  }

  // Cheap guard against numerically degenerate candidates (dm1 <= 0, Omega
  // not meaningfully PD) that would make E_step()/M_step() produce garbage
  // outright; the objective() comparison in try_squarem_step() is what
  // actually judges quality.
  bool state_is_feasible(const arma::vec& p) const override {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    if (arma::any(p.subvec(nB, nB + nd - 1) <= 0.0)) return false;
    arma::mat Omega_cand = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omega_.n_rows, Omega_.n_cols);
    arma::vec eigval;
    if (!arma::eig_sym(eigval, arma::symmatu(Omega_cand))) return false;
    return eigval.min() > 1e-8 * eigval.max();
  }

};

#endif // NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H
