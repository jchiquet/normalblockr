#ifndef NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <memory>
#include "normal_block_data.h"
#include "omega_estimation.h"

// Abstract base class, equivalent of the R6 class NormalBlockVarBase
// (R/NormalBlockVarBase.R) stripped of initialization/heuristics (done in R;
// resulting parameters are passed in at construction). Concrete subclasses
// implement E_step(), M_step() and objective(). See inst/normal_block_models.qmd.
class NormalBlockVarBase {
protected:
  const NormalBlockData& data_;
  int q_;
  arma::mat B_;       // d x p, regression coefficients
  arma::vec dm1_;     // p, inverse variance per variable (1 / diag(D))
  arma::mat Omegaq_;  // q x q, precision matrix of the blocks
  double sparsity_;          // sparsity penalty applied to Omegaq (0 = unpenalized)
  arma::mat sparsity_weights_; // q x q, per-pair penalty weights (see R/NormalBlockVarBase.R)

  virtual void E_step() = 0;
  virtual void M_step() = 0;

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

  // Opt-in gate for the SQUAREM-style acceleration in run_em(). Off by
  // default; each leaf class overrides it (typically sparsity_ <= 0 only --
  // see inst/normal_block_models.qmd for why sparsity excludes it).
  virtual bool supports_acceleration() const { return false; }

  // Copies back the B_/dm1_/Omegaq_ slice of `other`; used by every
  // subclass's restore_from() (can't just be `*this = other`, see clone()).
  void copy_tracked_state_from(const NormalBlockVarBase& other) {
    B_ = other.B_;
    dm1_ = other.dm1_;
    Omegaq_ = other.Omegaq_;
    XB_valid_ = false;
  }

public:
  NormalBlockVarBase(const NormalBlockData& data, int q,
                  const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                  double sparsity, const arma::mat& sparsity_weights) :
    data_(data), q_(q), B_(B0), dm1_(dm1_0), Omegaq_(Omegaq0),
    sparsity_(sparsity), sparsity_weights_(sparsity_weights) {}

  virtual ~NormalBlockVarBase() = default;

  // (V)EM criterion (log-likelihood or ELBO) at the current parameter values.
  virtual double objective() const = 0;

  // Deep copy/restore of the entire concrete object; used by run_em() to
  // snapshot a known-good iterate before a SQUAREM extrapolation and back
  // out of it if infeasible. clone() relies on the compiler-generated copy
  // constructor; restore_from() copies field-by-field since data_ (a
  // reference member) blocks copy assignment.
  virtual std::unique_ptr<NormalBlockVarBase> clone() const = 0;
  virtual void restore_from(const NormalBlockVarBase& other) = 0;

  // Alternates E_step()/M_step(), recording objective() each iteration,
  // stopping once the increment falls below `tol` or after `maxit`
  // iterations -- equivalent of `private$EM_optimize` in R/NormalBlockVarBase.R.
  // Each cycle also attempts a SQUAREM extrapolation (Varadhan & Roland,
  // 2008) on top of plain EM to accelerate the slow-linear-EM regime; see
  // inst/normal_block_models.qmd for the full derivation and acceptance test.
  void run_em(int maxit, double tol) {
    objective_trace_.clear();
    objective_trace_.push_back(objective());
    niter_ = 0;
    if (!supports_acceleration()) {
      for (int h = 0; h < maxit; ++h) {
        E_step(); M_step();
        objective_trace_.push_back(objective());
        niter_ = h + 1;
        if (std::abs(objective_trace_[h + 1] - objective_trace_[h]) < tol) break;
      }
      return;
    }

    while (niter_ < maxit) {
      arma::vec p0 = get_state();

      E_step(); M_step();
      double obj1 = objective();
      objective_trace_.push_back(obj1);
      ++niter_;
      if (std::abs(obj1 - objective_trace_[objective_trace_.size() - 2]) < tol) return;
      if (niter_ >= maxit) return;
      arma::vec p1 = get_state();

      E_step(); M_step();
      double obj2 = objective();
      objective_trace_.push_back(obj2);
      ++niter_;
      if (std::abs(obj2 - obj1) < tol) return;
      if (niter_ >= maxit) return;
      arma::vec p2 = get_state();

      if (try_squarem_step(p0, p1, p2, obj2)) {
        double obj_acc = objective();
        objective_trace_.push_back(obj_acc);
        ++niter_;
        if (std::abs(obj_acc - obj2) < tol) return;
      }
    }
  }

  const arma::mat& B() const { return B_; }
  const arma::vec& dm1() const { return dm1_; }
  const arma::mat& Omegaq() const { return Omegaq_; }
  const std::vector<double>& objective_trace() const { return objective_trace_; }
  int niter() const { return niter_; }

private:
  std::vector<double> objective_trace_;
  int niter_ = 0;
  mutable arma::mat XB_cache_;
  mutable bool XB_valid_ = false;

  // p = vec(B_, dm1_, Omegaq_): the parameters SQUAREM extrapolates over.
  arma::vec get_state() const {
    return arma::join_vert(arma::join_vert(arma::vectorise(B_), dm1_), arma::vectorise(Omegaq_));
  }

  void set_state(const arma::vec& p) {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    B_      = arma::reshape(p.subvec(0, nB - 1), B_.n_rows, B_.n_cols);
    dm1_    = p.subvec(nB, nB + nd - 1);
    Omegaq_ = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omegaq_.n_rows, Omegaq_.n_cols);
    XB_valid_ = false;
  }

  // Cheap guard against numerically degenerate candidates (dm1 <= 0, Omegaq
  // not meaningfully PD) that would make E_step()/M_step() produce garbage
  // outright; the objective() comparison in try_squarem_step() is what
  // actually judges quality.
  bool state_is_feasible(const arma::vec& p) const {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    if (arma::any(p.subvec(nB, nB + nd - 1) <= 0.0)) return false;
    arma::mat Omega_cand = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omegaq_.n_rows, Omegaq_.n_cols);
    arma::vec eigval;
    if (!arma::eig_sym(eigval, arma::symmatu(Omega_cand))) return false;
    return eigval.min() > 1e-8 * eigval.max();
  }

  // One SQUAREM cycle from the current (p2, objective() == obj2) state;
  // backtracks the steplength (Varadhan & Roland's scheme) until the
  // stabilized candidate is at least as good as obj2, or gives up and
  // restores p2. See inst/normal_block_models.qmd for the full derivation.
  bool try_squarem_step(const arma::vec& p0, const arma::vec& p1, const arma::vec& p2, double obj2) {
    arma::vec r = p1 - p0;
    arma::vec v = (p2 - p1) - r;
    double vv = arma::dot(v, v);
    if (vv < 1e-14) return false; // already essentially in the linear/converged regime

    double alpha = -std::sqrt(arma::dot(r, r) / vv); // steplength "S3" of Varadhan & Roland (2008)
    double accept_slack = 1e-8 * (1.0 + std::abs(obj2)); // tolerate floating-point-scale noise only

    std::unique_ptr<NormalBlockVarBase> snapshot = clone(); // == current (p2) state, to back out to
    for (int attempt = 0; attempt < 10 && alpha < -1.0; ++attempt) {
      arma::vec p_new = p0 - 2.0 * alpha * r + (alpha * alpha) * v;
      if (state_is_feasible(p_new)) {
        set_state(p_new);
        E_step();
        M_step();
        if (state_is_feasible(get_state()) && objective() >= obj2 - accept_slack) return true;
        restore_from(*snapshot);
      }
      alpha = 0.5 * (alpha - 1.0); // halve the distance to alpha = -1 (no extrapolation)
    }
    return false;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_VAR_BASE_H
