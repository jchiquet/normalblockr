#ifndef NORMALBLOCKR_NORMAL_BLOCK_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <memory>
#include "normal_block_data.h"
#include "omega_estimation.h"

// Abstract base class, equivalent of the R6 class NormalBlockBase (R/NormalBlockBase.R) stripped
// of everything related to initialization, heuristics and clustering
// approximation (all of that is done in R and the resulting parameters are
// passed in at construction time): only the data reference, the
// regression/precision parameters shared by every variant, and the generic
// (V)EM loop are kept here. Concrete subclasses implement E_step(), M_step()
// and objective().
class NormalBlockBase {
protected:
  const NormalBlockData& data_;
  int q_;
  arma::mat B_;       // d x p, regression coefficients
  arma::vec dm1_;     // p, inverse variance per variable (1 / diag(D))
  arma::mat Omegaq_;  // q x q, precision matrix of the blocks
  double sparsity_;          // sparsity penalty applied to Omegaq (0 = unpenalized)
  arma::mat sparsity_weights_; // q x q, per-pair penalty weights (see R/NormalBlockBase.R)

  virtual void E_step() = 0;
  virtual void M_step() = 0;

  // Equivalent of the shared private method `NormalBlockBase$get_Omegaq` (R/NormalBlockBase.R):
  // plain inversion when sparsity_ <= 0, graphical lasso (callback to R's
  // glassoFast) otherwise. See omega_estimation.h.
  arma::mat estimate_omega(const arma::mat& Sigma_hat) const {
    return nb_omega::estimate(Sigma_hat, sparsity_, sparsity_weights_);
  }

  // Opt-in gate for the SQUAREM-style acceleration in run_em() (see its
  // docstring). Defaults to off. Currently overridden to (conditionally)
  // `true` by NormalBlockKnownClusters and ZINormalBlockKnownClusters only
  // (both: sparsity_ <= 0; see each class's own override for why sparsity_ >
  // 0 is excluded). Both rely on objective() being the *general* (valid at
  // any theta) marginal log-likelihood, recomputing Mu_/Gamma_ fresh from
  // (B_, dm1_, Omegaq_) on every call rather than reading anything E_step()
  // cached -- without that, the objective-comparison acceptance test in
  // try_squarem_step() would be comparing apples to oranges (theta_old's
  // posterior plugged into theta_new's M-step output is not log p(Y;
  // theta_new), only a lower bound for it). NormalBlockUnknownClusters
  // additionally carries cluster-membership probabilities (tau) and prior
  // probabilities (alpha) that the extrapolated vector below doesn't track,
  // which go stale relative to an aggressively extrapolated theta even if
  // its own objective() were made general the same way -- not attempted.
  // ZINormalBlockUnknownClusters inherits that same problem. Override to
  // return `true` only once a subclass has been validated the same way (see
  // git history for both currently-accelerated classes' validation).
  virtual bool supports_acceleration() const { return false; }

  // Used by every subclass's restore_from() override to copy back the
  // B_/dm1_/Omegaq_ slice of `other` (see restore_from()'s docstring on
  // clone() for why this can't just be `*this = other`).
  void copy_tracked_state_from(const NormalBlockBase& other) {
    B_ = other.B_;
    dm1_ = other.dm1_;
    Omegaq_ = other.Omegaq_;
  }

public:
  NormalBlockBase(const NormalBlockData& data, int q,
                  const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                  double sparsity, const arma::mat& sparsity_weights) :
    data_(data), q_(q), B_(B0), dm1_(dm1_0), Omegaq_(Omegaq0),
    sparsity_(sparsity), sparsity_weights_(sparsity_weights) {}

  virtual ~NormalBlockBase() = default;

  // (V)EM criterion (log-likelihood or ELBO) at the current parameter values.
  virtual double objective() const = 0;

  // Deep copy / restore of the *entire* concrete object (every subclass's
  // extra state -- Gamma_, Mu_, C_, alpha_, M_, S_... -- not just the
  // B_/dm1_/Omegaq_ tracked by SQUAREM below). Used by run_em() to snapshot
  // a known-good iterate before trying a SQUAREM extrapolation, and to back
  // out of it exactly if it turns out infeasible. `clone()` is implemented
  // per concrete subclass as `return std::make_unique<ThisClass>(*this)`,
  // relying on the compiler-generated copy constructor (every member is a
  // plain value type, so a member-wise copy is exactly right). `restore_from()`
  // cannot likewise use the compiler-generated copy *assignment* operator --
  // `data_` is a reference member, so that operator is implicitly deleted --
  // hence the explicit field-by-field copy via `copy_tracked_state_from()`
  // (this class's own B_/dm1_/Omegaq_) plus each subclass's own extra fields.
  virtual std::unique_ptr<NormalBlockBase> clone() const = 0;
  virtual void restore_from(const NormalBlockBase& other) = 0;

  // Alternates E_step()/M_step(), recording objective() after each full
  // iteration, and stops once the increment falls below `tol` (or after
  // `maxit` iterations). Equivalent of `private$EM_optimize` in R/NormalBlockBase.R.
  //
  // Every cycle also attempts a SQUAREM extrapolation (Varadhan & Roland,
  // 2008) on top of plain EM: from three consecutive plain iterates
  // p0 -> p1 -> p2 (p = vec(B_, dm1_, Omegaq_); everything else is a pure
  // byproduct of E_step()/M_step() given those three, exactly as in plain
  // EM), it extrapolates a point much further along the p0->p1->p2
  // direction than two more EM steps would reach, then "stabilizes" it with
  // one more E_step()/M_step() pass. This is the fix for the slow-linear-EM
  // regime documented in inst/CSDA_analyses (real fits stuck at the `niter`
  // cap as the number of blocks grows): it reaches the same fixed point in
  // far fewer recorded iterations.
  //
  // Acceptance is the textbook SQUAREM test: compare objective() at the
  // stabilized candidate against objective() at p2, and back out (restoring
  // p2 exactly) if it isn't at least as good, backtracking the steplength
  // (Varadhan & Roland's own scheme: alpha <- (alpha - 1) / 2, halving the
  // distance to "no extrapolation") before giving up on this cycle. This
  // requires objective() to be a valid criterion value at *any* theta, not
  // just one reached by a real M-step from its own preceding E-step -- true
  // for NormalBlockKnownClusters (see its objective(), and the comment on
  // supports_acceleration() below for which subclasses this does *not* hold
  // for yet). state_is_feasible() is now only a cheap guard against
  // numerically degenerate candidates (dm1 <= 0, Omegaq not meaningfully
  // PD) that would make E_step()/M_step() produce garbage outright -- the
  // objective comparison is what actually judges quality.
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

  // p = vec(B_, dm1_, Omegaq_): the parameters SQUAREM extrapolates over.
  arma::vec get_state() const {
    return arma::join_vert(arma::join_vert(arma::vectorise(B_), dm1_), arma::vectorise(Omegaq_));
  }

  void set_state(const arma::vec& p) {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    B_      = arma::reshape(p.subvec(0, nB - 1), B_.n_rows, B_.n_cols);
    dm1_    = p.subvec(nB, nB + nd - 1);
    Omegaq_ = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omegaq_.n_rows, Omegaq_.n_cols);
  }

  // Cheap guard against numerically degenerate candidates -- not a quality
  // judgment (objective() comparison in try_squarem_step() is what actually
  // decides that now). dm1 must stay positive (it is an inverse variance)
  // and Omegaq must stay symmetric positive-definite with *some* margin:
  // a technically-PD but near-machine-singular candidate can make
  // log_det_sympd()/inv_sympd() return outright garbage (not just lose a
  // few digits) rather than erroring, which would corrupt the very
  // objective() comparison meant to catch a bad step. 1e-8 is far looser
  // than what an earlier, purely feasibility-gated version of this safety
  // net needed (it had no other way to judge quality, so it leaned on tight
  // conditioning as a proxy for it); here it only needs to rule out
  // candidates that are *unusable*, not merely poor.
  bool state_is_feasible(const arma::vec& p) const {
    arma::uword nB = B_.n_elem, nd = dm1_.n_elem;
    if (arma::any(p.subvec(nB, nB + nd - 1) <= 0.0)) return false;
    arma::mat Omega_cand = arma::reshape(p.subvec(nB + nd, p.n_elem - 1), Omegaq_.n_rows, Omegaq_.n_cols);
    arma::vec eigval;
    if (!arma::eig_sym(eigval, arma::symmatu(Omega_cand))) return false;
    return eigval.min() > 1e-8 * eigval.max();
  }

  // One SQUAREM cycle from the current (valid, == p2, objective() == obj2)
  // state. Varadhan & Roland's own backtracking rule: starting from the
  // data-driven steplength "S3", repeatedly halve the *distance to no
  // extrapolation* (alpha <- (alpha - 1) / 2) until the (feasible,
  // stabilized) candidate's objective() is at least as good as obj2, or
  // until alpha reaches -1 (no extrapolation at all -- not worth taking,
  // give up and keep plain EM's p2). Every attempt is judged on actual
  // objective() value, not a feasibility/conditioning proxy, since
  // objective() is valid at any theta for this subclass (see run_em()'s
  // docstring). On success, leaves the model in the accelerated+stabilized
  // state and returns true; on failure, leaves it exactly back at p2
  // (restored from a snapshot) and returns false.
  bool try_squarem_step(const arma::vec& p0, const arma::vec& p1, const arma::vec& p2, double obj2) {
    arma::vec r = p1 - p0;
    arma::vec v = (p2 - p1) - r;
    double vv = arma::dot(v, v);
    if (vv < 1e-14) return false; // already essentially in the linear/converged regime

    double alpha = -std::sqrt(arma::dot(r, r) / vv); // steplength "S3" of Varadhan & Roland (2008)
    double accept_slack = 1e-8 * (1.0 + std::abs(obj2)); // tolerate floating-point-scale noise only

    std::unique_ptr<NormalBlockBase> snapshot = clone(); // == current (p2) state, to back out to
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

#endif // NORMALBLOCKR_NORMAL_BLOCK_BASE_H
