#ifndef NORMALBLOCKR_NORMAL_BLOCK_EM_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_EM_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <memory>

// (V)EM driver shared by both model families (normal_block_var_base.h and
// normal_block_mean_base.h): the iteration loop, the objective trace, and the
// SQUAREM extrapolation (Varadhan & Roland, 2008) layered on top of plain EM
// to accelerate the slow-linear-EM regime.
//
// Everything family-specific is left to the derived classes: the order of the
// E/M pair (em_cycle()), which parameters SQUAREM extrapolates over
// (get_state()/set_state()) and what makes a candidate infeasible
// (state_is_feasible()). See inst/normal_block_models.qmd for the derivation
// and the acceptance test.
class NormalBlockEMBase {
public:
  virtual ~NormalBlockEMBase() = default;

  // (V)EM criterion (log-likelihood or ELBO) at the current parameter values.
  virtual double objective() const = 0;

  // Deep copy/restore of the entire concrete object, used to snapshot a
  // known-good iterate before an extrapolation and back out of it if
  // infeasible. clone() relies on the compiler-generated copy constructor;
  // restore_from() copies field-by-field since the data reference member
  // blocks copy assignment.
  virtual std::unique_ptr<NormalBlockEMBase> clone() const = 0;
  virtual void restore_from(const NormalBlockEMBase& other) = 0;

  void run_em(int maxit, double tol) {
    objective_trace_.clear();
    objective_trace_.push_back(objective());
    niter_ = 0;
    if (!supports_acceleration()) {
      for (int h = 0; h < maxit; ++h) {
        em_cycle();
        objective_trace_.push_back(objective());
        niter_ = h + 1;
        if (std::abs(objective_trace_[h + 1] - objective_trace_[h]) < tol) break;
      }
      return;
    }

    while (niter_ < maxit) {
      arma::vec p0 = get_state();

      em_cycle();
      double obj1 = objective();
      objective_trace_.push_back(obj1);
      ++niter_;
      if (std::abs(obj1 - objective_trace_[objective_trace_.size() - 2]) < tol) return;
      if (niter_ >= maxit) return;
      arma::vec p1 = get_state();

      em_cycle();
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

  const std::vector<double>& objective_trace() const { return objective_trace_; }
  int niter() const { return niter_; }

protected:
  // One E/M pair, in whichever order the family's recursion uses.
  virtual void em_cycle() = 0;

  // Parameters SQUAREM extrapolates over, and the cheap guard rejecting
  // numerically degenerate candidates outright (the objective() comparison in
  // try_squarem_step() is what actually judges quality).
  virtual arma::vec get_state() const = 0;
  virtual void set_state(const arma::vec& p) = 0;
  virtual bool state_is_feasible(const arma::vec& p) const = 0;

  // Opt-in gate for the acceleration; off by default.
  virtual bool supports_acceleration() const { return false; }

private:
  std::vector<double> objective_trace_;
  int niter_ = 0;

  // One SQUAREM cycle from the current (p2, objective() == obj2) state;
  // backtracks the steplength (Varadhan & Roland's scheme) until the
  // stabilized candidate is at least as good as obj2, or gives up and
  // restores p2.
  bool try_squarem_step(const arma::vec& p0, const arma::vec& p1, const arma::vec& p2, double obj2) {
    arma::vec r = p1 - p0;
    arma::vec v = (p2 - p1) - r;
    double vv = arma::dot(v, v);
    if (vv < 1e-14) return false; // already essentially in the linear/converged regime

    double alpha = -std::sqrt(arma::dot(r, r) / vv); // steplength "S3" of Varadhan & Roland (2008)
    double accept_slack = 1e-8 * (1.0 + std::abs(obj2)); // tolerate floating-point-scale noise only

    std::unique_ptr<NormalBlockEMBase> snapshot = clone(); // == current (p2) state, to back out to
    for (int attempt = 0; attempt < 10 && alpha < -1.0; ++attempt) {
      arma::vec p_new = p0 - 2.0 * alpha * r + (alpha * alpha) * v;
      if (state_is_feasible(p_new)) {
        set_state(p_new);
        em_cycle();
        if (state_is_feasible(get_state()) && objective() >= obj2 - accept_slack) return true;
        restore_from(*snapshot);
      }
      alpha = 0.5 * (alpha - 1.0); // halve the distance to alpha = -1 (no extrapolation)
    }
    return false;
  }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_EM_BASE_H
