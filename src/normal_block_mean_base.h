#ifndef NORMALBLOCKR_NORMAL_BLOCK_MEAN_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_MEAN_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include <memory>
#include "normal_block_data.h"
#include "omega_estimation.h"

// Abstract base class, equivalent of the R6 class NormalBlockMeanBase
// (R/NormalBlockMeanBase.R) stripped of initialization/heuristics (done in R;
// resulting parameters are passed in at construction).
//
// Mirrors NormalBlockVarBase (normal_block_var_base.h) but for the family
// where the clustering structures the *mean* (mu_i = C B' X_i), which changes
// the state it carries: B_ is d x q (one column per cluster) and Omega_ is
// p x p (variable-level), with no per-variable noise vector. run_em() also
// runs M *then* E, matching R/NormalBlockMeanUnknownClusters.R.
class NormalBlockMeanBase {
protected:
  const NormalBlockData& data_;
  int q_;
  arma::mat B_;      // d x q, regression coefficients (one column per cluster)
  arma::mat Omega_;  // p x p, precision matrix of the variables
  double sparsity_;
  arma::mat sparsity_weights_; // p x p

  virtual void M_step() = 0;
  virtual void E_step() = 0; // no-op for known clusterings

  arma::mat estimate_omega(const arma::mat& Sigma_hat) const {
    return nb_omega::estimate(Sigma_hat, sparsity_, sparsity_weights_);
  }

  // Hook for the SQUAREM acceleration, deliberately not implemented yet: the
  // port is kept a faithful translation of the R recursion so both can be
  // compared exactly (tests/testthat/test-cpp-normal-block-mean.R).
  virtual bool supports_acceleration() const { return false; }

public:
  NormalBlockMeanBase(const NormalBlockData& data, int q,
                      const arma::mat& B0, const arma::mat& Omega0,
                      double sparsity, const arma::mat& sparsity_weights) :
    data_(data), q_(q), B_(B0), Omega_(Omega0),
    sparsity_(sparsity), sparsity_weights_(sparsity_weights) {}

  virtual ~NormalBlockMeanBase() = default;

  // (V)EM criterion at the current parameter values. Unlike the
  // variance-block family, the lasso penalty is *not* subtracted here: the R
  // side adds it back through the `sparsity_term` active binding.
  virtual double objective() const = 0;

  void run_em(int maxit, double tol) {
    objective_trace_.clear();
    objective_trace_.push_back(objective());
    niter_ = 0;
    for (int h = 0; h < maxit; ++h) {
      M_step(); E_step();
      objective_trace_.push_back(objective());
      niter_ = h + 1;
      if (std::abs(objective_trace_[h + 1] - objective_trace_[h]) < tol) break;
    }
  }

  const arma::mat& B() const { return B_; }
  const arma::mat& Omega() const { return Omega_; }
  const std::vector<double>& objective_trace() const { return objective_trace_; }
  int niter() const { return niter_; }

private:
  std::vector<double> objective_trace_;
  int niter_ = 0;
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_MEAN_BASE_H
