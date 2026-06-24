#ifndef NORMALBLOCKR_NORMAL_BLOCK_BASE_H
#define NORMALBLOCKR_NORMAL_BLOCK_BASE_H

#include <RcppArmadillo.h>
#include <vector>
#include <cmath>
#include "normal_block_data.h"
#include "omega_estimation.h"

// Abstract base class, equivalent of the R6 class NB (R/NB-Class.R) stripped
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
  arma::mat sparsity_weights_; // q x q, per-pair penalty weights (see R/NB-Class.R)

  virtual void E_step() = 0;
  virtual void M_step() = 0;

  // Equivalent of the shared private method `NB$get_Omegaq` (R/NB-Class.R):
  // plain inversion when sparsity_ <= 0, graphical lasso (callback to R's
  // glassoFast) otherwise. See omega_estimation.h.
  arma::mat estimate_omega(const arma::mat& Sigma_hat) const {
    return nb_omega::estimate(Sigma_hat, sparsity_, sparsity_weights_);
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

  // Alternates E_step()/M_step(), recording objective() after each full
  // iteration, and stops once the increment falls below `tol` (or after
  // `maxit` iterations). Equivalent of `private$EM_optimize` in R/NB-Class.R.
  void run_em(int maxit, double tol) {
    objective_trace_.clear();
    objective_trace_.push_back(objective());
    niter_ = 0;
    for (int h = 0; h < maxit; ++h) {
      E_step();
      M_step();
      objective_trace_.push_back(objective());
      niter_ = h + 1;
      if (std::abs(objective_trace_[h + 1] - objective_trace_[h]) < tol) break;
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
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_BASE_H
