#ifndef NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H
#define NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H

#include <RcppArmadillo.h>
#include "normal_block_base.h"
#include "normal_block_data.h"
#include "utils_arma.h"

// Normal-block model with known clusters (fixed indicator matrix C), templated
// on the residual-noise policy (DiagonalNoise / SphericalNoise, see
// noise_models.h). Equivalent of the R6 class NormalBlockKnownClusters
// (R/NormalBlockKnownClusters.R), implementing the EM recursion of section 2
// (diagonal D) / section 3 (spherical D) of normal_block_calculations_v2.pdf,
// summarized in 6.1/6.2.
template <typename NoisePolicy>
class NormalBlockKnownClusters : public NormalBlockBase {
  arma::mat C_;      // p x q, fixed cluster-indicator matrix
  arma::mat Gamma_;  // q x q, posterior covariance of W | Y
  arma::mat Mu_;     // n x q, posterior mean of W | Y

  void E_step() override {
    // Gamma = (Omegaq + C^T diag(dm1) C)^{-1}, Mu = R diag(dm1) C Gamma
    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    Gamma_ = arma::inv_sympd(Omegaq_ + arma::diagmat(C_.t() * dm1_));
    arma::mat R = data_.Y - data_.X * B_;
    Mu_ = R * dm1C * Gamma_;
  }

  void M_step() override {
    arma::mat YmMuCT = data_.Y - Mu_ * C_.t();
    B_ = data_.XtXm1 * (data_.X.t() * YmMuCT);
    arma::mat resid = YmMuCT - data_.X * B_;
    arma::vec ddiag = arma::vectorise(arma::mean(arma::square(resid), 0));
    ddiag += C_ * Gamma_.diag();
    dm1_ = NoisePolicy::update_dm1(ddiag);
    arma::mat Sigma_hat = Mu_.t() * Mu_ / data_.n + Gamma_;
    Omegaq_ = estimate_omega(Sigma_hat);
  }

public:
  NormalBlockKnownClusters(const NormalBlockData& data, const arma::mat& C,
                            const arma::mat& B0, const arma::vec& dm1_0, const arma::mat& Omegaq0,
                            double sparsity, const arma::mat& sparsity_weights) :
    NormalBlockBase(data, C.n_cols, B0, dm1_0, Omegaq0, sparsity, sparsity_weights),
    C_(C), Gamma_(arma::eye(C.n_cols, C.n_cols)), Mu_(arma::zeros(data.n, C.n_cols)) {}

  // General (non-profiled) marginal log-likelihood of Y, valid at *any*
  // (B_, dm1_, Omegaq_), not just at an M-step optimum -- see "Criterion"
  // in section 6.1/6.2 of normal_block_calculations_v2.pdf for the
  // derivation. The model's exact marginal is Y_i ~ N(B^T X_i, Sigma_Y),
  // Sigma_Y = D + C Omegaq^{-1} C^T (D = diag(1/dm1_), p x p); rather than
  // forming that p x p matrix directly (expensive, and ill-conditioned
  // whenever q << p), |Sigma_Y| and Sigma_Y^{-1} are obtained from the
  // matrix determinant lemma and the Woodbury identity, both expressed
  // through the q x q posterior precision Gamma^{-1} = Omegaq + C^T
  // diag(dm1) C (the same matrix E_step() inverts to get Gamma/Mu --
  // recomputed fresh here rather than read from the cached Gamma_/Mu_,
  // which lag one E_step() behind whenever objective() is called right
  // after M_step(), as it always is in run_em()):
  //
  //   log|Sigma_Y|                = -sum(log(dm1)) - log|Omegaq| - log|Gamma|
  //   trace(Sigma_Y^{-1} R^T R)   = sum(dm1 % colSums(R^2)) - trace(Gamma^{-1} Mu^T Mu)
  //
  // where R = Y - XB and Mu = R diag(dm1) C Gamma is the exact posterior
  // mean *for this same theta* -- an identity that holds for any theta, not
  // just an M-step optimum (E_step()'s formula is ordinary Gaussian
  // conditioning, not an approximation, for this known-clusters model).
  //
  // This replaces an earlier "profiled" shortcut (valid only when dm1_ was
  // exactly the M-step's own argmax for the *cached* Mu_/Gamma_, collapsing
  // the quadratic residual term to a constant) that carried a sign bug --
  // it added log|Gamma^{-1}| where the algebra needs log|Gamma| -- silently
  // masked at every M-step-consistent point by a cancelling error in the
  // omitted quadratic term, and therefore wrong everywhere else (verified
  // against a brute-force multivariate-normal evaluation: the shortcut
  // could be off by thousands of log-lik units even at full convergence on
  // real data). The general formula here has no such restriction, and needs
  // no special-casing for sparsity > 0 either: it never assumes Omegaq is
  // Sigma_hat's exact inverse, so it stays exactly log p(Y; theta) however
  // Omegaq was obtained (plain inversion or graphical lasso).
  double objective() const override {
    arma::mat R = data_.Y - data_.X * B_;
    arma::mat Gamma_inv = Omegaq_ + arma::diagmat(C_.t() * dm1_);
    auto Gamma_fresh = nb_utils::inv_and_log_det_sympd(Gamma_inv);
    arma::mat dm1C = C_;
    dm1C.each_col() %= dm1_;
    arma::mat Mu_fresh = R * dm1C * Gamma_fresh.inv;

    double log_det_Omegaq = arma::log_det_sympd(Omegaq_);
    double sum_log_dm1    = arma::sum(arma::log(dm1_));
    double SSQ_w = arma::dot(dm1_, arma::vectorise(arma::sum(arma::square(R), 0)));
    double quad  = SSQ_w - arma::trace(Gamma_inv * (Mu_fresh.t() * Mu_fresh));

    double J = -0.5 * data_.n * data_.p * std::log(2.0 * arma::datum::pi);
    J += 0.5 * data_.n * sum_log_dm1;
    J += 0.5 * data_.n * log_det_Omegaq;
    J -= 0.5 * data_.n * Gamma_fresh.log_det; // -log|Gamma^{-1}| = +log|Gamma|
    J -= 0.5 * quad;
    return J;
  }

  const arma::mat& Gamma() const { return Gamma_; }
  const arma::mat& Mu() const { return Mu_; }

  std::unique_ptr<NormalBlockBase> clone() const override {
    return std::make_unique<NormalBlockKnownClusters<NoisePolicy>>(*this);
  }
  void restore_from(const NormalBlockBase& other) override {
    copy_tracked_state_from(other);
    const auto& o = static_cast<const NormalBlockKnownClusters<NoisePolicy>&>(other);
    Gamma_ = o.Gamma_;
    Mu_ = o.Mu_;
  }
  // Validated on real data (see git history) for sparsity_ <= 0. Excluding
  // sparsity_ > 0 is no longer about objective() validity (it is the
  // general marginal log-likelihood for *any* PD Omegaq, sparse or not, so
  // the objective-comparison gate in try_squarem_step() is just as sound
  // here in principle) -- it was re-tested directly and rejected on cost/
  // benefit: estimate_omega() calls glassoFast, an *approximate* iterative
  // solver, so even *plain* EM's M-step isn't a guaranteed ascent step when
  // sparsity_ > 0 (a separate, pre-existing issue -- plain EM alone already
  // shows tiny non-monotonicity here, see git history). Letting SQUAREM
  // jump into a region of parameter space plain EM would never have
  // visited amplifies that pre-existing instability (observed: real
  // increment violations up to a few units, ~2x worse than plain EM's own,
  // on `brca_rppa`), while the achieved speedup is much more modest than
  // the unpenalized case (the regularization already smooths out most of
  // the slow-linear-EM regime this acceleration targets -- observed ~1.1x-
  // 1.8x here vs 5x-10x+ for sparsity_ <= 0). Not worth the added risk for
  // the gain; kept disabled.
  bool supports_acceleration() const override { return sparsity_ <= 0.0; }
};

#endif // NORMALBLOCKR_NORMAL_BLOCK_KNOWN_CLUSTERS_H
