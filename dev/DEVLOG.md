# normalblockr development log

Design decisions and the evidence behind them, for the developers. Users
should not need this file: the models and the algorithmic choices are
documented in `inst/normal_block_models.qmd` (published as the technical
report on the pkgdown site), the user-visible changes in `NEWS.md`.

The design studies cited below live in the article repository
(`normalblockr_jss`, `normalblockr-jss-article-quarto/experiments/design/`,
subfolders `mean_block/` and `initialization/`), with their results.

## Open issues, for a patch release (1.0.x)

Found during the 1.0.0 review, deliberately left unchanged in 1.0.0 (they
change numerical results, and the article's replication code had to stay
stable).

1. **Penalized variance-block objective not monotone.** With `sparsity > 0`
   the (V)EM objective can decrease between iterations, by up to about 0.025
   (spherical noise, known clusters, `sparsity = 0.05`, on the
   `tests/testthat/testdata/testdata_normal.RDS` fixture), with or without
   SQUAREM. Likely cause: the graphical-lasso M-step is an inexact,
   warm-started iterative solve (`nb_omega::estimate()`), so the M-step is not
   an exact argmax. Leads: compare the penalized objective before/after each
   M-step and fall back to a cold or tighter solve when it decreases; check
   that the objective traced is exactly the one the M-step maximizes
   (penalty included, same threshold).
2. **Small sparse variance-block paths converge slowly.** In the StARS test
   (small n and p, a few penalties), some models hit the 500-iteration cap
   with increments hovering around the 1e-4 threshold (the test now uses
   `threshold = 1e-3`). Probably the same root cause as 1 (the objective
   oscillates at the solver's precision). Leads: as in 1, and a stopping rule
   on the relative change of the parameters rather than of the objective
   when `sparsity > 0`.
3. **`refine()` does not remove every deviance monotonicity violation.**
   The models of a collection are nested, so the deviance should not
   increase with q. After `refine()`, about 1.3-1.8 violations per
   collection remain on simulated data (`var_init_study.R`), one on
   `brca_rppa` for the mean-block family (20 before), and on `onema` the
   `spectral` initialization goes from 1.7 to 2.6 violations (improving one
   q can expose its neighbour). The documentation only claims that
   `refine()` repairs poor local optima, not that it removes every
   violation. Leads: iterate `refine()` until no model improves, or retry
   the remaining violating q from both neighbours with more candidates.

Minor, same family of issues: SQUAREM leaves a decrease of about 4e-6 on an
unknown-cluster spherical variance-block fit (an accepted extrapolation
slightly below the previous objective).

(Fixed before the release: zero-inflated spherical models started from
per-variable variances; they now start from the pooled masked estimate.)

## Model-selection criteria (1.0.0, #55, checked in #59)

- **Per-group BIC penalties** `sum nu_g log N_g`, after Daudin, Picard &
  Robin (2008) and Keribin et al. (2015): alpha `log p`; mean-block profiles
  and a spherical noise variance `log(np)`; the rest `log n`. The `log(np)`
  on the mean-block profiles is slightly conservative with weakly separated
  clusters (under-selection), kept for consistency with the latent block
  model literature.
- **ICL = BIC + 2 H(tau)** only. The former ICL also added the Gaussian
  entropy of W, which made it differ from the BIC for a known clustering,
  depend on the scale of the data, and select q much worse (variance-block
  ICL 63% -> 89% of correct selections on separable simulated data with
  1.0.0; the exact ICL, W integrated out, does no better: 78%,
  `criteria_selection_study.R`).
- **EBIC** `BIC + 2|E| log(dim)`, dim = q (var) or p (mean, which used q
  before 1.0.0).
- **Open question, deliberately unchanged**: under `sparsity > 0` with a full
  mean-block Sigma, `nb_param` counts the non-zero off-diagonal terms of Omega
  (lasso df), which makes the BIC non-monotone in q and over-select; counting
  the dense parameters would be monotone but over-charge
  (`nestor_protocol.R`, unpenalised arm). The diagonal default is unaffected.

## Rescaling of Y per family (1.0.0, #60)

Data simulated from either model are misspecified once rescaled: every
criterion over-selects (mean-block 76-92% over-selection on separable
datasets, variance-block 27%). On real data the exact ICL on the original
scale (Jacobian included) favors rescaling for the variance-block family
(`brca_rppa`, `onema`; also the cross-validated log-likelihood,
`var_scale_cv_real_data.R`, rerun with 1.0.0: +2,074 (se 124) on brca,
+83 (se 66, n.s.) on onema, n.s. on university) and raw data for the mean-block family. Hence `scale = NULL`,
resolved per family.

## Initialization: kmeans + refine() for both families (1.0.0, #65)

- Mean-block: kmeans is the best of the four heuristics on `brca_rppa`,
  `onema`, a simulation, and Nestor's M2 grid (mean ARI 0.87 at SNR 0.3,
  0.98 at SNR 1.5 with 1.0.0; sbm close, ward2/spectral behind) (`nestor_protocol.R`,
  `diagnostic_mean_block_real_data.R`).
- Variance-block (`var_init_study.R`): heuristics equivalent on simulations;
  on real data kmeans + refine() best on brca (deviance gap to the best fit
  summed over q: 2,409 vs 8,624 for ward2 with 1.0.0), close on university, but ward2
  better on zero-inflated onema (197 vs 590), where best_of_inits is best. The 0.3.0 default was ward2,
  chosen on cold-started collections, where kmeans had more deviance
  monotonicity violations and a seed-dependent ICL choice: refine() is what
  fixes that, so do not turn refine off by default without revisiting the
  initialization.
- `kmeansvar` (ClustOfVar) ranked worst everywhere and was dropped, with the
  ClustOfVar dependency.
- A single SBM exploration over the whole q range is 15-90x cheaper than one
  per q, for the same ICL after the (V)EM.

## refine() (1.0.0, #56)

- Variance-block merges ranked by the signed correlation of the latent
  signals (was |Omega_kl|, which tried anticorrelated pairs first); splits by
  a 2-means on the residuals (was the median noise rule). On brca, refine()
  lowered the deviance summed over q = 1..30 by 15,758 (old rules), 18,588
  (signed merges), 22,166 (signed merges + kmeans splits), selected q 27 ->
  22, and ran faster.
- Mean-block split rule kept at "noise": 1.7% more deviance removed on brca
  (5 seeds), < 1% on onema (3 seeds), tie on simulations, same selected q and
  ARI (`split_rule_comparison.R`).
- Cost: refine() is 3 to 6 times the cold-started collection.

## Mean-block noise covariance: diagonal by default

A full Sigma's p(p+1)/2 parameters drown the mean structure: BIC picks the
true q 9-10/12 with a diagonal Sigma against 5/12 with a full one at
n/p = 1.3 (1.0.0), even on data generated with a full Sigma; the ARI at the true q is
the same (`covariance_shape_selection_study.R`). On brca, q = 1 -> 15 buys a
full Sigma 514 units of deviance against 97,056 for a diagonal one
(`deviance_monotonicity.R`). `sparsity > 0` implies "full"; zero-inflated
mean-block models have no full Sigma (each row would need its own submatrix
inverse), hence no sparsity.

## SQUAREM policy (1.0.0)

- Variance-block: decisive at large q (q = 40, p = 200: 182 iterations and
  2.7 s, against the 5,000-iteration cap without convergence in 26 s);
  roughly 2-3x fewer iterations for q <= 10, 5-9x at q = 25-30.
- Mean-block with a diagonal or spherical Sigma: never reduces the number of
  iterations (exact EM, a few dozen iterations), and costs 1.1 to 2.2x the
  time through rejected extrapolations: off by default.
- Never with `sparsity > 0`: the graphical lasso is an approximate M-step,
  and extrapolation amplified its non-monotonicity for a modest gain.

## Graphical lasso in-package (0.3.0) and shared with PLNmodels (1.0.0)

- Replaced `glassoFast`, called back from C++ into R once per M-step. That
  boundary was behind an intermittent heap corruption (CRAN
  r-devel-fedora-gcc on 0.2.1, ~10% of CI runs, never reproduced locally).
  Confidence in the fix is structural: the mechanism is gone.
- Two glassoFast bugs fixed: with no off-diagonal mass (every q = 1 problem)
  it returned ~9e15 instead of 1/S_ii; a NaN input hung its inner loop.
  Agreement elsewhere: 477 random problems (worst relative difference
  2.2e-13) and 40 (Sigma, rho) pairs captured from a real path (identical
  number of inner passes).
- Warm starts between M-steps: on 124 consecutive M-steps of a real path,
  1.355 s at 6.1e-5 from the exact solution cold (threshold 1e-4), 0.589 s
  at 8.9e-7 warm (threshold 1e-6). An earlier rejection of warm starts came
  from a synthetic micro-benchmark. A warm start can diverge on an
  ill-conditioned Sigma: `nb_omega::estimate()` retries cold before falling
  back to the unpenalized inverse (e591f28); keep that retry.
- Shared verbatim with PLNmodels (1.0.0): solved scaled to a unit diagonal
  (median number of sweeps from 646 to 4 on PLN residual covariances),
  indefinite results shifted on the diagonal rather than eigenvalue-clamped
  by `ensure_pd()` (which densified them). Change it in both packages
  together.

## Closed-form ZI updates and other implementation history

- The ZI updates of B and M were first implemented with `nloptr`'s L-BFGS
  inside the (V)EM; they are quadratic, now solved in closed form, and the
  nloptr dependency is gone.
- SPD solves with a reused buffer: 1.7-2x on a whole zero-inflated fit at
  q = 4..25.
- Earlier R6 revisions reused a one-iteration-stale residual in three of the
  four variance-block variants; fixed on both sides.
- `objective()` used to read back the previous E-step's Sigma-hat and
  log|Gamma|, one iteration behind; harmless for the profiled ELBO, wrong
  for the general criterion at a SQUAREM candidate.
- The mean-block tau update was first simultaneous over the rows: it fell
  into a period-2 cycle; it is now a sequential sweep.
- The singular weighted normal equations surfaced on `onema`, whose rarest
  species (seen at one station) drove dm1 to infinity; those species are now
  filtered out of the dataset, and the SSQ is floored.

## Sequential mean-then-variance fits

On a positive control the two structures are recovered exactly (ARI = 1 on
each); on brca the two partitions are unrelated (ARI 0.017, q = 50 in the
mean, 19 in the covariance), while the covariance clustering is close to the
direct one (ARI 0.715; with 1.0.0: ARI 0.018 and 0.644, q = 65 and 20) (`sequential_mean_then_variance.R`).

## Profiling (2026-09, -O2 builds)

Two regimes: on large real data (brca) the C++ (V)EM dominates; on small
problems initialization dominates (up to 90%, ZI mean-block). A mean-block
sparsity path over a full p x p Sigma is ~96% graphical lasso: the only lever
is calling it less (a 5-penalty path on brca takes ~36 s, never in a
vignette). The ZI variance-block E-step does one q x q `inv_sympd` per
observation (intrinsic). kmeans `nstart = 30` is a quality/speed dial.
