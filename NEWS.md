# normalblockr (development version)

## One initialization for both families, and refine() by default

Fits over a range of cluster counts change: they are refined, and the
variance-block ones start from another heuristic.

* `NB_control(refine = TRUE)` is now the default, for collections of both
  families: every collection over a range of `q` is polished by `refine()`
  after its cold-started fit (typically 3 to 6 times the cost). A one-line
  report says how many models it improved. `refine = FALSE` restores the
  previous behaviour; heuristic fits (`heuristic = TRUE`) are never refined,
  and `refine()` called on one now says why it cannot run.
* The default `clustering_init` of variance-block models is now `"kmeans"`,
  as for mean-block ones. On simulated data the heuristics are equivalent;
  on real data `kmeans` followed by `refine()` gives the best fits on
  `brca_rppa` and close to the best on the university data, while `"ward2"`
  remains better on zero-inflated abundances such as `onema` (see
  `inst/clustering_initialization_benchmark/var_init_study.R`).

## Rescaling of Y now depends on the model family

The fits of mean-block models on default data change.

* `NormalBlockData$new()` gains `scale = NULL` as its default: the model
  family decides when a model is fitted. Variance-block models keep
  rescaling each column of Y by its standard deviation; mean-block models no
  longer do. `scale = TRUE` or `FALSE` forces either for both families.
* Rescaling breaks the mean-block model, whose variables share a mean profile
  within a block: on `brca_rppa` and `onema` the unscaled mean-block fit has
  the better exact ICL on the original scale, while the rescaled one wins for
  the variance-block model. On data simulated from either model, rescaling
  makes every criterion over-select the number of clusters.
* `normal_block_sequential()` hands the residuals of its mean-block stage to
  the variance-block stage with the same `scale` request, so each stage
  applies its own default.
* New method `NormalBlockData$for_family()` returns the data a model of a
  given family is fitted on.

## Diagnostics

* `within_block_correlations()` on variance-block models: for each block, the
  empirical correlations of the residuals `Y - XB` (mean, minimum, proportion
  of negative ones) against the correlation the fitted model implies. Within a
  block two variables covary positively by construction (their covariance is
  `Sigma_kk > 0`), so a block holding anticorrelated variables, typically from
  a known clustering taken from an annotation, is misspecified.
* Collections over a range of q gain a `deviance_violations` field, the
  values of q whose deviance exceeds the one at the previous q, and `print()`
  warns about them. The models are nested, so every violation signals a poor
  local optimum, which `refine()` is designed to repair.
* Vignettes and README build a known clustering with
  `model.matrix(~ 0 + factor(labels))` rather than the internal
  `normalblockr:::as_indicator()`.

## SQUAREM can be turned off

* `NB_control(acceleration = FALSE)` runs plain (V)EM instead of its
  SQUAREM-accelerated version, in both families and for every variant (the
  acceleration is on by default, and never attempted when `sparsity > 0`).
  Useful to measure what the acceleration buys: on a variance-block model with
  8 unknown clusters, plain VEM needs about ten times as many iterations.

## Model-selection criteria aligned with Biernacki and Daudin

The criteria change, so the number of clusters they select can change too.

* **BIC: one penalty per group of parameters.** Each group of parameters is
  charged the log of the number of independent observations informing it,
  following the ICL of Daudin, Picard & Robin (2008) for SBMs and Keribin et
  al. (2015) for latent block models, instead of `log(n)` for all of them:
  the cluster proportions `log(p)` (the variables are the units being
  clustered), the mean-block cluster profiles and a spherical noise variance
  `log(n * p)` (shared by every cell of a block), the rest `log(n)`. The new
  fields `penalty_terms` (one row per group: `nu`, `N`) and `penalty` expose
  the computation; `nb_param` is unchanged.
* **ICL completes the partition only**, as in Biernacki, Celeux & Govaert
  (2000): `ICL = BIC + 2 * entropy` with `entropy = -sum(tau * log(tau))`.
  The variance-block models used to add the Gaussian entropy of the latent
  vectors W, which made the ICL differ from the BIC even for a known
  clustering and depend on the scale of the data. It now equals the BIC for a
  known clustering, in both families.
* **EBIC on the right network.** The EBIC term is `2 * n_edges * log(dim)`,
  with `dim = q` for the variance-block family as before, but `dim = p` for
  the mean-block family, whose network relates the p variables (it used `q`).
* On `brca_rppa`, the mean-block ICL now selects 50 clusters instead of 70
  over `seq(5, 120, 5)`; the variance-block BIC selection is unchanged.

## Better split and merge candidates in `refine()`

* Variance-block merge candidates are ranked by the **signed** correlation
  between the two blocks' latent signals instead of `|Omega_kl|`. Within a
  block, two variables necessarily covary positively (their covariance is
  `Sigma_kk > 0`), so merging two anticorrelated blocks is never a sensible
  candidate; the former ranking tried those first.
* `split()` gains a `rule` argument, and `refine()` / `candidates_split()` a
  `split_rule` one: `"kmeans"` divides a cluster by a 2-means on its members'
  columns of residuals (variance-block) or least-squares trajectories
  (mean-block), `"noise"` at the median of their noise variances (the former
  rule). The default is family-specific: `"kmeans"` for the variance-block
  family, `"noise"` for the mean-block one, where the two rules are close
  and `"noise"` is slightly ahead on real data
  (`inst/mean_block_analyses/split_rule_comparison.R`).
* On `brca_rppa`, `refine()` now lowers the variance-block deviance summed
  over `q = 1:30` by 22,166 instead of 15,758, and runs faster; the selected
  `q` moves from 27 to 22.

## Graphical lasso shared with PLNmodels

* `src/graphical_lasso.h` is now the very same file as in PLNmodels, which
  ported the same glassoFast Fortran: the two copies had drifted apart. It is
  shared verbatim (namespace `graphical_lasso`), and brings to normalblockr:
  * **the problem is solved scaled to a unit diagonal**, an exact change of
    variables. On a covariance whose variances span orders of magnitude, the
    unscaled descent could stall in a limit cycle, fail in its inner loop, or
    return an indefinite precision matrix; scaled, it converges in a few sweeps
    (median number of sweeps from 646 to 4 on PLN residual covariances). On
    ordinary input the result is unchanged up to the tolerance;
  * **an indefinite precision matrix is shifted on its diagonal** to positive
    definite, which keeps its support, where `ensure_pd()` clamped its
    eigenvalues and made it dense;
  * the detection of limit cycles, as a safeguard, and interrupt checks, so
    that a long solve can be stopped from R (by the user or by a time limit).
* `graphical_lasso_fit()` also returns `status`, `delta` and `shift`.

# normalblockr 0.3.0

## New model family: clustering in the mean

The package now fits two complementary families. The original Normal-Block
model clusters variables by how they *covary*; the new mean-block family
clusters them by how they *respond* to the covariates, constraining the mean
as `mu_i = C B' X_i` with one regression profile per cluster.

* `normal_block(..., model = "mean")` fits it, for a known clustering
  (`NormalBlockMeanKnownClusters`) or an unknown one inferred by variational
  EM (`NormalBlockMeanUnknownClusters`), over a range of cluster counts
  (`NormalBlockMeanCollectionClusters`), a penalty path
  (`NormalBlockMeanCollectionSparsity`) or both
  (`NormalBlockMeanCollectionClustersSparsity`).
* Zero-inflated counterparts (`ZINormalBlockMeanKnownClusters`,
  `ZINormalBlockMeanUnknownClusters`), and collections of those over a range
  of cluster counts.
* `noise_covariance` gains a `"full"` shape for this family, whose Sigma is
  the p x p residual covariance rather than a diagonal noise term.
  `"diagonal"` is the default: a full Sigma spends p(p+1)/2 parameters that
  drown the mean structure the criteria are weighing, and selects the number
  of clusters markedly worse as p approaches n. Asking for `sparsity > 0`
  implies `"full"`, since a diagonal precision matrix has nothing for the
  graphical lasso to penalize.
* `normal_block_sequential()` chains the two families: a mean-block fit, then
  a variance-block fit on its residuals. A heuristic two-stage estimator, not
  a joint model.
* `NormalBlockData$new()` gains a `zeros` argument, to carry an explicit
  zero-inflation mask when the matrix handed to a model is no longer the one
  carrying the zeros (the residuals of a first stage, say).

## The graphical lasso is now in-package

This also fixes the segfault CRAN reported on
`r-devel-linux-x86_64-fedora-gcc` when re-building
`breast-cancer-proteomics.Rmd` (`address 0x580, cause 'memory not mapped'`),
which our own CI had been hitting intermittently on Linux.

* `glassoFast` is no longer a dependency. The C++ (V)EM used to call back into
  R once per M-step to reach it -- thousands of round trips from inside a
  single `.Call` for one sparsity path, and a recurring source of intermittent
  crashes. The solver now lives in `src/graphical_lasso.h`.
* Two bugs in `glassoFast` are fixed on the way. When the empirical covariance
  has no off-diagonal mass the problem separates and the solution is diagonal;
  it returned `1 / max(rho_ii, eps)` there, dropping the variance term, so with
  an unpenalized diagonal it returned about 9e15 instead of `1 / S_ii` -- and
  every `q = 1` problem takes that branch. Separately, its inner loop could not
  terminate on non-finite input.
* Each M-step warm-starts from the previous one, at a tightened threshold:
  faster *and* closer to the optimum than the previous cold starts.

## Bug fixes

* Zero-inflation with more than one zero-inflation covariate (`X0` with two or
  more columns) failed with "non-conformable arguments": `B0` was built
  transposed. A second bug behind it made `kappa` non-finite when the response
  had no zeros at all. Both were unreachable before, so no working fit changes.
* Fixed two crashes on the penalized path, one of them a long-standing
  memory-safety bug (an R object cached in a C++ `static`).

## Other changes

* A collection over cluster counts no longer repeats work that does not depend
  on the cluster count: the OLS fit, the zero-inflation logistic regressions,
  and the shareable part of each clustering heuristic (one hierarchical tree
  cut at each q, one eigendecomposition, one lossless row compression) are
  computed once. Fitted values are unchanged.
* Printed model descriptions now name the residual-covariance shape for every
  family, and the zero-inflated variance-block models name their family.

# normalblockr 0.2.1

First CRAN submission

* S3 methods `print()`, `summary()`, `plot()`, `logLik()` and `BIC()` for fitted models (any `NormalBlockVarBase` subclass), and `print()`/`summary()`/`logLik()`/`BIC()` for collections of models; accessing `$loglik` on a collection now raises an informative error instead of silently returning `NULL`.
* Shortened/title-cased man page titles, added missing `@examples`, cross-referenced `$plot_network()`/`$latent_network()` in each other's documentation, and replaced a few inefficient matrix operations (`solve()` on symmetric positive-definite matrices, `M %*% t(C)`) with `chol2inv(chol())`/`tcrossprod()`.

# normalblockr 0.2.0

## New features

* Zero-inflation extension (`ZINormalBlockVarKnownClusters`/`ZINormalBlockVarUnknownClusters`) for data with an excess of exact zeros.
* Sparsity path on the cluster-level precision matrix (graphical lasso), with warm-starting across penalties (`NormalBlockVarCollectionSparsity`, `sparsity = TRUE` in `normal_block()`).
* Accelerated variational EM (SQUAREM-style extrapolation) for both known- and unknown-clustering models.
* Several clustering-initialization heuristics (`ward2`, `kmeans`, `spectral`, `sbm`, selectable via `NB_control(clustering_init = )`), and `best_of_inits()` to try several and keep the best-ELBO fit.
* `refine()` on `NormalBlockVarCollectionClusters`: post-hoc split/merge search seeded from neighboring cluster counts, to escape mediocre local optima left by independent per-q cold starts.
* New real datasets: `brca_rppa` (breast cancer proteomics), `onema` (French stream fish biomass, zero-inflated), `university` (WebKB text data), each with a dedicated vignette.

## Other changes

* Renamed the model classes (`NormalBlock*` to `NormalBlockVar*`) for consistency; `NormalBlockData` is unchanged.
* `NormalBlockData` rescales columns of `Y` by default; fitted values and regression coefficients are reported back on the original scale.
* Removed the `ClustOfVar` dependency (the `kmeansvar` clustering heuristic was dropped after benchmarking showed it was both the worst-ranked and least reliable of the available heuristics).
* Package cleanup for CRAN submission: license, documentation, and package structure.

# normalblockr 0.1.0

* Initial implementation of the Normal-Block model: a Gaussian graphical model with a latent clustering structure, for known or unknown clusterings, fit by variational EM.
