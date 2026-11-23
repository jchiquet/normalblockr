# normalblockr 1.0.0

## Changes of default behavior since 0.3.0

Users of 0.3.0 will see different results on the same calls:

* **Model-selection criteria.** The BIC charges each group of parameters the
  log of the number of independent observations informing it (Daudin, Picard
  & Robin, 2008; Keribin et al., 2015), instead of `log(n)` for all: the
  cluster proportions `log(p)`, the mean-block cluster profiles and a
  spherical noise variance `log(n * p)`, the rest `log(n)`. The ICL completes
  the partition only (Biernacki, Celeux & Govaert, 2000), `ICL = BIC + 2 *
  entropy` with `entropy = -sum(tau * log(tau))`, so that it equals the BIC
  for a known clustering in both families. The EBIC of mean-block models
  counts its edges among the p variables (it used q). The selected number of
  clusters can change.
* **Rescaling of Y depends on the model family.** `NormalBlockData$new()`
  gains `scale = NULL` as its default: variance-block models keep rescaling
  each column of Y by its standard deviation, mean-block models no longer
  do (rescaling breaks the mean profile shared within a block). `TRUE` or
  `FALSE` forces either for both families.
* **One initialization for both families, followed by `refine()`.** The
  default `clustering_init` is `"kmeans"` for both families (it was
  `"ward2"` for variance-block models), and collections over a range of q
  are polished by the split and merge search of `refine()` after their
  cold-started fit (`NB_control(refine = TRUE)`; `FALSE` restores the former
  behavior). `"ward2"` remains a good choice for zero-inflated abundances.
* **SQUAREM off for mean-block models with a diagonal or spherical Sigma.**
  Their EM is exact and converges in a few dozen iterations, so the
  extrapolations only cost time. `NB_control(acceleration = NULL)`, the new
  default, keeps SQUAREM for variance-block models and full-Sigma mean-block
  models; `TRUE` or `FALSE` forces it.
* **Extractors on the original scale of Y, with names.** `coef()` of
  variance-block models returns B on the original scale of Y (it returned it
  on the rescaled scale of the fit, unlike `fitted()` and `predict()`), and
  so does `sigma()` for the p x p covariance of mean-block models. Quantities
  defined per cluster (variance-block `sigma()`, mean-block `coef()`) stay on
  the scale of the fit. `coef()`, `sigma()`, `fitted()`, `predict()` and
  `latent_network()` keep the names of the variables and covariates, and
  name the clusters `Cluster_1`, ... (also in the network plots and the
  StARS stability path).

## New features

* **Held-out log-likelihood and cross-validation.** `heldout_loglik(newdata)`
  on every fitted model gives the exact log-likelihood of new observations at
  the fitted parameters and the MAP partition, on Y's original scale (for
  variance-block models, the marginal over the latent W, with q x q matrices
  only). `normal_block_cv()` computes the K-fold cross-validated
  log-likelihood of a collection over q, a penalty-free criterion that
  compares families, shapes of Sigma and scalings on the same footing, and
  selects q by the one-standard-error rule; it has `print()` and `plot()`
  methods.
* **StARS for mean-block sparsity paths.** `get_best_model("StARS")` and
  `stability_selection()` apply to both families, the clustering being held
  fixed on the subsamples.
* **Sequential fits from a fitted mean-block model.**
  `normal_block_sequential()` accepts, as `blocks_mean`, a mean-block model or
  collection already fitted on the same data, which is then not refitted.
* **Diagnostics.** `within_block_correlations()` compares, for each block of a
  variance-block model, the empirical correlations of the residuals with the
  positive correlation the model implies: a block holding anticorrelated
  variables, typically from a known clustering taken from an annotation, is
  misspecified. Collections over q gain a `deviance_violations` field (values
  of q whose deviance exceeds the previous one, a sign of a poor local
  optimum), reported by `print()`.
* **Split rules in `refine()`.** `split()` gains a `rule` argument and
  `refine()` / `candidates_split()` a `split_rule` one: `"kmeans"` (a 2-means
  on the members' residuals or least-squares trajectories, the default for
  the variance-block family) or `"noise"` (a cut at the median of the noise
  variances, the default for the mean-block family). Variance-block merge
  candidates are ranked by the signed correlation between the latent block
  signals, since a block can only gather positively co-varying variables.
* `plot()` method for collections of models: the criteria against q or the
  penalty, or a heatmap over both (it wraps the collections' `$plot()`). The
  model selected by each criterion is marked on its curve and labelled with
  the criterion and the selected q or penalty, so that it is clear which
  dashed line belongs to which criterion.
* New fields `penalty_terms` and `penalty` (the BIC penalty, group by
  group), `NormalBlockData$for_family()` (the data a family is fitted on),
  `NormalBlockData$subset()` and `$Y_original`, and `sparsity_details` on
  mean-block sparsity paths.
* The graphical lasso (`src/graphical_lasso.h`) is shared verbatim with
  PLNmodels. It solves the problem scaled to a unit diagonal, which makes it
  converge on covariances whose variances span orders of magnitude, shifts an
  indefinite precision matrix on its diagonal (keeping its support), detects
  limit cycles and can be interrupted. `graphical_lasso_fit()` also returns
  `status`, `delta` and `shift`.

## Bug fixes

* StARS: the subsample paths reuse the collection's penalties (they derived
  their own grid from each subsample, so the selection probabilities were
  filed under penalties that had not produced them), and keep the
  zero-inflation design `X0` and the `scale` request.
* `plot_network()` computes the layout without the signed edge weights, so
  that force-directed layouts such as `igraph::layout_with_fr` work; `layout`
  may also be a matrix of coordinates. Mean-block networks are labeled with
  the variable names.
* `model_par` exposes the cluster proportions `alpha` of mean-block models
  with an unknown clustering, and `B0` of zero-inflated mean-block models, as
  in the variance-block family.
* `NB_control(clustering_init = "best_of_inits")` with a known clustering
  crashed: the known clustering is now simply fitted.
* A sparsity path of variance-block models failed when `sparsity_weights`
  left some pairs unpenalized; such pairs are now ignored by the default
  penalty grid.
* A constant column of Y made the variance-block fit diverge, and `refine()`
  then failed: `NormalBlockData$new()` now refuses such columns, naming them
  (`normal_block_cv()` names the fold where it happens). A split or merge
  candidate whose fit fails is discarded by `refine()` instead of aborting
  it.
* `generate_normal_block_var_data()` and `generate_normal_block_mean_data()`
  honor `alpha`, and could put a variable in two clusters when the first
  draw left a cluster empty. Data simulated for a given seed are unchanged
  otherwise.
* StARS subsamples keep an explicit `zeros` mask and the formula.
* A collection over an unsorted range of q is ordered by q.
* Known clusterings are checked to have one row per variable.
* StARS on a variance-block sparsity path with an unknown clustering
  failed when the least penalized model had a nearly empty cluster: an
  initial clustering given as probabilities now only needs some mass in
  every cluster.
* Zero-inflated models with a spherical noise start from the spherical
  estimate of the noise variance, pooled over the non-zero cells (the
  variance-block ones started from per-variable variances).

## Other changes

* The help pages render their links, code and emphasis (roxygen markdown
  was not enabled: links to other pages showed as plain text in brackets);
  the links of the online reference to inherited R6 methods reach the
  method.
* Every help page has a running example; help titles follow the two families
  of the article (Normal-Block-Var, Normal-Block-Mean).
* `summary()` of a model and `print()` of a sequential fit wrap their long
  lines at the console width; the descriptions of the collections are the
  same for both families.
* The progress line of a collection, rewritten in place, is only printed in
  an interactive session (it made one long line in documents and logs).
  `print()` of a collection no longer suggests `refine()` for a deviance
  increase that `refine()` has already failed to repair.
* Printed outputs speak of clusters of variables rather than blocks.
* Vignettes and README build a known clustering with
  `model.matrix(~ 0 + factor(labels))`.
* The update of the mean-block cluster profiles solves a linear system
  instead of inverting a matrix. Results are unchanged.
* `NB_control(threshold = )` negative runs exactly `niter` iterations,
  without the warning about the iteration cap.
* `fixed_tau` is read from the control list by both families (the public
  `$fixed_tau` field of variance-block models is removed); `update()` loses
  the arguments `Psi`, `Phi`, `Lambda` and `kappa`, which it stored but never
  used.
* The packages dplyr, tidyr, tibble and scales are no longer imported, and
  Metrics is no longer suggested. The code shared by the two families and by
  their variants lives in their base classes; results are unchanged.

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
