# normalblockr 0.2.1

* S3 methods `print()`, `summary()`, `plot()`, `logLik()` and `BIC()` for fitted models (any `NormalBlockVarBase` subclass), and `print()`/`summary()`/`logLik()`/`BIC()` for collections of models; accessing `$loglik` on a collection now raises an informative error instead of silently returning `NULL`.
* Addressed reviewer feedback ahead of CRAN submission: shortened/title-cased man page titles, added missing `@examples`, cross-referenced `$plot_network()`/`$latent_network()` in each other's documentation, and replaced a few inefficient matrix operations (`solve()` on symmetric positive-definite matrices, `M %*% t(C)`) with `chol2inv(chol())`/`tcrossprod()`.

# normalblockr 0.2.0

First CRAN submission

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
