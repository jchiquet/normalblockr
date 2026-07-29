# Changelog

## normalblockr 0.2.0

First CRAN submission

### New features

- Zero-inflation extension
  (`ZINormalBlockVarKnownClusters`/`ZINormalBlockVarUnknownClusters`)
  for data with an excess of exact zeros.
- Sparsity path on the cluster-level precision matrix (graphical lasso),
  with warm-starting across penalties
  (`NormalBlockVarCollectionSparsity`, `sparsity = TRUE` in
  [`normal_block()`](../reference/normal_block.md)).
- Accelerated variational EM (SQUAREM-style extrapolation) for both
  known- and unknown-clustering models.
- Several clustering-initialization heuristics (`ward2`, `kmeans`,
  `spectral`, `sbm`, selectable via `NB_control(clustering_init = )`),
  and `best_of_inits()` to try several and keep the best-ELBO fit.
- `refine()` on `NormalBlockVarCollectionClusters`: post-hoc split/merge
  search seeded from neighboring cluster counts, to escape mediocre
  local optima left by independent per-q cold starts.
- New real datasets: `brca_rppa` (breast cancer proteomics), `onema`
  (French stream fish biomass, zero-inflated), `university` (WebKB text
  data), each with a dedicated vignette.

### Other changes

- Renamed the model classes (`NormalBlock*` to `NormalBlockVar*`) for
  consistency; `NormalBlockData` is unchanged.
- `NormalBlockData` rescales columns of `Y` by default; fitted values
  and regression coefficients are reported back on the original scale.
- Removed the `ClustOfVar` dependency (the `kmeansvar` clustering
  heuristic was dropped after benchmarking showed it was both the
  worst-ranked and least reliable of the available heuristics).
- Package cleanup for CRAN submission: license, documentation, and
  package structure.

## normalblockr 0.1.0

- Initial implementation of the Normal-Block model: a Gaussian graphical
  model with a latent clustering structure, for known or unknown
  clusterings, fit by variational EM.
