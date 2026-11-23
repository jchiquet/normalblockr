# Plot a Collection of Normal-Block Models

Plots the model-selection criteria of a collection of normal-block
models: against the number of clusters q for a collection over q
([NormalBlockVarCollectionClusters](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionClusters.md),
[NormalBlockMeanCollectionClusters](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionClusters.md)),
against the penalty for a sparsity path
([NormalBlockVarCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionSparsity.md),
[NormalBlockMeanCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionSparsity.md)),
with a dashed line at the model selected by each criterion; as a heatmap
over q and the penalty for a collection over both
([NormalBlockVarCollectionClustersSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionClustersSparsity.md),
[NormalBlockMeanCollectionClustersSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionClustersSparsity.md)).

## Usage

``` r
# S3 method for class 'NormalBlockCollection'
plot(x, ...)
```

## Arguments

- x:

  An object inheriting from NormalBlockCollection.

- ...:

  passed to the collection's `$plot()` method: `criteria` (the criteria
  to plot, among `"deviance"`, `"ICL"`, `"BIC"` and `"EBIC"`) for a
  collection over q or a sparsity path, plus `log.x` for the latter;
  `criterion` (a single one) and `n_intervals` for a collection over
  both.

## Value

Invisibly returns the
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object; called for its side effect of plotting.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
plot(models)

plot(models, criteria = c("BIC", "ICL"))
```
