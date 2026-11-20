# Collection of Normal-Block-Var Models over a Range of Cluster Counts

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level. The q-independent part of
the initial clustering heuristic is computed once for the whole
collection.

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.md)
-\> `NormalBlockVarCollectionClusters`

## Methods

### Public methods

- [`NormalBlockVarCollectionClusters$clone()`](#method-NormalBlockVarCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClusters$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-get_best_model)
- [`NormalBlockCollectionClusters$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-get_model)
- [`NormalBlockCollectionClusters$initialize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-initialize)
- [`NormalBlockCollectionClusters$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-optimize)
- [`NormalBlockCollectionClusters$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-plot)
- [`NormalBlockCollectionClusters$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-print)
- [`NormalBlockCollectionClusters$refine()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-refine)

------------------------------------------------------------------------

### `NormalBlockVarCollectionClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:4, control = NB_control(verbose = FALSE))
models$get_best_model("ICL")$q
#> [1] 4
```
