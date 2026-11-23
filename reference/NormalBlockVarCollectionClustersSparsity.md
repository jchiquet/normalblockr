# Collection of Normal-Block-Var Models over Cluster Counts and Sparsity Levels

R6 class for a collection of normal-block models with different number
of clusters (q) and different sparsity levels, i.e. one
[NormalBlockVarCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionSparsity.md)
per q.

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionClustersSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.md)
-\> `NormalBlockVarCollectionClustersSparsity`

## Methods

### Public methods

- [`NormalBlockVarCollectionClustersSparsity$clone()`](#method-NormalBlockVarCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClustersSparsity$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-get_best_model)
- [`NormalBlockCollectionClustersSparsity$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-get_model)
- [`NormalBlockCollectionClustersSparsity$initialize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-initialize)
- [`NormalBlockCollectionClustersSparsity$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-plot)

------------------------------------------------------------------------

### `NormalBlockVarCollectionClustersSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionClustersSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:3, sparsity = TRUE,
                       control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
models$get_best_model("BIC")$q
#> [1] 3
```
