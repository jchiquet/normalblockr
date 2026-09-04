# Collection of Normal-Block Models over Cluster Counts and Sparsity Levels

R6 class for a collection of normal-block models with different number
of clusters (q) and different sparsity levels.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionClustersSparsity`](NormalBlockCollectionClustersSparsity.md)
-\> `NormalBlockVarCollectionClustersSparsity`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionClustersSparsity$new()`](#method-NormalBlockVarCollectionClustersSparsity-initialize)

- [`NormalBlockVarCollectionClustersSparsity$clone()`](#method-NormalBlockVarCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClustersSparsity$get_best_model()`](NormalBlockCollectionClustersSparsity.html#method-get_best_model)
- [`NormalBlockCollectionClustersSparsity$get_model()`](NormalBlockCollectionClustersSparsity.html#method-get_model)
- [`NormalBlockCollectionClustersSparsity$plot()`](NormalBlockCollectionClustersSparsity.html#method-plot)

------------------------------------------------------------------------

### `NormalBlockVarCollectionClustersSparsity$new()`

Create a new \[\`NormalBlockVarCollectionClustersSparsity\`\] object.

#### Usage

    NormalBlockVarCollectionClustersSparsity$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `q_list`:

  list of q values (number of groups) in the collection

- `zero_inflation`:

  boolean to specify whether data is zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockVarCollectionClustersSparsity\`\] object

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
