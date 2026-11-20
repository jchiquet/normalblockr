# Collection of Normal-Block-Var Models over a Sparsity Path

R6 class for a collection of normal-block models with a fixed clustering
(blocks) and different sparsity levels. The default penalty grid is read
off the latent block covariance of a short unpenalized fit.

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.md)
-\> `NormalBlockVarCollectionSparsity`

## Methods

### Public methods

- [`NormalBlockVarCollectionSparsity$clone()`](#method-NormalBlockVarCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-get_best_model)
- [`NormalBlockCollectionSparsity$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$initialize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-initialize)
- [`NormalBlockCollectionSparsity$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-plot)
- [`NormalBlockCollectionSparsity$stability_selection()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-stability_selection)

------------------------------------------------------------------------

### `NormalBlockVarCollectionSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = ex$parameters$C, sparsity = TRUE,
                       control = NB_control(verbose = FALSE, n_sparsity_penalties = 5))
models$get_best_model("BIC")$sparsity
#> [1] 0.06648302
```
