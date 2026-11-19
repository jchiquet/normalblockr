# Collection of Normal-Block-Var Models over a Sparsity Path

R6 class for a collection of normal-block models with a fixed clustering
(blocks) and different sparsity levels.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionSparsity`](NormalBlockCollectionSparsity.md) -\>
`NormalBlockVarCollectionSparsity`

## Active bindings

- `sparsity_details`:

  list of information about model's penalties

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionSparsity$new()`](#method-NormalBlockVarCollectionSparsity-initialize)

- [`NormalBlockVarCollectionSparsity$clone()`](#method-NormalBlockVarCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_best_model()`](NormalBlockCollectionSparsity.html#method-get_best_model)
- [`NormalBlockCollectionSparsity$get_model()`](NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$optimize()`](NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](NormalBlockCollectionSparsity.html#method-plot)
- [`NormalBlockCollectionSparsity$stability_selection()`](NormalBlockCollectionSparsity.html#method-stability_selection)

------------------------------------------------------------------------

### `NormalBlockVarCollectionSparsity$new()`

Create a new \[\`NormalBlockVarCollectionSparsity\`\] object.

#### Usage

    NormalBlockVarCollectionSparsity$new(
      mydata,
      blocks,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `blocks`:

  either a clustering matrix (known, fixed clustering) or a single
  integer (number of blocks to infer)

- `zero_inflation`:

  boolean to specify whether data is zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockVarCollectionSparsity\`\] object

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
