# Collection of Normal-Block Models over a Sparsity Path

R6 class for a collection of normal-block models with a fixed clustering
(blocks) and different sparsity levels.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionSparsity`](NormalBlockCollectionSparsity.md) -\>
`NormalBlockVarCollectionSparsity`

## Active bindings

- `sparsity_details`:

  list of information about model's penalties

- `criteria`:

  a data frame with the values of some criteria ((approximated)
  log-likelihood, BIC) for the collection of models

- `stability_path`:

  measure of edges stability based on StARS method

- `stability`:

  mean edge stability along the sparsity penalties path

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionSparsity$new()`](#method-NormalBlockVarCollectionSparsity-initialize)

- [`NormalBlockVarCollectionSparsity$get_best_model()`](#method-NormalBlockVarCollectionSparsity-get_best_model)

- [`NormalBlockVarCollectionSparsity$stability_selection()`](#method-NormalBlockVarCollectionSparsity-stability_selection)

- [`NormalBlockVarCollectionSparsity$clone()`](#method-NormalBlockVarCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_model()`](NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$optimize()`](NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](NormalBlockCollectionSparsity.html#method-plot)

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

### `NormalBlockVarCollectionSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockVarCollectionSparsity$get_best_model(
      crit = c("BIC", "EBIC", "ICL", "StARS"),
      stability = 0.9
    )

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection.

- `stability`:

  if criterion = "StARS" gives level of stability required. Either
  "BIC", "EBIC", "ICL" or "StARS". Default is BIC

#### Returns

a \[\`NormalBlockVarUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockVarCollectionSparsity$stability_selection()`

Compute the stability path by stability selection

#### Usage

    NormalBlockVarCollectionSparsity$stability_selection(
      subsamples = NULL,
      n_subsamples = 10
    )

#### Arguments

- `subsamples`:

  a list of vectors describing the subsamples. The number of vectors (or
  list length) determines the number of subsamples used in the stability
  selection. Automatically set to 20 subsamples with size
  \`10\*sqrt(n)\` if \`n \>= 144\` and \`0.8\*n\` otherwise following
  Liu et al. (2010) recommendations.

- `n_subsamples`:

  number of subsamples to create if the subsamples are not given

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
#> [1] 0.001058849
```
