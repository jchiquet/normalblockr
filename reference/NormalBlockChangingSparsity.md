# R6 class for a collection of normal-block models with a fixed clustering (blocks) and different sparsity levels.

R6 class for a collection of normal-block models with a fixed clustering
(blocks) and different sparsity levels.

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockChangingSparsity`

## Public fields

- `data`:

  object of NormalBlockData class, with responses and design matrix

## Active bindings

- `q`:

  number of blocks

- `blocks`:

  group matrix or number of blocks.

- `sparsity`:

  list of sparsity penalties

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

- [`NormalBlockChangingSparsity$new()`](#method-NormalBlockChangingSparsity-initialize)

- [`NormalBlockChangingSparsity$optimize()`](#method-NormalBlockChangingSparsity-optimize)

- [`NormalBlockChangingSparsity$get_model()`](#method-NormalBlockChangingSparsity-get_model)

- [`NormalBlockChangingSparsity$get_best_model()`](#method-NormalBlockChangingSparsity-get_best_model)

- [`NormalBlockChangingSparsity$plot()`](#method-NormalBlockChangingSparsity-plot)

- [`NormalBlockChangingSparsity$stability_selection()`](#method-NormalBlockChangingSparsity-stability_selection)

- [`NormalBlockChangingSparsity$clone()`](#method-NormalBlockChangingSparsity-clone)

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$new()`

Create a new \[\`NormalBlockChangingSparsity\`\] object.

#### Usage

    NormalBlockChangingSparsity$new(
      mydata,
      blocks,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `blocks`:

  group matrix or number of blocks.

- `zero_inflation`:

  boolean to specify whether data is zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockChangingSparsity\`\] object

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$optimize()`

optimizes every model in the sparsity path, warm-starting each one
(after the first) from the previous, adjacent penalty's converged
parameters (see \[NormalBlockBase\]'s \`warm_start_from()\`) instead of
re-deriving everything from the heuristic clustering, the way the
generic \[NormalBlockCollection\] \`optimize()\` would. \`blocks\`
(hence q) is fixed across the whole path, only the sparsity penalty
changes, so the warm start is always between models of matching shape.

#### Usage

    NormalBlockChangingSparsity$optimize(
      control = list(niter = 100, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$get_model()`

returns the NormalBlockKnownClusters model corresponding to given
penalty

#### Usage

    NormalBlockChangingSparsity$get_model(sparsity)

#### Arguments

- `sparsity`:

  sparsity penalty asked by user

#### Returns

A NormalBlockKnownClusters (sparse) object with given value penalty

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockChangingSparsity$get_best_model(
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

a \[\`NormalBlockUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NormalBlockChangingSparsity$plot(
      criteria = c("deviance", "BIC", "EBIC", "ICL"),
      log.x = TRUE
    )

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in \`c("deviance", BIC",
  "EBIC", "ICL")\`. Defaults to all of them.

- `log.x`:

  logical: should the x-axis be represented in log-scale? Default is
  \`TRUE\`.

#### Returns

a \[\`ggplot\`\] graph

------------------------------------------------------------------------

### `NormalBlockChangingSparsity$stability_selection()`

Compute the stability path by stability selection

#### Usage

    NormalBlockChangingSparsity$stability_selection(
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

### `NormalBlockChangingSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockChangingSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
