# R6 class for normal-block model with unknown q (number of groups)

R6 class for normal-block model with unknown q (number of groups)

R6 class for normal-block model with unknown q (number of groups)

## Public fields

- `models`:

  list of NB_fixed_q models corresponding to each nb_block value

- `control`:

  store the list of user-defined model settings and optimization
  parameters

- `data`:

  object of NB_data class, with responses and design matrix

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

- [`NB_changing_sparsity$new()`](#method-NB_changing_sparsity-new)

- [`NB_changing_sparsity$optimize()`](#method-NB_changing_sparsity-optimize)

- [`NB_changing_sparsity$get_model()`](#method-NB_changing_sparsity-get_model)

- [`NB_changing_sparsity$get_best_model()`](#method-NB_changing_sparsity-get_best_model)

- [`NB_changing_sparsity$plot()`](#method-NB_changing_sparsity-plot)

- [`NB_changing_sparsity$stability_selection()`](#method-NB_changing_sparsity-stability_selection)

- [`NB_changing_sparsity$clone()`](#method-NB_changing_sparsity-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`NB_changing_sparsity\`\] object.

#### Usage

    NB_changing_sparsity$new(
      mydata,
      blocks,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NB_data class, with responses and design matrix

- `blocks`:

  group matrix or number of blocks.

- `zero_inflation`:

  boolean to specify whether data is zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NB_changing_sparsity\`\] object

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

optimizes a model for each penalty

#### Usage

    NB_changing_sparsity$optimize(
      control = list(niter = 100, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### Method [`get_model()`](get_model.md)

returns the NB_fixed_block model corresponding to given penalty

#### Usage

    NB_changing_sparsity$get_model(sparsity)

#### Arguments

- `sparsity`:

  sparsity penalty asked by user

#### Returns

A NB_fixed_blocks_sparse object with given value penalty

------------------------------------------------------------------------

### Method `get_best_model()`

Extract best model in the collection

#### Usage

    NB_changing_sparsity$get_best_model(
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

a \[\`NB_fixed_q\`\] object

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NB_changing_sparsity$plot(
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

### Method `stability_selection()`

Compute the stability path by stability selection

#### Usage

    NB_changing_sparsity$stability_selection(subsamples = NULL, n_subsamples = 10)

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NB_changing_sparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
