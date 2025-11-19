# R6 class for a collection of normal-block models with different number of clusters (q) and different sparsity levels.

R6 class for a collection of normal-block models with different number
of clusters (q) and different sparsity levels.

R6 class for a collection of normal-block models with different number
of clusters (q) and different sparsity levels.

## Public fields

- `models`:

  list of NB_fixed_q models corresponding to each nb_block value

- `control`:

  store the list of user-defined model settings and optimization
  parameters

## Active bindings

- `q_list`:

  number of blocks

- `criteria`:

  a data frame with the values of some criteria ((approximated)
  log-likelihood, BIC, ICL) for the collection of models

- `sparsity`:

  list of penalties used for each q

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NB_unknown_q_changing_sparsity$new()`](#method-NB_unknown_q_changing_sparsity-new)

- [`NB_unknown_q_changing_sparsity$optimize()`](#method-NB_unknown_q_changing_sparsity-optimize)

- [`NB_unknown_q_changing_sparsity$get_model()`](#method-NB_unknown_q_changing_sparsity-get_model)

- [`NB_unknown_q_changing_sparsity$get_best_model()`](#method-NB_unknown_q_changing_sparsity-get_best_model)

- [`NB_unknown_q_changing_sparsity$plot()`](#method-NB_unknown_q_changing_sparsity-plot)

- [`NB_unknown_q_changing_sparsity$clone()`](#method-NB_unknown_q_changing_sparsity-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`NB_unknown_q_changing_sparsity\`\] object.

#### Usage

    NB_unknown_q_changing_sparsity$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NB_data class, with responses and design matrix

- `q_list`:

  list of q values (number of groups) in the collection

- `zero_inflation`:

  boolean to specify whether data is zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NB_unknown_q_changing_sparsity\`\] object

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

optimizes an NB_changing_sparsity object for each penalty value

#### Usage

    NB_unknown_q_changing_sparsity$optimize(
      control = list(niter = 100, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### Method [`get_model()`](get_model.md)

returns a collection of NB_unknown models corresponding to given q or
one single model if penalty is also given

#### Usage

    NB_unknown_q_changing_sparsity$get_model(q, sparsity = NA)

#### Arguments

- `q`:

  number of blocks asked by user.

- `sparsity`:

  sparsity penalty penalty asked by user

#### Returns

either a NB_changing_sparsity or a NB_fixed_q object

------------------------------------------------------------------------

### Method `get_best_model()`

Extract best model in the collection

#### Usage

    NB_unknown_q_changing_sparsity$get_best_model(crit = c("ICL", "BIC", "EBIC"))

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "BIC", "EBIC" or "ICL. "ICL" is the default criterion

#### Returns

a \[\`NB_unknown\`\] object

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NB_unknown_q_changing_sparsity$plot(
      criterion = c("deviance", "ICL", "BIC", "EBIC"),
      n_intervals = NULL
    )

#### Arguments

- `criterion`:

  The criteria to plot in \`c("deviance", BIC", "EBIC", "ICL")\`.
  Defaults deviance.

- `n_intervals`:

  number of intervals into which the penalties range should be splitted

#### Returns

a \[\`ggplot\`\] heatmap

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NB_unknown_q_changing_sparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
