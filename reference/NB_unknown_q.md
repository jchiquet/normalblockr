# R6 class for a collection of normal-block models with different number of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

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

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NB_unknown_q$new()`](#method-NB_unknown_q-new)

- [`NB_unknown_q$optimize()`](#method-NB_unknown_q-optimize)

- [`NB_unknown_q$get_model()`](#method-NB_unknown_q-get_model)

- [`NB_unknown_q$get_best_model()`](#method-NB_unknown_q-get_best_model)

- [`NB_unknown_q$plot()`](#method-NB_unknown_q-plot)

- [`NB_unknown_q$clone()`](#method-NB_unknown_q-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`NB_unknown_q\`\] object.

#### Usage

    NB_unknown_q$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NB_data class, with responses and design matrix

- `q_list`:

  list of q values (number of groups) in the collection

- `zero_inflation`:

  whether the models in the collection should be zero-inflated or not

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NB_unknown_q\`\] object

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

optimizes an NB_fixed_q object for each value of q

#### Usage

    NB_unknown_q$optimize(
      control = list(niter = 100, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### Method [`get_model()`](get_model.md)

returns the NB_fixed_q model corresponding to given q

#### Usage

    NB_unknown_q$get_model(q)

#### Arguments

- `q`:

  number of blocks asked by user

#### Returns

A NB_fixed_q object with given value q

------------------------------------------------------------------------

### Method `get_best_model()`

Extract best model in the collection

#### Usage

    NB_unknown_q$get_best_model(crit = c("ICL", "BIC", "EBIC", "deviance"))

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "ICL" or "BIC". "ICL" is the default criterion

#### Returns

a \[\`NB_fixed_q\`\] object

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NB_unknown_q$plot(criteria = c("deviance", "ICL", "BIC", "EBIC"))

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in \`c("deviance", "BIC",
  "ICL")\`. Defaults to all of them.

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NB_unknown_q$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
