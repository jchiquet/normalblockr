# R6 class for a collection of normal-block models with different number of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockUnknownQ`

## Active bindings

- `q_list`:

  number of blocks

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockUnknownQ$new()`](#method-NormalBlockUnknownQ-initialize)

- [`NormalBlockUnknownQ$get_model()`](#method-NormalBlockUnknownQ-get_model)

- [`NormalBlockUnknownQ$get_best_model()`](#method-NormalBlockUnknownQ-get_best_model)

- [`NormalBlockUnknownQ$plot()`](#method-NormalBlockUnknownQ-plot)

- [`NormalBlockUnknownQ$clone()`](#method-NormalBlockUnknownQ-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](NormalBlockCollection.html#method-optimize)

------------------------------------------------------------------------

### `NormalBlockUnknownQ$new()`

Create a new \[\`NormalBlockUnknownQ\`\] object.

#### Usage

    NormalBlockUnknownQ$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

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

A new \[\`NormalBlockUnknownQ\`\] object

------------------------------------------------------------------------

### `NormalBlockUnknownQ$get_model()`

returns the NormalBlockUnknownClusters model corresponding to given q

#### Usage

    NormalBlockUnknownQ$get_model(q)

#### Arguments

- `q`:

  number of blocks asked by user

#### Returns

A NormalBlockUnknownClusters object with given value q

------------------------------------------------------------------------

### `NormalBlockUnknownQ$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockUnknownQ$get_best_model(crit = c("ICL", "BIC", "EBIC", "deviance"))

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "ICL" or "BIC". "ICL" is the default criterion

#### Returns

a \[\`NormalBlockUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockUnknownQ$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NormalBlockUnknownQ$plot(criteria = c("deviance", "ICL", "BIC", "EBIC"))

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in \`c("deviance", "BIC",
  "ICL")\`. Defaults to all of them.

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### `NormalBlockUnknownQ$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockUnknownQ$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
