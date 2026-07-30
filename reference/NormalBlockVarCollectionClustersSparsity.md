# Collection of Normal-Block Models over Cluster Counts and Sparsity Levels

R6 class for a collection of normal-block models with different number
of clusters (q) and different sparsity levels.

## Super class

[`NormalBlockVarCollection`](NormalBlockVarCollection.md) -\>
`NormalBlockVarCollectionClustersSparsity`

## Active bindings

- `q_list`:

  number of blocks

- `sparsity`:

  list of penalties used for each q

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionClustersSparsity$new()`](#method-NormalBlockVarCollectionClustersSparsity-initialize)

- [`NormalBlockVarCollectionClustersSparsity$get_model()`](#method-NormalBlockVarCollectionClustersSparsity-get_model)

- [`NormalBlockVarCollectionClustersSparsity$get_best_model()`](#method-NormalBlockVarCollectionClustersSparsity-get_best_model)

- [`NormalBlockVarCollectionClustersSparsity$plot()`](#method-NormalBlockVarCollectionClustersSparsity-plot)

- [`NormalBlockVarCollectionClustersSparsity$clone()`](#method-NormalBlockVarCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockVarCollection$optimize()`](NormalBlockVarCollection.html#method-optimize)
- [`NormalBlockVarCollection$print()`](NormalBlockVarCollection.html#method-print)
- [`NormalBlockVarCollection$summary()`](NormalBlockVarCollection.html#method-summary)

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

### `NormalBlockVarCollectionClustersSparsity$get_model()`

returns a collection of models corresponding to given q or one single
model if penalty is also given

#### Usage

    NormalBlockVarCollectionClustersSparsity$get_model(q, sparsity = NA)

#### Arguments

- `q`:

  number of blocks asked by user.

- `sparsity`:

  sparsity penalty asked by user

#### Returns

either a NormalBlockVarCollectionSparsity or a
NormalBlockVarUnknownClusters object

------------------------------------------------------------------------

### `NormalBlockVarCollectionClustersSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockVarCollectionClustersSparsity$get_best_model(
      crit = c("ICL", "BIC", "EBIC")
    )

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "BIC", "EBIC" or "ICL". "ICL" is the default criterion

#### Returns

a \[\`NormalBlockVarUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockVarCollectionClustersSparsity$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NormalBlockVarCollectionClustersSparsity$plot(
      criterion = c("deviance", "ICL", "BIC", "EBIC"),
      n_intervals = NULL
    )

#### Arguments

- `criterion`:

  The criteria to plot in \`c("deviance", BIC", "EBIC", "ICL")\`.
  Defaults deviance.

- `n_intervals`:

  number of intervals into which the penalties range should be split

#### Returns

a \[\`ggplot\`\] heatmap

------------------------------------------------------------------------

### `NormalBlockVarCollectionClustersSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionClustersSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
