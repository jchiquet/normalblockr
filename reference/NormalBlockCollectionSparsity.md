# Base Class for a Collection of Models over a Sparsity Path

Shared scaffolding for
[NormalBlockVarCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionSparsity.md)
and
[NormalBlockMeanCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionSparsity.md):
the penalty grid, the warm-started path traversal, the penalty lookup
and the criteria plot. Model selection, StARS included (stability
selection, Liu et al. 2010), is shared: a family only says on which
covariance estimate the default penalty grid is read
(`private$penalty_scale()`), how to refit its path on a subsample with
the clustering held fixed (`private$subsample_path()`) and how to name
the nodes of its network (`private$network_nodes()`).

## Super class

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\> `NormalBlockCollectionSparsity`

## Public fields

- `data`:

  object of NormalBlockData class, with responses and design matrix

## Active bindings

- `sparsity_details`:

  number, range and ratio of the penalties of the path

- `q`:

  number of blocks

- `blocks`:

  group matrix or number of blocks

- `sparsity`:

  list of sparsity penalties

- `criteria`:

  a data frame with the criteria of every model of the path, plus its
  StARS stability once `stability_selection()` has run

- `stability_path`:

  selection probability of every edge at every penalty (StARS), `NULL`
  until `stability_selection()` has run

- `stability`:

  mean edge stability along the penalty path, `1 - mean(4 p (1 - p))`
  over the edges; `NA` until `stability_selection()` has run

## Methods

### Public methods

- [`NormalBlockCollectionSparsity$new()`](#method-NormalBlockCollectionSparsity-initialize)

- [`NormalBlockCollectionSparsity$optimize()`](#method-NormalBlockCollectionSparsity-optimize)

- [`NormalBlockCollectionSparsity$get_model()`](#method-NormalBlockCollectionSparsity-get_model)

- [`NormalBlockCollectionSparsity$get_best_model()`](#method-NormalBlockCollectionSparsity-get_best_model)

- [`NormalBlockCollectionSparsity$stability_selection()`](#method-NormalBlockCollectionSparsity-stability_selection)

- [`NormalBlockCollectionSparsity$plot()`](#method-NormalBlockCollectionSparsity-plot)

- [`NormalBlockCollectionSparsity$clone()`](#method-NormalBlockCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$new()`

Create a collection over a path of sparsity penalties. By default,
`control$n_sparsity_penalties` values on a log scale, from the largest
absolute (weighted) entry of a covariance estimate, which removes every
edge, down to `control$min_ratio` times it.

#### Usage

    NormalBlockCollectionSparsity$new(
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

  whether the models of the path are zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new collection object

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$optimize()`

optimizes every model in the sparsity path, warm-starting each one
(after the first) from the previous, adjacent penalty's converged
parameters (see the family's base class's `warm_start_from()`) instead
of re-deriving everything from the heuristic clustering, the way the
generic
[NormalBlockCollection](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
[`optimize()`](https://rdrr.io/r/stats/optimize.html) would. `blocks`
(hence q) is fixed across the whole path, only the sparsity penalty
changes, so the warm start is always between models of matching shape.

#### Usage

    NormalBlockCollectionSparsity$optimize(
      control = list(niter = 500, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$get_model()`

returns the NormalBlockVarKnownClusters model corresponding to given
penalty

#### Usage

    NormalBlockCollectionSparsity$get_model(sparsity)

#### Arguments

- `sparsity`:

  sparsity penalty asked by user

#### Returns

A NormalBlockVarKnownClusters (sparse) object with given value penalty

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockCollectionSparsity$get_best_model(
      crit = c("BIC", "EBIC", "ICL", "StARS"),
      stability = 0.9
    )

#### Arguments

- `crit`:

  a character for the criterion used to perform the selection: "BIC"
  (default), "EBIC", "ICL" or "StARS".

- `stability`:

  if `crit = "StARS"`, the level of stability required.

#### Returns

a model of the collection's family

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$stability_selection()`

Computes the stability path (StARS, Liu et al. 2010): the penalty path
is refitted on subsamples of the observations, with the clustering held
fixed, and the selection frequency of each edge is recorded at every
penalty.

#### Usage

    NormalBlockCollectionSparsity$stability_selection(
      subsamples = NULL,
      n_subsamples = 10
    )

#### Arguments

- `subsamples`:

  a list of vectors describing the subsamples. The number of vectors (or
  list length) determines the number of subsamples used in the stability
  selection. Automatically set to 10 subsamples with size `10*sqrt(n)`
  if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010)
  recommendations.

- `n_subsamples`:

  number of subsamples to create if the subsamples are not given

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$plot()`

Plot the criteria of the models of the path against the penalty, with a
dashed line at the penalty selected by each criterion (see also
[`plot.NormalBlockCollection()`](https://jchiquet.github.io/normalblockr/reference/plot.NormalBlockCollection.md)).

#### Usage

    NormalBlockCollectionSparsity$plot(
      criteria = c("deviance", "BIC", "EBIC", "ICL"),
      log.x = TRUE
    )

#### Arguments

- `criteria`:

  the criteria to plot, among `"deviance"`, `"BIC"`, `"EBIC"` and
  `"ICL"`; all of them by default.

- `log.x`:

  logical: should the x-axis be represented in log-scale? Default is
  `TRUE`.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockCollectionSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- see
# normal_block() for how collections are created and fitted.
```
