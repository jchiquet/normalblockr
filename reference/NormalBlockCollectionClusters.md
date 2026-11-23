# Base Class for a Collection of Models over a Range of Cluster Counts

Shared scaffolding for
[NormalBlockVarCollectionClusters](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionClusters.md)
and
[NormalBlockMeanCollectionClusters](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionClusters.md):
the models of the collection, one per q in increasing order, model
lookup and selection, the criteria plot and the split/merge `refine()`
search. Concrete subclasses only name the model family.

## Super class

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\> `NormalBlockCollectionClusters`

## Active bindings

- `q_list`:

  number of blocks

- `deviance_violations`:

  the values of q whose deviance is higher than at the previous q of the
  collection. The models being nested, each one signals a poor local
  optimum, which `refine()` can repair.

## Methods

### Public methods

- [`NormalBlockCollectionClusters$new()`](#method-NormalBlockCollectionClusters-initialize)

- [`NormalBlockCollectionClusters$get_model()`](#method-NormalBlockCollectionClusters-get_model)

- [`NormalBlockCollectionClusters$get_best_model()`](#method-NormalBlockCollectionClusters-get_best_model)

- [`NormalBlockCollectionClusters$plot()`](#method-NormalBlockCollectionClusters-plot)

- [`NormalBlockCollectionClusters$optimize()`](#method-NormalBlockCollectionClusters-optimize)

- [`NormalBlockCollectionClusters$refine()`](#method-NormalBlockCollectionClusters-refine)

- [`NormalBlockCollectionClusters$print()`](#method-NormalBlockCollectionClusters-print)

- [`NormalBlockCollectionClusters$clone()`](#method-NormalBlockCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$new()`

Create a collection over a range of q, sorted increasingly.

#### Usage

    NormalBlockCollectionClusters$new(
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

  values of q (number of clusters) in the collection

- `zero_inflation`:

  whether the models of the collection are zero-inflated

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new collection object

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$get_model()`

The model of the collection with q clusters.

#### Usage

    NormalBlockCollectionClusters$get_model(q)

#### Arguments

- `q`:

  a number of clusters among those of the collection (`$q_list`)

#### Returns

a fitted model, of the class of the collection's models

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$get_best_model()`

The best model of the collection for a criterion.

#### Usage

    NormalBlockCollectionClusters$get_best_model(
      crit = c("ICL", "BIC", "EBIC", "deviance")
    )

#### Arguments

- `crit`:

  the criterion minimized: "ICL" (the default), "BIC", "EBIC", or
  "deviance" (which always selects the largest q)

#### Returns

a copy of the selected fitted model

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$plot()`

Plot the criteria of the models of the collection against q, with a
dashed line at the q selected by each criterion (see also
[`plot.NormalBlockCollection()`](https://jchiquet.github.io/normalblockr/reference/plot.NormalBlockCollection.md)).

#### Usage

    NormalBlockCollectionClusters$plot(
      criteria = c("deviance", "ICL", "BIC", "EBIC")
    )

#### Arguments

- `criteria`:

  the criteria to plot, among `"deviance"`, `"ICL"`, `"BIC"` and
  `"EBIC"`; all of them by default.

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$optimize()`

optimizes every model in the collection, then – if `control$refine` is
`TRUE` (see
[`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md),
the default) – calls `$refine()` automatically, reporting in one line
how many models it improved.

#### Usage

    NormalBlockCollectionClusters$optimize(
      control = list(niter = 500, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter, threshold, verbose)

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$refine()`

Tries to improve every model in the collection with a short
split-and-reoptimize trial seeded from its smaller-q neighbor
(`"split"`), a short merge-and-reoptimize trial seeded from its larger-q
neighbor (`"merge"`), or both (the default); a candidate replaces the
original only if it strictly lowers the deviance, so this can only
improve (or leave unchanged) each model it touches. Only contiguous q
pairs (`q` and `q -/+ 1`, both present in the collection) are refined.
Run by default after fitting a collection, see
[`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)'s
`refine` for the rationale.

#### Usage

    NormalBlockCollectionClusters$refine(
      trial_niter = 2,
      max_candidates = 30,
      directions = c("split", "merge"),
      verbose = self$control$verbose,
      split_rule = NULL
    )

#### Arguments

- `trial_niter`:

  number of EM iterations used for the cheap trial candidates (passed to
  `candidates_split()`/`candidates_merge()`) before fully re-optimizing
  only the best one.

- `max_candidates`:

  passed to `candidates_merge()` (ignored for `"split"`, which is never
  combinatorial in q) – see its documentation.

- `directions`:

  which neighbor(s) to seed refinement candidates from: `"split"`
  (smaller-q neighbor), `"merge"` (larger-q neighbor), or both (the
  default).

- `verbose`:

  whether to print, for each q attempted, whether the candidate from
  that neighbor improved on it. Defaults to `control$verbose` (the value
  set at construction, see
  [`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)).

- `split_rule`:

  how a cluster is split to seed the q + 1 model, see the `rule`
  argument of the models'
  [`split()`](https://rdrr.io/r/base/split.html): `"kmeans"`, `"noise"`,
  or `NULL` (default) for the family's own default (`"kmeans"` for the
  variance-block family, `"noise"` for the mean-block one).

#### Returns

invisibly returns `self`; improved models replace the originals in
`$models` in place.

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$print()`

User-friendly print method: as for any collection, plus a warning line
when the deviance increases with q somewhere (see
`deviance_violations`).

#### Usage

    NormalBlockCollectionClusters$print()

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- see
# normal_block() for how collections are created and fitted.
```
