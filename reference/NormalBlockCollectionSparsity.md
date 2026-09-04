# Base Class for a Collection of Models over a Sparsity Path

Shared scaffolding for \[NormalBlockVarCollectionSparsity\] and
\[NormalBlockMeanCollectionSparsity\]: the warm-started path traversal,
the penalty lookup and the criteria plot. Concrete subclasses derive the
penalty grid in their \`initialize()\` and provide their own
\`get_best_model()\` (only the variance-block family offers StARS).

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockCollectionSparsity`

## Public fields

- `data`:

  object of NormalBlockData class, with responses and design matrix

## Active bindings

- `q`:

  number of blocks

- `blocks`:

  group matrix or number of blocks

- `sparsity`:

  list of sparsity penalties

## Methods

### Public methods

- [`NormalBlockCollectionSparsity$optimize()`](#method-NormalBlockCollectionSparsity-optimize)

- [`NormalBlockCollectionSparsity$get_model()`](#method-NormalBlockCollectionSparsity-get_model)

- [`NormalBlockCollectionSparsity$plot()`](#method-NormalBlockCollectionSparsity-plot)

- [`NormalBlockCollectionSparsity$clone()`](#method-NormalBlockCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionSparsity$optimize()`

optimizes every model in the sparsity path, warm-starting each one
(after the first) from the previous, adjacent penalty's converged
parameters (see the family's base class's \`warm_start_from()\`) instead
of re-deriving everything from the heuristic clustering, the way the
generic \[NormalBlockCollection\] \`optimize()\` would. \`blocks\`
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

### `NormalBlockCollectionSparsity$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a collection)

#### Usage

    NormalBlockCollectionSparsity$plot(
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

a \[\`ggplot2::ggplot\`\] graph

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
