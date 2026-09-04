# Base Class for a Collection of Models over a Range of Cluster Counts

Shared scaffolding for \[NormalBlockVarCollectionClusters\] and
\[NormalBlockMeanCollectionClusters\]: everything that does not depend
on the model family (model lookup, model selection, the criteria plot
and the split/merge \`refine()\` search). Concrete subclasses only build
\`self\$models\` in their \`initialize()\` and name themselves through
\`who_am_I\`.

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockCollectionClusters`

## Active bindings

- `q_list`:

  number of blocks

## Methods

### Public methods

- [`NormalBlockCollectionClusters$get_model()`](#method-NormalBlockCollectionClusters-get_model)

- [`NormalBlockCollectionClusters$get_best_model()`](#method-NormalBlockCollectionClusters-get_best_model)

- [`NormalBlockCollectionClusters$plot()`](#method-NormalBlockCollectionClusters-plot)

- [`NormalBlockCollectionClusters$optimize()`](#method-NormalBlockCollectionClusters-optimize)

- [`NormalBlockCollectionClusters$refine()`](#method-NormalBlockCollectionClusters-refine)

- [`NormalBlockCollectionClusters$clone()`](#method-NormalBlockCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$get_model()`

returns the unknown-clusters model corresponding to given q

#### Usage

    NormalBlockCollectionClusters$get_model(q)

#### Arguments

- `q`:

  number of blocks asked by user

#### Returns

A unknown-clusters object with given value q

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockCollectionClusters$get_best_model(
      crit = c("ICL", "BIC", "EBIC", "deviance")
    )

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "ICL" or "BIC". "ICL" is the default criterion

#### Returns

a \[\`unknown-clusters\`\] object

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NormalBlockCollectionClusters$plot(
      criteria = c("deviance", "ICL", "BIC", "EBIC")
    )

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in \`c("deviance", "BIC",
  "ICL")\`. Defaults to all of them.

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$optimize()`

optimizes every model in the collection, then – if \`control\$refine\`
is \`TRUE\` (see \[NB_control()\], default \`FALSE\`) – calls
\[refine()\] automatically.

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
(\`"split"\`), a short merge-and-reoptimize trial seeded from its
larger-q neighbor (\`"merge"\`), or both (the default); a candidate
replaces the original only if it strictly lowers the deviance, so this
can only improve (or leave unchanged) each model it touches. Only
contiguous q pairs (\`q\` and \`q -/+ 1\`, both present in the
collection) are refined. See
\`inst/methods_initialization_and_refine.md\` for the rationale and
empirical evidence.

#### Usage

    NormalBlockCollectionClusters$refine(
      trial_niter = 2,
      max_candidates = 30,
      directions = c("split", "merge"),
      verbose = self$control$verbose
    )

#### Arguments

- `trial_niter`:

  number of EM iterations used for the cheap trial candidates (passed to
  \`candidates_split()\`/\`candidates_merge()\`) before fully
  re-optimizing only the best one.

- `max_candidates`:

  passed to \`candidates_merge()\` (ignored for \`"split"\`, which is
  never combinatorial in q) – see its documentation.

- `directions`:

  which neighbor(s) to seed refinement candidates from: \`"split"\`
  (smaller-q neighbor), \`"merge"\` (larger-q neighbor), or both (the
  default).

- `verbose`:

  whether to print, for each q attempted, whether the candidate from
  that neighbor improved on it. Defaults to \`control\$verbose\` (the
  value set at construction, see \[NB_control()\]).

#### Returns

invisibly returns \`self\`; improved models replace the originals in
\`\$models\` in place.

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
