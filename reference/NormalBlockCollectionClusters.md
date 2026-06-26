# R6 class for a collection of normal-block models with different number of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockCollectionClusters`

## Active bindings

- `q_list`:

  number of blocks

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockCollectionClusters$new()`](#method-NormalBlockCollectionClusters-initialize)

- [`NormalBlockCollectionClusters$get_model()`](#method-NormalBlockCollectionClusters-get_model)

- [`NormalBlockCollectionClusters$get_best_model()`](#method-NormalBlockCollectionClusters-get_best_model)

- [`NormalBlockCollectionClusters$plot()`](#method-NormalBlockCollectionClusters-plot)

- [`NormalBlockCollectionClusters$optimize()`](#method-NormalBlockCollectionClusters-optimize)

- [`NormalBlockCollectionClusters$refine()`](#method-NormalBlockCollectionClusters-refine)

- [`NormalBlockCollectionClusters$clone()`](#method-NormalBlockCollectionClusters-clone)

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$new()`

Create a new \[\`NormalBlockCollectionClusters\`\] object.

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

  list of q values (number of groups) in the collection

- `zero_inflation`:

  whether the models in the collection should be zero-inflated or not

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockCollectionClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockCollectionClusters$get_model()`

returns the NormalBlockUnknownClusters model corresponding to given q

#### Usage

    NormalBlockCollectionClusters$get_model(q)

#### Arguments

- `q`:

  number of blocks asked by user

#### Returns

A NormalBlockUnknownClusters object with given value q

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

a \[\`NormalBlockUnknownClusters\`\] object

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
replaces the original only if it strictly lowers the deviance. Useful
because every model in the collection is cold-started independently (a
heuristic clustering on the residuals, see
\`private\$clustering_methods\` in \[NormalBlockBase\]) and can settle
into a milder local optimum than an adjacent q's solution would, once
split/merged by one cluster – a gap that "deviance must not increase in
q" alone cannot reliably catch (see Details).

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

#### Details

Two cheaper attempts at targeting \*which\* q to refine – flagging
outright deviance increases, then flagging drops unusually small
relative to neighboring q's – found nothing to fix on a real dataset
(\`brca_rppa\`) where a full split/merge search was known to do
systematically better at almost every q: the gap was spread evenly
across the whole range, leaving no local anomaly to detect. \`refine()\`
instead tries unconditionally everywhere and discards whatever doesn't
help. On \`brca_rppa\`, \`"split"\` alone already recovers part of that
gap (~1.5-2x the cost of fitting the collection alone); adding
\`"merge"\` closes substantially more of it – the two directions catch
different cases, since split-from-below and merge-from-above are
different, independently cold-started starting points – for about twice
the cost again, still well under \[SelectionNClusters\]'s full
forward/backward search. Because nothing is kept unless strictly better,
this never makes a model worse: on a second dataset (\`university\`)
where independent cold starts were already better than chaining models
together, both directions found close to nothing to fix rather than
degrading anything.

Only contiguous q pairs (\`q\` and \`q -/+ 1\`, both present in the
collection) are refined; gaps in \`q_list\` are left untouched on that
side.

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
