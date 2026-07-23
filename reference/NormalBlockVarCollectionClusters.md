# R6 class for a collection of normal-block models with different number of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

## Super class

[`normalblockr::NormalBlockVarCollection`](NormalBlockVarCollection.md)
-\> `NormalBlockVarCollectionClusters`

## Active bindings

- `q_list`:

  number of blocks

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionClusters$new()`](#method-NormalBlockVarCollectionClusters-new)

- [`NormalBlockVarCollectionClusters$get_model()`](#method-NormalBlockVarCollectionClusters-get_model)

- [`NormalBlockVarCollectionClusters$get_best_model()`](#method-NormalBlockVarCollectionClusters-get_best_model)

- [`NormalBlockVarCollectionClusters$plot()`](#method-NormalBlockVarCollectionClusters-plot)

- [`NormalBlockVarCollectionClusters$optimize()`](#method-NormalBlockVarCollectionClusters-optimize)

- [`NormalBlockVarCollectionClusters$refine()`](#method-NormalBlockVarCollectionClusters-refine)

- [`NormalBlockVarCollectionClusters$clone()`](#method-NormalBlockVarCollectionClusters-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`NormalBlockVarCollectionClusters\`\] object.

#### Usage

    NormalBlockVarCollectionClusters$new(
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

A new \[\`NormalBlockVarCollectionClusters\`\] object

------------------------------------------------------------------------

### Method [`get_model()`](get_model.md)

returns the NormalBlockVarUnknownClusters model corresponding to given q

#### Usage

    NormalBlockVarCollectionClusters$get_model(q)

#### Arguments

- `q`:

  number of blocks asked by user

#### Returns

A NormalBlockVarUnknownClusters object with given value q

------------------------------------------------------------------------

### Method `get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockVarCollectionClusters$get_best_model(
      crit = c("ICL", "BIC", "EBIC", "deviance")
    )

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "ICL" or "BIC". "ICL" is the default criterion

#### Returns

a \[\`NormalBlockVarUnknownClusters\`\] object

------------------------------------------------------------------------

### Method [`plot()`](https://rdrr.io/r/graphics/plot.default.html)

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits (a
\[\`Networkfamily\`\])

#### Usage

    NormalBlockVarCollectionClusters$plot(
      criteria = c("deviance", "ICL", "BIC", "EBIC")
    )

#### Arguments

- `criteria`:

  vector of characters. The criteria to plot in \`c("deviance", "BIC",
  "ICL")\`. Defaults to all of them.

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

optimizes every model in the collection, then – if \`control\$refine\`
is \`TRUE\` (see \[NB_control()\], default \`FALSE\`) – calls
\[refine()\] automatically.

#### Usage

    NormalBlockVarCollectionClusters$optimize(
      control = list(niter = 500, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter, threshold, verbose)

------------------------------------------------------------------------

### Method `refine()`

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

    NormalBlockVarCollectionClusters$refine(
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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
