# R6 abstract class for a collection of normal-block models

R6 abstract class for a collection of normal-block models

R6 abstract class for a collection of normal-block models

## Details

Shared scaffolding for the collections explored by \[get_model()\]/
\[normal_block()\]: a sweep over sparsity penalties
(\[\`NormalBlockVarCollectionSparsity\`\]), over the number of clusters
(\[\`NormalBlockVarCollectionClusters\`\]), or over both
(\[\`NormalBlockVarCollectionClustersSparsity\`\]). Concrete subclasses
set \`private\$progress_field\`/\`private\$progress_label\` in their
\`initialize()\` and provide their own \`get_best_model()\`, delegating
the (row of \`self\$criteria\` minimizing a criterion) lookup to
\`private\$best_id()\`.

## Public fields

- `models`:

  list of models (or sub-collections) explored by the collection

- `control`:

  store the list of user-defined model settings and optimization
  parameters

## Active bindings

- `criteria`:

  a data frame with the values of some criteria for the collection of
  models

## Methods

### Public methods

- [`NormalBlockVarCollection$optimize()`](#method-NormalBlockVarCollection-optimize)

- [`NormalBlockVarCollection$clone()`](#method-NormalBlockVarCollection-clone)

------------------------------------------------------------------------

### Method [`optimize()`](https://rdrr.io/r/stats/optimize.html)

optimizes every model (or sub-collection) in the collection

#### Usage

    NormalBlockVarCollection$optimize(
      control = list(niter = 500, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold). When
  \`control\$clustering_init\` is \`"best_of_inits"\`, each leaf model
  is fit via its own \`best_of_inits()\` instead of a plain
  \`optimize()\`.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
