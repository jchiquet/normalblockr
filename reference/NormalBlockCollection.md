# R6 abstract class for a collection of normal-block models

Shared scaffolding for the collections explored by \[get_model()\]/
\[normal_block()\]: a sweep over sparsity penalties
(\[\`NormalBlockChangingSparsity\`\]), over the number of clusters
(\[\`NormalBlockUnknownQ\`\]), or over both
(\[\`NormalBlockUnknownQChangingSparsity\`\]). Concrete subclasses set
\`private\$progress_field\`/\`private\$progress_label\` in their
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

- [`NormalBlockCollection$optimize()`](#method-NormalBlockCollection-optimize)

- [`NormalBlockCollection$clone()`](#method-NormalBlockCollection-clone)

------------------------------------------------------------------------

### `NormalBlockCollection$optimize()`

optimizes every model (or sub-collection) in the collection

#### Usage

    NormalBlockCollection$optimize(
      control = list(niter = 100, threshold = 1e-04, verbose = TRUE)
    )

#### Arguments

- `control`:

  optimization parameters (niter and threshold)

------------------------------------------------------------------------

### `NormalBlockCollection$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockCollection$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
