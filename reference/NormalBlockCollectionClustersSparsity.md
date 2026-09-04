# Base Class for a Collection over Cluster Counts and Sparsity Levels

Shared scaffolding for \[NormalBlockVarCollectionClustersSparsity\] and
\[NormalBlockMeanCollectionClustersSparsity\]: a collection of sparsity
sub-collections, one per q. Everything family-agnostic (two-key model
lookup, model selection over both axes, the criteria heatmap) lives
here; subclasses only build \`self\$models\` and name themselves.

## Super class

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
`NormalBlockCollectionClustersSparsity`

## Active bindings

- `q_list`:

  number of blocks

- `sparsity`:

  list of penalties used for each q

## Methods

### Public methods

- [`NormalBlockCollectionClustersSparsity$get_model()`](#method-NormalBlockCollectionClustersSparsity-get_model)

- [`NormalBlockCollectionClustersSparsity$get_best_model()`](#method-NormalBlockCollectionClustersSparsity-get_best_model)

- [`NormalBlockCollectionClustersSparsity$plot()`](#method-NormalBlockCollectionClustersSparsity-plot)

- [`NormalBlockCollectionClustersSparsity$clone()`](#method-NormalBlockCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionClustersSparsity$get_model()`

returns a collection of models corresponding to given q or one single
model if penalty is also given

#### Usage

    NormalBlockCollectionClustersSparsity$get_model(q, sparsity = NA)

#### Arguments

- `q`:

  number of blocks asked by user.

- `sparsity`:

  sparsity penalty asked by user

#### Returns

either a sparsity sub-collection or a single model object

------------------------------------------------------------------------

### `NormalBlockCollectionClustersSparsity$get_best_model()`

Extract best model in the collection

#### Usage

    NormalBlockCollectionClustersSparsity$get_best_model(
      crit = c("ICL", "BIC", "EBIC")
    )

#### Arguments

- `crit`:

  a character for the criterion used to performed the selection. Either
  "BIC", "EBIC" or "ICL". "ICL" is the default criterion

#### Returns

a single fitted model object

------------------------------------------------------------------------

### `NormalBlockCollectionClustersSparsity$plot()`

Display various outputs (goodness-of-fit criteria, robustness,
diagnostic) associated with a collection of network fits

#### Usage

    NormalBlockCollectionClustersSparsity$plot(
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

a \[\`ggplot2::ggplot\`\] heatmap

------------------------------------------------------------------------

### `NormalBlockCollectionClustersSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockCollectionClustersSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- see
# normal_block() for how collections are created and fitted.
```
