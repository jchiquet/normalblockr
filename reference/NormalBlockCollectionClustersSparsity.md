# Base Class for a Collection over Cluster Counts and Sparsity Levels

Shared scaffolding for \[NormalBlockVarCollectionClustersSparsity\] and
\[NormalBlockMeanCollectionClustersSparsity\]: a collection of sparsity
sub-collections, one per q in increasing order. Everything lives here
(two-key model lookup, model selection over both axes, the criteria
heatmap); subclasses only name the model family.

## Super class

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\> `NormalBlockCollectionClustersSparsity`

## Active bindings

- `q_list`:

  number of blocks

- `sparsity`:

  list of penalties used for each q

## Methods

### Public methods

- [`NormalBlockCollectionClustersSparsity$new()`](#method-NormalBlockCollectionClustersSparsity-initialize)

- [`NormalBlockCollectionClustersSparsity$get_model()`](#method-NormalBlockCollectionClustersSparsity-get_model)

- [`NormalBlockCollectionClustersSparsity$get_best_model()`](#method-NormalBlockCollectionClustersSparsity-get_best_model)

- [`NormalBlockCollectionClustersSparsity$plot()`](#method-NormalBlockCollectionClustersSparsity-plot)

- [`NormalBlockCollectionClustersSparsity$clone()`](#method-NormalBlockCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)

------------------------------------------------------------------------

### `NormalBlockCollectionClustersSparsity$new()`

Create a collection over a range of q, sorted increasingly, with a
sparsity path for each q.

#### Usage

    NormalBlockCollectionClustersSparsity$new(
      mydata,
      q_list,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `q_list`:

  values of q (number of clusters) in the collection

- `zero_inflation`:

  whether the models of the collection are zero-inflated

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new collection object

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
