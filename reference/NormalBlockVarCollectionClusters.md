# Collection of Normal-Block Models over a Range of Cluster Counts

R6 class for a collection of normal-block models with different number
of clusters (q) and a fixed sparsity level.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionClusters`](NormalBlockCollectionClusters.md) -\>
`NormalBlockVarCollectionClusters`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarCollectionClusters$new()`](#method-NormalBlockVarCollectionClusters-initialize)

- [`NormalBlockVarCollectionClusters$clone()`](#method-NormalBlockVarCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClusters$get_best_model()`](NormalBlockCollectionClusters.html#method-get_best_model)
- [`NormalBlockCollectionClusters$get_model()`](NormalBlockCollectionClusters.html#method-get_model)
- [`NormalBlockCollectionClusters$optimize()`](NormalBlockCollectionClusters.html#method-optimize)
- [`NormalBlockCollectionClusters$plot()`](NormalBlockCollectionClusters.html#method-plot)
- [`NormalBlockCollectionClusters$refine()`](NormalBlockCollectionClusters.html#method-refine)

------------------------------------------------------------------------

### `NormalBlockVarCollectionClusters$new()`

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

### `NormalBlockVarCollectionClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:4, control = NB_control(verbose = FALSE))
models$get_best_model("ICL")$q
#> [1] 4
```
