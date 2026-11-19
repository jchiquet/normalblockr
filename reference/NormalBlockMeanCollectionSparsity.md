# Collection of Mean-Block Models over a Sparsity Path

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
with a fixed clustering (or a fixed number of blocks) and different
sparsity levels applied to the p x p precision matrix of the variables.
Mirrors \[NormalBlockVarCollectionSparsity\], StARS included.

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionSparsity`](NormalBlockCollectionSparsity.md) -\>
`NormalBlockMeanCollectionSparsity`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanCollectionSparsity$new()`](#method-NormalBlockMeanCollectionSparsity-initialize)

- [`NormalBlockMeanCollectionSparsity$clone()`](#method-NormalBlockMeanCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_best_model()`](NormalBlockCollectionSparsity.html#method-get_best_model)
- [`NormalBlockCollectionSparsity$get_model()`](NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$optimize()`](NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](NormalBlockCollectionSparsity.html#method-plot)
- [`NormalBlockCollectionSparsity$stability_selection()`](NormalBlockCollectionSparsity.html#method-stability_selection)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$new()`

Create a new \[\`NormalBlockMeanCollectionSparsity\`\] object.

#### Usage

    NormalBlockMeanCollectionSparsity$new(mydata, blocks, control = NB_control())

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `blocks`:

  either a clustering matrix (known, fixed clustering) or a single
  integer (number of blocks to infer)

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockMeanCollectionSparsity\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanCollectionSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 3, sparsity = TRUE, model = "mean",
                       control = NB_control(n_sparsity_penalties = 5))
#> Fitting a normal-block-mean model with sparsity path 
#>   penalty = 2.722819           
     penalty = 0.861031           
     penalty = 0.2722819           
     penalty = 0.0861031           
     penalty = 0.02722819           

#> DONE
models$plot(c("BIC", "EBIC"))
```
