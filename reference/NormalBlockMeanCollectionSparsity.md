# Collection of Normal-Block-Mean Models over a Sparsity Path

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
with a fixed clustering (or a fixed number of blocks) and different
sparsity levels applied to the p x p precision matrix of the variables.
Mirrors \[NormalBlockVarCollectionSparsity\], StARS included; the
default penalty grid is read off the covariance of the least-squares
residuals.

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.md)
-\> `NormalBlockMeanCollectionSparsity`

## Methods

### Public methods

- [`NormalBlockMeanCollectionSparsity$new()`](#method-NormalBlockMeanCollectionSparsity-initialize)

- [`NormalBlockMeanCollectionSparsity$clone()`](#method-NormalBlockMeanCollectionSparsity-clone)

Inherited methods

- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionSparsity$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-get_best_model)
- [`NormalBlockCollectionSparsity$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-get_model)
- [`NormalBlockCollectionSparsity$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-optimize)
- [`NormalBlockCollectionSparsity$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-plot)
- [`NormalBlockCollectionSparsity$stability_selection()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.html#method-stability_selection)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionSparsity$new()`

Create a new \[\`NormalBlockMeanCollectionSparsity\`\] object.

#### Usage

    NormalBlockMeanCollectionSparsity$new(
      mydata,
      blocks,
      zero_inflation = FALSE,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `blocks`:

  either a clustering matrix (known, fixed clustering) or a single
  integer (number of blocks to infer)

- `zero_inflation`:

  must be \`FALSE\`: the path penalizes a full Sigma, which
  zero-inflated mean-block models do not carry

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
#> Fitting a collection of full normal-block-mean models over sparsity penalties, with fixed q 
#>   penalty = 5.060078           
     penalty = 1.600137           
     penalty = 0.5060078           
     penalty = 0.1600137           
     penalty = 0.05060078           

#> DONE
models$plot(c("BIC", "EBIC"))
```
