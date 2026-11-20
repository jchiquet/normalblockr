# Normal-Block-Mean Model with Unknown Clustering

R6 class for a Normal-Block-Mean model with a fixed number of clusters
(but unknown clustering), inferred by variational EM.

## Super classes

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\>
[`NormalBlockMeanBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.md)
-\> `NormalBlockMeanUnknownClusters`

## Active bindings

- `model_par`:

  the parameters of \[NormalBlockMeanBase\] (B, Omega), plus alpha
  (cluster proportions)

- `var_par`:

  a list with the variational parameter: tau (posterior group
  probabilities)

## Methods

### Public methods

- [`NormalBlockMeanUnknownClusters$new()`](#method-NormalBlockMeanUnknownClusters-initialize)

- [`NormalBlockMeanUnknownClusters$clone()`](#method-NormalBlockMeanUnknownClusters-clone)

Inherited methods

- [`NormalBlockBase$best_of_inits()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-best_of_inits)
- [`NormalBlockBase$candidates_merge()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-candidates_merge)
- [`NormalBlockBase$candidates_split()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-candidates_split)
- [`NormalBlockBase$heldout_loglik()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-heldout_loglik)
- [`NormalBlockBase$latent_network()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-latent_network)
- [`NormalBlockBase$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-optimize)
- [`NormalBlockBase$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-plot)
- [`NormalBlockBase$plot_loglik()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-plot_loglik)
- [`NormalBlockBase$plot_network()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-plot_network)
- [`NormalBlockBase$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-update)
- [`NormalBlockBase$warm_start_from()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-warm_start_from)
- [`NormalBlockMeanBase$merge()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.html#method-merge)
- [`NormalBlockMeanBase$predict()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.html#method-predict)
- [`NormalBlockMeanBase$split()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.html#method-split)

------------------------------------------------------------------------

### `NormalBlockMeanUnknownClusters$new()`

Create a new \[\`NormalBlockMeanUnknownClusters\`\] object.

#### Usage

    NormalBlockMeanUnknownClusters$new(
      data,
      q,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `q`:

  number of clusters

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockMeanUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanUnknownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = 3, model = "mean")
#> Fitting a diagonal normal-block-mean model with 3 unknown blocks 
#> 
#> DONE
model$clustering
#>  [1] 1 2 1 2 2 2 2 2 1 1 2 2 2 2 2 1 1 2 2 1
```
