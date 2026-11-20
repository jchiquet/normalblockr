# Zero-Inflated Normal-Block-Mean Model with Known Clustering

R6 class for a zero-inflated Normal-Block-Mean model with a known
clustering. Sigma is diagonal or spherical here. See
\[NormalBlockMeanBase\] for why a full one is out of reach under a mask.

## Super classes

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\>
[`NormalBlockMeanBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.md)
-\> `ZINormalBlockMeanKnownClusters`

## Active bindings

- `model_par`:

  the parameters of \[NormalBlockMeanBase\] (B, Omega), plus B0
  (zero-inflation coefficients) and kappa (zero-inflation probabilities)

## Methods

### Public methods

- [`ZINormalBlockMeanKnownClusters$new()`](#method-ZINormalBlockMeanKnownClusters-initialize)

- [`ZINormalBlockMeanKnownClusters$clone()`](#method-ZINormalBlockMeanKnownClusters-clone)

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

### `ZINormalBlockMeanKnownClusters$new()`

Create a new \[\`ZINormalBlockMeanKnownClusters\`\] object.

#### Usage

    ZINormalBlockMeanKnownClusters$new(
      data,
      C,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `C`:

  p x q clustering matrix, C_jk = 1 if variable j belongs to cluster k

- `sparsity`:

  unused here, kept for signature symmetry (must be 0)

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`ZINormalBlockMeanKnownClusters\`\] object

------------------------------------------------------------------------

### `ZINormalBlockMeanKnownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINormalBlockMeanKnownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
Y <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
data <- NormalBlockData$new(Y, ex$X)
model <- normal_block(data, blocks = ex$parameters$C, model = "mean",
                      zero_inflation = TRUE)
#> Fitting a zero-inflated diagonal normal-block-mean model with fixed blocks 
#> 
#> DONE
model$clustering
#>  [1] 3 3 2 2 1 2 2 2 2 2 1 2 2 2 3 2 2 2 1 2
```
