# Zero-Inflated Normal-Block-Var Model with Known Clustering

R6 class for a zero-inflated normal-block model with a known clustering.

## Super classes

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\>
[`NormalBlockVarBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.md)
-\> `ZINormalBlockVarKnownClusters`

## Active bindings

- `posterior_par`:

  a list with the parameters of posterior distribution W \| Y

- `model_par`:

  the parameters of \[NormalBlockVarBase\] (B, B0, dm1, Omega), plus
  kappa (zero-inflation probabilities)

## Methods

### Public methods

- [`ZINormalBlockVarKnownClusters$new()`](#method-ZINormalBlockVarKnownClusters-initialize)

- [`ZINormalBlockVarKnownClusters$clone()`](#method-ZINormalBlockVarKnownClusters-clone)

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
- [`NormalBlockBase$predict()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-predict)
- [`NormalBlockBase$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-update)
- [`NormalBlockBase$warm_start_from()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-warm_start_from)
- [`NormalBlockVarBase$merge()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.html#method-merge)
- [`NormalBlockVarBase$split()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.html#method-split)
- [`NormalBlockVarBase$within_block_correlations()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.html#method-within_block_correlations)

------------------------------------------------------------------------

### `ZINormalBlockVarKnownClusters$new()`

Create a new \[\`ZINormalBlockVarKnownClusters\`\] object.

#### Usage

    ZINormalBlockVarKnownClusters$new(
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

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`ZINormalBlockVarKnownClusters\`\] object

------------------------------------------------------------------------

### `ZINormalBlockVarKnownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINormalBlockVarKnownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3, kappa = rep(0.3, 20))
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = ex$parameters$C, zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE))
model$clustering
#>  [1] 3 3 1 3 2 3 3 2 2 3 2 3 1 2 1 3 2 3 1 1
```
