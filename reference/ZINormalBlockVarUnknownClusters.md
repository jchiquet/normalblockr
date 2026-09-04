# Zero-Inflated Normal-Block Model with Unknown Clustering

R6 class for a zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

## Super classes

[`NormalBlockBase`](NormalBlockBase.md) -\>
[`NormalBlockVarBase`](NormalBlockVarBase.md) -\>
`ZINormalBlockVarUnknownClusters`

## Public fields

- `fixed_tau`:

  whether tau should be fixed at clustering_init during optimization,
  useful for stability selection

## Active bindings

- `nb_param`:

  number of parameters in the model

- `var_par`:

  a list with variational parameters

- `model_par`:

  a list with model parameters: B (covariates), dm1 (species variance),
  Omega (blocks precision matrix), kappa (zero-inflation probabilities)

- `entropy`:

  Entropy of the conditional distribution

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`ZINormalBlockVarUnknownClusters$new()`](#method-ZINormalBlockVarUnknownClusters-initialize)

- [`ZINormalBlockVarUnknownClusters$clone()`](#method-ZINormalBlockVarUnknownClusters-clone)

Inherited methods

- [`NormalBlockBase$best_of_inits()`](NormalBlockBase.html#method-best_of_inits)
- [`NormalBlockBase$candidates_merge()`](NormalBlockBase.html#method-candidates_merge)
- [`NormalBlockBase$candidates_split()`](NormalBlockBase.html#method-candidates_split)
- [`NormalBlockBase$latent_network()`](NormalBlockBase.html#method-latent_network)
- [`NormalBlockBase$optimize()`](NormalBlockBase.html#method-optimize)
- [`NormalBlockBase$plot()`](NormalBlockBase.html#method-plot)
- [`NormalBlockBase$plot_loglik()`](NormalBlockBase.html#method-plot_loglik)
- [`NormalBlockBase$plot_network()`](NormalBlockBase.html#method-plot_network)
- [`NormalBlockBase$predict()`](NormalBlockBase.html#method-predict)
- [`NormalBlockBase$print()`](NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](NormalBlockBase.html#method-update)
- [`NormalBlockVarBase$merge()`](NormalBlockVarBase.html#method-merge)
- [`NormalBlockVarBase$split()`](NormalBlockVarBase.html#method-split)
- [`NormalBlockVarBase$warm_start_from()`](NormalBlockVarBase.html#method-warm_start_from)

------------------------------------------------------------------------

### `ZINormalBlockVarUnknownClusters$new()`

Create a new \[\`ZINormalBlockVarUnknownClusters\`\] object.

#### Usage

    ZINormalBlockVarUnknownClusters$new(
      data,
      q,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `data`:

  object of NormalBlockVarData class, with responses and design matrix

- `q`:

  required number of groups

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters

#### Returns

A new \[\`ZINormalBlockVarUnknownClusters\`\] object

------------------------------------------------------------------------

### `ZINormalBlockVarUnknownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINormalBlockVarUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3, kappa = rep(0.3, 20))
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = 3, zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE))
model$clustering
#>  [1] 1 2 3 3 2 1 3 3 1 2 2 2 1 1 2 3 3 3 3 2
```
