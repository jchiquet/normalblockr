# Normal-Block Model with Known Clustering

R6 class for a normal-block model with known clustering.

## Super classes

[`NormalBlockBase`](NormalBlockBase.md) -\>
[`NormalBlockVarBase`](NormalBlockVarBase.md) -\>
`NormalBlockVarKnownClusters`

## Active bindings

- `posterior_par`:

  a list with the parameters of posterior distribution W \| Y

- `entropy`:

  Entropy of the conditional distribution

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarKnownClusters$new()`](#method-NormalBlockVarKnownClusters-initialize)

- [`NormalBlockVarKnownClusters$clone()`](#method-NormalBlockVarKnownClusters-clone)

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

### `NormalBlockVarKnownClusters$new()`

Create a new \[\`NormalBlockVarKnownClusters\`\] object.

#### Usage

    NormalBlockVarKnownClusters$new(data, C, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NormalBlockVarData class, with responses and design matrix

- `C`:

  clustering matrix C_jk = 1 if species j belongs to cluster k

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockVarKnownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockVarKnownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarKnownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = ex$parameters$C, control = NB_control(verbose = FALSE))
model$clustering
#>  [1] 2 3 2 1 1 3 2 3 3 3 3 1 3 2 1 1 3 1 3 2
```
