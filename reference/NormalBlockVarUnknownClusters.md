# Normal-Block-Var Model with Unknown Clustering

R6 class for a normal-block model with a fixed number of clusters (but
unknown clustering).

## Super classes

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\>
[`NormalBlockVarBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.md)
-\> `NormalBlockVarUnknownClusters`

## Active bindings

- `model_par`:

  the parameters of \[NormalBlockVarBase\] (B, B0, dm1, Omega), plus
  alpha (cluster proportions)

- `var_par`:

  a list with the matrices of the variational parameters: M (means), S
  (variances), tau (posterior group probabilities)

## Methods

### Public methods

- [`NormalBlockVarUnknownClusters$new()`](#method-NormalBlockVarUnknownClusters-initialize)

- [`NormalBlockVarUnknownClusters$clone()`](#method-NormalBlockVarUnknownClusters-clone)

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

### `NormalBlockVarUnknownClusters$new()`

Create a new \[\`NormalBlockVarUnknownClusters\`\] object.

#### Usage

    NormalBlockVarUnknownClusters$new(
      data,
      q,
      sparsity = 0,
      control = NB_control()
    )

#### Arguments

- `data`:

  contains the matrix of responses (Y) and the design matrix (X).

- `q`:

  required number of groups

- `sparsity`:

  sparsity penalty to add on blocks precision matrix for sparsity

- `control`:

  structured list for specific parameters

#### Returns

A new \[\`NormalBlockVarUnknownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockVarUnknownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
model$clustering
#>  [1] 1 1 2 1 3 2 2 1 2 1 1 3 2 3 1 2 2 2 2 1
```
