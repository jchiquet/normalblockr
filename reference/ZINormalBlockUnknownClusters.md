# R6 class for zero-inflated normal-block model with a fixed number of clusters (but unknown clustering).

R6 class for zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

## Super class

[`NormalBlockBase`](NormalBlockBase.md) -\>
`ZINormalBlockUnknownClusters`

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
  Omegaq (blocks precision matrix), kappa (zero-inflation probabilities)

- `entropy`:

  Entropy of the conditional distribution

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`ZINormalBlockUnknownClusters$new()`](#method-ZINormalBlockUnknownClusters-initialize)

- [`ZINormalBlockUnknownClusters$clone()`](#method-ZINormalBlockUnknownClusters-clone)

Inherited methods

- [`NormalBlockBase$candidates_merge()`](NormalBlockBase.html#method-candidates_merge)
- [`NormalBlockBase$candidates_split()`](NormalBlockBase.html#method-candidates_split)
- [`NormalBlockBase$latent_network()`](NormalBlockBase.html#method-latent_network)
- [`NormalBlockBase$merge()`](NormalBlockBase.html#method-merge)
- [`NormalBlockBase$optimize()`](NormalBlockBase.html#method-optimize)
- [`NormalBlockBase$plot()`](NormalBlockBase.html#method-plot)
- [`NormalBlockBase$plot_loglik()`](NormalBlockBase.html#method-plot_loglik)
- [`NormalBlockBase$plot_network()`](NormalBlockBase.html#method-plot_network)
- [`NormalBlockBase$predict()`](NormalBlockBase.html#method-predict)
- [`NormalBlockBase$print()`](NormalBlockBase.html#method-print)
- [`NormalBlockBase$split()`](NormalBlockBase.html#method-split)
- [`NormalBlockBase$update()`](NormalBlockBase.html#method-update)
- [`NormalBlockBase$warm_start_from()`](NormalBlockBase.html#method-warm_start_from)

------------------------------------------------------------------------

### `ZINormalBlockUnknownClusters$new()`

Create a new \[\`ZINormalBlockUnknownClusters\`\] object.

#### Usage

    ZINormalBlockUnknownClusters$new(data, q, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `q`:

  required number of groups

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters

#### Returns

A new \[\`ZINormalBlockUnknownClusters\`\] object

------------------------------------------------------------------------

### `ZINormalBlockUnknownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINormalBlockUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
