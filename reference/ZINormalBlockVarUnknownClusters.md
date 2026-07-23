# R6 class for zero-inflated normal-block model with a fixed number of clusters (but unknown clustering).

R6 class for zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

R6 class for zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

## Super class

[`normalblockr::NormalBlockVarBase`](NormalBlockVarBase.md) -\>
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
  Omegaq (blocks precision matrix), kappa (zero-inflation probabilities)

- `entropy`:

  Entropy of the conditional distribution

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`ZINormalBlockVarUnknownClusters$new()`](#method-ZINormalBlockVarUnknownClusters-new)

- [`ZINormalBlockVarUnknownClusters$clone()`](#method-ZINormalBlockVarUnknownClusters-clone)

Inherited methods

- [`normalblockr::NormalBlockVarBase$best_of_inits()`](NormalBlockVarBase.html#method-best_of_inits)
- [`normalblockr::NormalBlockVarBase$candidates_merge()`](NormalBlockVarBase.html#method-candidates_merge)
- [`normalblockr::NormalBlockVarBase$candidates_split()`](NormalBlockVarBase.html#method-candidates_split)
- [`normalblockr::NormalBlockVarBase$latent_network()`](NormalBlockVarBase.html#method-latent_network)
- [`normalblockr::NormalBlockVarBase$merge()`](NormalBlockVarBase.html#method-merge)
- [`normalblockr::NormalBlockVarBase$optimize()`](NormalBlockVarBase.html#method-optimize)
- [`normalblockr::NormalBlockVarBase$plot()`](NormalBlockVarBase.html#method-plot)
- [`normalblockr::NormalBlockVarBase$plot_loglik()`](NormalBlockVarBase.html#method-plot_loglik)
- [`normalblockr::NormalBlockVarBase$plot_network()`](NormalBlockVarBase.html#method-plot_network)
- [`normalblockr::NormalBlockVarBase$predict()`](NormalBlockVarBase.html#method-predict)
- [`normalblockr::NormalBlockVarBase$print()`](NormalBlockVarBase.html#method-print)
- [`normalblockr::NormalBlockVarBase$split()`](NormalBlockVarBase.html#method-split)
- [`normalblockr::NormalBlockVarBase$update()`](NormalBlockVarBase.html#method-update)
- [`normalblockr::NormalBlockVarBase$warm_start_from()`](NormalBlockVarBase.html#method-warm_start_from)

------------------------------------------------------------------------

### Method `new()`

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

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINormalBlockVarUnknownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
