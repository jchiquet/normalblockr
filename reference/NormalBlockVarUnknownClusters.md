# R6 class for a normal-block model with fixed number of clusters (but unknown clustering).

R6 class for a normal-block model with fixed number of clusters (but
unknown clustering).

## Super class

[`NormalBlockVarBase`](NormalBlockVarBase.md) -\>
`NormalBlockVarUnknownClusters`

## Public fields

- `fixed_tau`:

  whether tau should be fixed at clustering_init during optimization,
  useful for stability selection

## Active bindings

- `model_par`:

  a list with the matrices of the model parameters: B (covariates), dm1
  (species variance), Omegaq (groups precision matrix))

- `nb_param`:

  number of parameters in the model

- `var_par`:

  a list with the matrices of the variational parameters: M (means), S
  (variances), tau (posterior group probabilities)

- `entropy`:

  Entropy of the conditional distribution

- `fitted`:

  Y values predicted by the model, in Y's original units

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockVarUnknownClusters$new()`](#method-NormalBlockVarUnknownClusters-initialize)

- [`NormalBlockVarUnknownClusters$clone()`](#method-NormalBlockVarUnknownClusters-clone)

Inherited methods

- [`NormalBlockVarBase$best_of_inits()`](NormalBlockVarBase.html#method-best_of_inits)
- [`NormalBlockVarBase$candidates_merge()`](NormalBlockVarBase.html#method-candidates_merge)
- [`NormalBlockVarBase$candidates_split()`](NormalBlockVarBase.html#method-candidates_split)
- [`NormalBlockVarBase$latent_network()`](NormalBlockVarBase.html#method-latent_network)
- [`NormalBlockVarBase$merge()`](NormalBlockVarBase.html#method-merge)
- [`NormalBlockVarBase$optimize()`](NormalBlockVarBase.html#method-optimize)
- [`NormalBlockVarBase$plot()`](NormalBlockVarBase.html#method-plot)
- [`NormalBlockVarBase$plot_loglik()`](NormalBlockVarBase.html#method-plot_loglik)
- [`NormalBlockVarBase$plot_network()`](NormalBlockVarBase.html#method-plot_network)
- [`NormalBlockVarBase$predict()`](NormalBlockVarBase.html#method-predict)
- [`NormalBlockVarBase$print()`](NormalBlockVarBase.html#method-print)
- [`NormalBlockVarBase$split()`](NormalBlockVarBase.html#method-split)
- [`NormalBlockVarBase$update()`](NormalBlockVarBase.html#method-update)
- [`NormalBlockVarBase$warm_start_from()`](NormalBlockVarBase.html#method-warm_start_from)

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
