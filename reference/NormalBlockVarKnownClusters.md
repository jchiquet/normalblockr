# R6 class for a normal-block model with known clustering.

R6 class for a normal-block model with known clustering.

## Super class

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
