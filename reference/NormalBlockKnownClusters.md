# R6 class for a normal-block model with known clustering.

R6 class for a normal-block model with known clustering.

## Super class

[`NormalBlockBase`](NormalBlockBase.md) -\> `NormalBlockKnownClusters`

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

- [`NormalBlockKnownClusters$new()`](#method-NormalBlockKnownClusters-initialize)

- [`NormalBlockKnownClusters$clone()`](#method-NormalBlockKnownClusters-clone)

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

### `NormalBlockKnownClusters$new()`

Create a new \[\`NormalBlockKnownClusters\`\] object.

#### Usage

    NormalBlockKnownClusters$new(data, C, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `C`:

  clustering matrix C_jk = 1 if species j belongs to cluster k

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`NormalBlockKnownClusters\`\] object

------------------------------------------------------------------------

### `NormalBlockKnownClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockKnownClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
