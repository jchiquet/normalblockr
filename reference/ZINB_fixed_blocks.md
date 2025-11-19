# R6 class for a Zero-Inflated normal-block model with a known clustering.

R6 class for a Zero-Inflated normal-block model with a known clustering.

R6 class for a Zero-Inflated normal-block model with a known clustering.

## Super class

[`normalblockr::NB`](NB.md) -\> `ZINB_fixed_blocks`

## Active bindings

- `posterior_par`:

  a list with the parameters of posterior distribution W \| Y

- `entropy`:

  Entropy of the conditional distribution

- `nb_param`:

  number of parameters in the model

- `model_par`:

  a list with model parameters: B (covariates), dm1 (species variance),
  Omegaq (groups precision matrix), kappa (zero-inflation probabilities)

- `fitted`:

  Y values predicted by the model

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`ZINB_fixed_blocks$new()`](#method-ZINB_fixed_blocks-new)

- [`ZINB_fixed_blocks$clone()`](#method-ZINB_fixed_blocks-clone)

Inherited methods

- [`normalblockr::NB$candidates_merge()`](NB.html#method-candidates_merge)
- [`normalblockr::NB$candidates_split()`](NB.html#method-candidates_split)
- [`normalblockr::NB$latent_network()`](NB.html#method-latent_network)
- [`normalblockr::NB$merge()`](NB.html#method-merge)
- [`normalblockr::NB$optimize()`](NB.html#method-optimize)
- [`normalblockr::NB$plot()`](NB.html#method-plot)
- [`normalblockr::NB$plot_loglik()`](NB.html#method-plot_loglik)
- [`normalblockr::NB$plot_network()`](NB.html#method-plot_network)
- [`normalblockr::NB$predict()`](NB.html#method-predict)
- [`normalblockr::NB$print()`](NB.html#method-print)
- [`normalblockr::NB$split()`](NB.html#method-split)
- [`normalblockr::NB$update()`](NB.html#method-update)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`ZINB_fixed_blocks_fixed_sparsity\`\] object.

#### Usage

    ZINB_fixed_blocks$new(data, C, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NB_data class, with responses and design matrix

- `C`:

  clustering matrix C_jk = 1 if species j belongs to cluster k

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

#### Returns

A new \[\`ZINB_fixed_blocks_fixed_sparsity\`\] object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINB_fixed_blocks$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
