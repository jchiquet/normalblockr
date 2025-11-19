# R6 class for zero-inflated normal-block model with a fixed number of clusters (but unknown clustering).

R6 class for zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

R6 class for zero-inflated normal-block model with a fixed number of
clusters (but unknown clustering).

## Super class

[`normalblockr::NB`](NB.md) -\> `ZINB_fixed_q`

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

  Y values predicted by the model Y values predicted by the model

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`ZINB_fixed_q$new()`](#method-ZINB_fixed_q-new)

- [`ZINB_fixed_q$clone()`](#method-ZINB_fixed_q-clone)

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

Create a new \[\`ZINB_fixed_q\`\] object.

#### Usage

    ZINB_fixed_q$new(data, q, sparsity = 0, control = NB_control())

#### Arguments

- `data`:

  object of NB_data class, with responses and design matrix

- `q`:

  required number of groups

- `sparsity`:

  to apply on variance matrix when calling GLASSO

- `control`:

  structured list of more specific parameters

#### Returns

A new \[\`ZINB_fixed_q\`\] object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    ZINB_fixed_q$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
