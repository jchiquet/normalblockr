# R6 class for a normal-block model with fixed number of clusters (but unknown clustering).

R6 class for a normal-block model with fixed number of clusters (but
unknown clustering).

R6 class for a normal-block model with fixed number of clusters (but
unknown clustering).

## Super class

[`normalblockr::NB`](NB.md) -\> `NB_fixed_q`

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

  Y values predicted by the model

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NB_fixed_q$new()`](#method-NB_fixed_q-new)

- [`NB_fixed_q$clone()`](#method-NB_fixed_q-clone)

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

Create a new \[\`NB_fixed_q\`\] object.

#### Usage

    NB_fixed_q$new(data, q, sparsity = 0, control = NB_control())

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

A new \[\`NB_fixed_q\`\] object

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NB_fixed_q$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
