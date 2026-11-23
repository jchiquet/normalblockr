# Base Class for Normal-Block-Mean Models

R6 abstract class for the Normal-Block models where the clustering
structures the mean (mu_i = C B' X_i).

## Super class

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\> `NormalBlockMeanBase`

## Active bindings

- `model_par`:

  a list with the model parameters: B (d x q, one profile per cluster)
  and Omega (p x p precision matrix of the residuals), on the scale of
  the fit (see [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`sigma()`](https://rdrr.io/r/stats/sigma.html)).

- `fitted`:

  Y values predicted by the model, in Y's original units: `X B C'` (with
  the posterior cluster probabilities as C when the clustering is
  unknown), with the observed zeros of zero-inflated models set to 0

- `sparsity_weights`:

  (weights associated to each pair of groups)

## Methods

### Public methods

- [`NormalBlockMeanBase$new()`](#method-NormalBlockMeanBase-initialize)

- [`NormalBlockMeanBase$predict()`](#method-NormalBlockMeanBase-predict)

- [`NormalBlockMeanBase$split()`](#method-NormalBlockMeanBase-split)

- [`NormalBlockMeanBase$merge()`](#method-NormalBlockMeanBase-merge)

- [`NormalBlockMeanBase$clone()`](#method-NormalBlockMeanBase-clone)

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
- [`NormalBlockBase$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-update)
- [`NormalBlockBase$warm_start_from()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.html#method-warm_start_from)

------------------------------------------------------------------------

### `NormalBlockMeanBase$new()`

Create a new `NormalBlockMeanBase` object.

#### Usage

    NormalBlockMeanBase$new(
      data,
      q,
      sparsity = 0,
      control = NB_control(),
      zero_inflation = FALSE
    )

#### Arguments

- `data`:

  object of NormalBlockData class, with responses and design matrix

- `q`:

  number of block/cluster

- `sparsity`:

  sparsity penalty on the network density

- `control`:

  structured list of more specific parameters, to generate with
  NB_control

- `zero_inflation`:

  whether the concrete subclass models zero-inflation; set by the ZI
  subclasses themselves, not meant to be set by the end user.

#### Returns

A new `NormalBlockMeanBase` object

------------------------------------------------------------------------

### `NormalBlockMeanBase$predict()`

Predicts the mean of Y for new covariates X, in Y's original units: mu_i
= C B' X_i, with the posterior cluster probabilities as C when the
clustering is unknown (as in
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)).

#### Usage

    NormalBlockMeanBase$predict(new_X)

#### Arguments

- `new_X`:

  new set of covariates.

#### Returns

A n\*p prediction matrix for new observations

------------------------------------------------------------------------

### `NormalBlockMeanBase$split()`

Create a clone of the current `NormalBlockMeanBase` object after
splitting cluster `index`. Unlike the variance-block family, Omega and
the sparsity weights are p x p here and do not depend on q, so they
carry over unchanged; only C (tau) and B (one column per cluster) are
affected.

#### Usage

    NormalBlockMeanBase$split(index, in_place = FALSE, rule = NULL)

#### Arguments

- `index`:

  index (integer) of the cluster to split

- `in_place`:

  should the split be applied to the object itself, or should a copy be
  sent? Default FALSE (send a copy)

- `rule`:

  how the members of the cluster are divided: `"kmeans"` runs a 2-means
  on their fitted mean trajectories `X B_ols`; `"noise"` cuts them at
  the median of their noise variance `1 / diag(Omega)`. Default `NULL`,
  i.e. `"noise"` for this family.

#### Returns

A new `NormalBlockMeanBase` object

------------------------------------------------------------------------

### `NormalBlockMeanBase$merge()`

Create a clone of the current `NormalBlockMeanBase` object after merging
clusters `indices`

#### Usage

    NormalBlockMeanBase$merge(indices, in_place = FALSE)

#### Arguments

- `indices`:

  indices (couple of integer) of the clusters to merge

- `in_place`:

  should the merge be applied to the object itself, or should a copy be
  sent? Default FALSE (send a copy)

#### Returns

A new `NormalBlockMeanBase` object

------------------------------------------------------------------------

### `NormalBlockMeanBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- use
# NormalBlockMeanKnownClusters / NormalBlockMeanUnknownClusters.
```
