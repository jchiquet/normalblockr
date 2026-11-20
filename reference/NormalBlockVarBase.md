# Base Class for Normal-Block-Var Models

R6 abstract class for the sparse Normal-Block models, where the
clustering structures the latent covariance.

## Super class

[`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
-\> `NormalBlockVarBase`

## Active bindings

- `B_original`:

  regression coefficients (d x p), converted back to Y's original units
  (undoing \`NormalBlockData(scale = TRUE)\`'s column-wise rescaling, if
  any). Use \`model_par\$B\` instead for the coefficients on the
  internal fitting scale.

- `fitted`:

  Y values predicted by the model, in Y's original units: \`X B\`, plus
  the posterior mean of the latent block signals for (V)EM fits, with
  the observed zeros of zero-inflated models set to 0

- `model_par`:

  a list with the model parameters: B (d x p regression coefficients),
  B0 (zero-inflation coefficients, if any), dm1 (inverse noise variances
  of the variables) and Omega (q x q precision matrix of the latent
  block signals), on the scale of the fit. Use
  \`\$B_original\`/\`\$dm1_original\` for B and dm1 on Y's original
  scale.

- `sparsity_weights`:

  (weights associated to each pair of groups)

- `dm1_original`:

  inverse residual variance per variable (1 / Var(Y_j)), converted back
  to Y's original units. Use \`model_par\$dm1\` instead for the internal
  fitting scale. With \`noise_covariance = "spherical"\`,
  \`model_par\$dm1\` is a single value repeated p times; once converted
  back per-variable, the p values returned here generally differ from
  one another whenever Y's columns were rescaled by different factors.

## Methods

### Public methods

- [`NormalBlockVarBase$new()`](#method-NormalBlockVarBase-initialize)

- [`NormalBlockVarBase$within_block_correlations()`](#method-NormalBlockVarBase-within_block_correlations)

- [`NormalBlockVarBase$split()`](#method-NormalBlockVarBase-split)

- [`NormalBlockVarBase$merge()`](#method-NormalBlockVarBase-merge)

- [`NormalBlockVarBase$clone()`](#method-NormalBlockVarBase-clone)

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

------------------------------------------------------------------------

### `NormalBlockVarBase$new()`

Create a new \[\`NormalBlockVarBase\`\] object.

#### Usage

    NormalBlockVarBase$new(
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
  subclasses themselves, not meant to be set by the end user. When
  \`FALSE\`, the (costly) zero-inflation probability fit
  (\`kappa\`/\`B0\`) is skipped entirely, since it would otherwise never
  be used downstream.

#### Returns

A new \[\`NormalBlockVarBase\`\] object

------------------------------------------------------------------------

### `NormalBlockVarBase$within_block_correlations()`

Checks that the blocks gather positively co-varying variables, as the
model requires (two variables of a block have a positive correlation).
Compares, for each block, the empirical correlations of the residuals
\`Y - XB\` with those implied by the model.

#### Usage

    NormalBlockVarBase$within_block_correlations()

#### Returns

a data frame with one row per block: \`block\`, \`size\`, the mean and
minimum of the empirical within-block correlations (\`mean_cor\`,
\`min_cor\`), the proportion of negative ones (\`prop_negative\`), and
the mean correlation implied by the model (\`model_cor\`). \`NA\` for
blocks with fewer than two variables.

------------------------------------------------------------------------

### `NormalBlockVarBase$split()`

Create a clone of the current \[\`NormalBlockVarBase\`\] object after
splitting cluster \`index\`

#### Usage

    NormalBlockVarBase$split(index, in_place = FALSE, rule = NULL)

#### Arguments

- `index`:

  index (integer) of the cluster to split

- `in_place`:

  should the split applied to the object itself, or should a copy be
  sent? default FALSE (send a copy)

- `rule`:

  how the members of the cluster are divided: \`"kmeans"\` runs a
  2-means on their columns of residuals \`Y - XB\`; \`"noise"\` cuts
  them at the median of their noise variance. Default \`NULL\`, i.e.
  \`"kmeans"\` for this family.

#### Returns

A new \[\`NormalBlockVarBase\`\] object

------------------------------------------------------------------------

### `NormalBlockVarBase$merge()`

Create a clone of the current \[\`NormalBlockVarBase\`\] object after
merging clusters \`cl1\` and \`cl2\`

#### Usage

    NormalBlockVarBase$merge(indices, in_place = FALSE)

#### Arguments

- `indices`:

  indices (couple of integer) of the clusters to merge

- `in_place`:

  should the split applied to the object itself, or should a copy be
  sent? default FALSE (send a copy)

#### Returns

A new \[\`NormalBlockVarBase\`\] object

------------------------------------------------------------------------

### `NormalBlockVarBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockVarBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly. See
# normal_block() for how concrete models (NormalBlockVarKnownClusters,
# NormalBlockVarUnknownClusters, and their zero-inflated variants) are
# actually created and fitted.
```
