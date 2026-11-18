# Base Class for Variance-Block Models

R6 abstract class for the sparse Normal-Block models, where the
clustering structures the latent covariance.

## Super class

[`NormalBlockBase`](NormalBlockBase.md) -\> `NormalBlockVarBase`

## Active bindings

- `B_original`:

  regression coefficients (d x p), converted back to Y's original units
  (undoing \`NormalBlockData(scale = TRUE)\`'s column-wise rescaling, if
  any). Use \`model_par\$B\` instead for the coefficients on the
  internal fitting scale.

- `d0`:

  number of zi variables (dimensions in X0)

- `model_par`:

  a list with the matrices of the model parameters: B (covariates), dm1
  (species variance), Omega (groups precision matrix)). On the internal
  fitting scale (\`self\$data\$Y\`, possibly column-rescaled by
  \`NormalBlockData(scale = TRUE)\`) – use
  \`\$B_original\`/\`\$dm1_original\` for the same quantities converted
  back to Y's original units.

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

- [`NormalBlockVarBase$warm_start_from()`](#method-NormalBlockVarBase-warm_start_from)

- [`NormalBlockVarBase$within_block_correlations()`](#method-NormalBlockVarBase-within_block_correlations)

- [`NormalBlockVarBase$split()`](#method-NormalBlockVarBase-split)

- [`NormalBlockVarBase$merge()`](#method-NormalBlockVarBase-merge)

- [`NormalBlockVarBase$clone()`](#method-NormalBlockVarBase-clone)

Inherited methods

- [`NormalBlockBase$best_of_inits()`](NormalBlockBase.html#method-best_of_inits)
- [`NormalBlockBase$candidates_merge()`](NormalBlockBase.html#method-candidates_merge)
- [`NormalBlockBase$candidates_split()`](NormalBlockBase.html#method-candidates_split)
- [`NormalBlockBase$heldout_loglik()`](NormalBlockBase.html#method-heldout_loglik)
- [`NormalBlockBase$latent_network()`](NormalBlockBase.html#method-latent_network)
- [`NormalBlockBase$optimize()`](NormalBlockBase.html#method-optimize)
- [`NormalBlockBase$plot()`](NormalBlockBase.html#method-plot)
- [`NormalBlockBase$plot_loglik()`](NormalBlockBase.html#method-plot_loglik)
- [`NormalBlockBase$plot_network()`](NormalBlockBase.html#method-plot_network)
- [`NormalBlockBase$predict()`](NormalBlockBase.html#method-predict)
- [`NormalBlockBase$print()`](NormalBlockBase.html#method-print)
- [`NormalBlockBase$update()`](NormalBlockBase.html#method-update)

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

### `NormalBlockVarBase$warm_start_from()`

Seed this model's starting parameters from another, already-optimized
model with the same q, instead of a fresh heuristic clustering. Used by
\[NormalBlockVarCollectionSparsity\] to warm-start each penalty in a
sparsity path from the previous one's solution.

#### Usage

    NormalBlockVarBase$warm_start_from(other)

#### Arguments

- `other`:

  a \[NormalBlockVarBase\] object, already optimized

#### Returns

Update the current object in place with \`other\`'s parameters

------------------------------------------------------------------------

### `NormalBlockVarBase$within_block_correlations()`

Diagnose whether the blocks are made of positively co-varying variables,
as the model requires. Within block k, two variables j and l are noisy
copies of the same latent signal, so their covariance is \`Sigma_kk \>
0\` and their correlation \`Sigma_kk / sqrt((Sigma_kk + d_j) (Sigma_kk +
d_l))\` lies in (0, 1): the model cannot put two anticorrelated
variables in the same block. This compares, for each block, the
empirical correlations of the variables' residuals \`Y - XB\` (non-zero
cells only for zero-inflated models) with the correlations the fitted
model implies. Mostly useful for a known clustering (e.g. from an
annotation): a block with many negative empirical correlations is
misspecified, and should rather be split into blocks linked by a
negative \`Sigma_kl\`.

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
  2-means on their columns of residuals \`Y - XB\`, the geometry of the
  model (members of a block are noisy copies of a common latent signal,
  so they lie close to each other); \`"noise"\` cuts them at the median
  of their noise variance (the rule used before 0.4.0). Default
  \`NULL\`, resolved to \`"kmeans"\` for this family: on \`brca_rppa\`
  it lets \`refine()\` lower the deviance summed over q = 1..30 by
  22,166 against 18,588 for \`"noise"\`.

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
