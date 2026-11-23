# Root Base Class for Normal-Block Models

R6 abstract class shared by the variance-block
([NormalBlockVarBase](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.md))
and mean-block
([NormalBlockMeanBase](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.md))
model families.

## Public fields

- `data`:

  object of NormalBlockData class, with responses and design matrix

## Active bindings

- `inference_method`:

  inference procedure used (heuristic or integrated with EM)

- `n`:

  number of samples

- `p`:

  number of responses per sample

- `d`:

  number of variables (dimensions in X)

- `d0`:

  number of zero-inflation covariates (dimensions in X0)

- `who_am_I`:

  a description of the model

- `q`:

  number of blocks

- `n_edges`:

  number of edges of the network (non null coefficient of the sparse
  precision matrix Omega)

- `objective`:

  evolution of the objective function during (V)EM algorithm

- `loglik`:

  (or its variational lower bound)

- `deviance`:

  (or its variational lower bound)

- `entropy`:

  entropy of the (variational) posterior distribution of the clustering,
  `-sum(tau * log(tau))`: 0 when the clustering is known. Only the
  discrete partition is completed, as in the ICL of Biernacki, Celeux
  and Govaert (2000); the latent vectors W of the variance-block family
  are integrated out, not completed.

- `penalty_terms`:

  data frame with one row per group of parameters: its name (`group`),
  its number of free parameters (`nu`) and the number of independent
  observations informing it (`N`), so that the group costs `nu * log(N)`
  in the BIC. As in the ICL of Daudin, Picard and Robin (2008) for SBMs
  and Keribin et al. (2015) for latent block models, the cluster
  proportions are charged `log(p)` (the variables are the units being
  clustered), parameters shared by every cell of a block `log(n * p)`,
  and the remaining ones `log(n)`.

- `nb_param`:

  number of free parameters in the model, `sum(penalty_terms$nu)`

- `penalty`:

  BIC penalty, `sum(nu * log(N))` over `penalty_terms`

- `BIC`:

  Bayesian information criterion, `deviance + penalty` (lower is
  better); built on the ELBO when the clustering is unknown and the
  likelihood is not available in closed form.

- `ICL`:

  integrated completed likelihood criterion, `BIC + 2 * entropy` (lower
  is better); equal to the BIC when the clustering is known.

- `EBIC`:

  extended BIC of Foygel and Drton (2010) with gamma = 1/2,
  `BIC + 2 * n_edges * log(dim)`, where `dim` is the dimension of the
  network: q for the variance-block family (network between clusters), p
  for the mean-block family (network between variables).

- `criteria`:

  a vector with loglik, BIC and number of parameters

- `sparsity`:

  (overall sparsity parameter)

- `sparsity_term`:

  (sparsity_term term in log-likelihood due to sparsity)

- `memberships`:

  the p x q matrix of cluster memberships: indicators for a known
  clustering, posterior probabilities (tau) otherwise

- `clustering`:

  the cluster label of each variable (MAP for an unknown clustering),
  named after the columns of Y

- `cluster_sizes`:

  given as a vector of cluster sizes

- `elements_per_cluster`:

  given as the list of elements contained in each cluster

## Methods

### Public methods

- [`NormalBlockBase$new()`](#method-NormalBlockBase-initialize)

- [`NormalBlockBase$update()`](#method-NormalBlockBase-update)

- [`NormalBlockBase$optimize()`](#method-NormalBlockBase-optimize)

- [`NormalBlockBase$best_of_inits()`](#method-NormalBlockBase-best_of_inits)

- [`NormalBlockBase$warm_start_from()`](#method-NormalBlockBase-warm_start_from)

- [`NormalBlockBase$candidates_split()`](#method-NormalBlockBase-candidates_split)

- [`NormalBlockBase$candidates_merge()`](#method-NormalBlockBase-candidates_merge)

- [`NormalBlockBase$heldout_loglik()`](#method-NormalBlockBase-heldout_loglik)

- [`NormalBlockBase$predict()`](#method-NormalBlockBase-predict)

- [`NormalBlockBase$latent_network()`](#method-NormalBlockBase-latent_network)

- [`NormalBlockBase$plot_loglik()`](#method-NormalBlockBase-plot_loglik)

- [`NormalBlockBase$plot_network()`](#method-NormalBlockBase-plot_network)

- [`NormalBlockBase$plot()`](#method-NormalBlockBase-plot)

- [`NormalBlockBase$print()`](#method-NormalBlockBase-print)

- [`NormalBlockBase$clone()`](#method-NormalBlockBase-clone)

------------------------------------------------------------------------

### `NormalBlockBase$new()`

Create a new `NormalBlockBase` object.

#### Usage

    NormalBlockBase$new(data, q, sparsity = 0, control, zero_inflation = FALSE)

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
  `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`) is
  skipped entirely, since it would otherwise never be used downstream.

#### Returns

A new `NormalBlockBase` object

------------------------------------------------------------------------

### `NormalBlockBase$update()`

Update the parameters of a model in place. Each argument applies to the
families or variants named in parentheses; `NA` (the default) leaves the
parameter unchanged.

#### Usage

    NormalBlockBase$update(
      B = NA,
      dm1 = NA,
      C = NA,
      Omega = NA,
      gamma = NA,
      mu = NA,
      alpha = NA,
      M = NA,
      S = NA,
      ll_list = NA,
      warm_started = NA,
      clustering_init = NA
    )

#### Arguments

- `B`:

  regression coefficients (all models)

- `dm1`:

  inverse noise variances of the variables (variance-block)

- `C`:

  cluster memberships: indicators, or posterior probabilities (all
  models)

- `Omega`:

  precision matrix: between clusters (variance-block) or between
  variables (mean-block)

- `gamma`:

  posterior variance of W (variance-block, known clusters)

- `mu`:

  posterior mean of W (variance-block, known clusters)

- `alpha`:

  cluster proportions (unknown clusters)

- `M`:

  variational means of W (variance-block, unknown clusters)

- `S`:

  variational variances of W (variance-block, unknown clusters)

- `ll_list`:

  list of log-lik (elbo) values

- `warm_started`:

  whether `optim_initialize()` should treat the model as already
  initialized (reuse B/Omega/dm1/C/alpha/M/S as they stand) rather than
  recomputing a fresh heuristic initialization. Set by
  `$warm_start_from()` and by `$split()`/`$merge()`.

- `clustering_init`:

  initial clustering

#### Returns

Update the current object in place

------------------------------------------------------------------------

### `NormalBlockBase$optimize()`

calls optimization (EM or heuristic) and updates relevant fields

#### Usage

    NormalBlockBase$optimize(
      control = list(niter = 500, threshold = 1e-04, fixed_point_niter = 5),
      warn = TRUE
    )

#### Arguments

- `control`:

  a list for controlling the optimization process

- `warn`:

  whether to warn when the (V)EM stops at the `niter` cap without
  reaching `threshold` (see `private$warn_if_not_converged()`). Set to
  `FALSE` for deliberately-truncated trial fits (cheap candidate scoring
  in `candidates_split()`/`candidates_merge()`, the sparsity-path
  warm-start probe in
  [NormalBlockCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.md))
  where stopping at the cap is expected and not a sign of trouble.

#### Returns

optimizes the model and updates its parameters

------------------------------------------------------------------------

### `NormalBlockBase$best_of_inits()`

Try several clustering-initialization heuristics and keep the best-ELBO
converged fit (see
[`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)'s
`clustering_init` for the heuristics and the rationale behind the
default). Every candidate is first screened with a short `trial_niter`
run (same idea as `candidates_split()`/`candidates_merge()`), and only
the `max_training` best-screened ones are fully retrained with
`control`.

#### Usage

    NormalBlockBase$best_of_inits(
      inits = private$default_inits,
      trial_niter = 10,
      max_training = 2,
      control = list(niter = 500, threshold = 1e-04, fixed_point_niter = 5)
    )

#### Arguments

- `inits`:

  vector of clustering-heuristic names to try; defaults to the model
  family's own preferred order (`private$default_inits`)

- `trial_niter`:

  number of (V)EM iterations used to cheaply screen every candidate in
  `inits` before fully retraining the best few

- `max_training`:

  how many of the screened candidates (best `loglik` after `trial_niter`
  iterations) get fully retrained with `control`

- `control`:

  [`optimize()`](https://rdrr.io/r/stats/optimize.html) control list
  (`niter`/`threshold`) used for the final full retraining of the
  `max_training` best candidates

#### Returns

a new, already-optimized model. Does not mutate the current object;
reassign the result (`model <- model$best_of_inits()`).

------------------------------------------------------------------------

### `NormalBlockBase$warm_start_from()`

Seed this model's starting parameters from another, already-optimized
model of the same family and q, instead of a fresh heuristic clustering.
Used by the sparsity paths to warm-start each penalty from the previous
one's solution.

#### Usage

    NormalBlockBase$warm_start_from(other)

#### Arguments

- `other`:

  a model of the same class, already optimized

#### Returns

Update the current object in place with `other`'s parameters

------------------------------------------------------------------------

### `NormalBlockBase$candidates_split()`

generate and select a set of candidate models by splitting the clusters
of the current model

#### Usage

    NormalBlockBase$candidates_split(trial_niter = 5, split_rule = NULL)

#### Arguments

- `trial_niter`:

  number of (V)EM iterations used to cheaply score each candidate before
  the caller fully re-optimizes the best few – kept short on purpose.

- `split_rule`:

  passed to [`split()`](https://rdrr.io/r/base/split.html) as its
  `rule`: `"kmeans"`, `"noise"`, or `NULL` for the family's default

------------------------------------------------------------------------

### `NormalBlockBase$candidates_merge()`

generate and select a set of candidate models by merging the clusters of
the current model

#### Usage

    NormalBlockBase$candidates_merge(max_candidates = 30, trial_niter = 2)

#### Arguments

- `max_candidates`:

  merge candidates are, unlike split's, quadratic in q (`choose(q, q-2)`
  pairs): beyond `max_candidates` pairs, only the most promising ones
  are actually built and trial-optimized, ranked by the family's own
  `private$merge_score()`. Set to `Inf` to always try every pair.

- `trial_niter`:

  see `$candidates_split()`

------------------------------------------------------------------------

### `NormalBlockBase$heldout_loglik()`

Log-likelihood of new observations at the fitted parameters and, for an
unknown clustering, at the MAP partition (see
[`normal_block_cv()`](https://jchiquet.github.io/normalblockr/reference/normal_block_cv.md)).
It is exact, computed for variance-block models with q x q matrices only
(Woodbury identity), and given on Y's original scale.

#### Usage

    NormalBlockBase$heldout_loglik(newdata, per_observation = FALSE)

#### Arguments

- `newdata`:

  a
  [NormalBlockData](https://jchiquet.github.io/normalblockr/reference/NormalBlockData.md)
  object with the same variables and covariates (and zero-inflation
  covariates) as the training data.

- `per_observation`:

  whether to return one value per observation rather than their sum.

#### Returns

the held-out log-likelihood, or a vector of `newdata$n` values

------------------------------------------------------------------------

### `NormalBlockBase$predict()`

Predicts the mean of Y for new covariates X, in Y's original units. For
the variance-block family, `X B`: the latent block signals average out,
so that on the training X it differs from
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), which adds
their posterior mean. The mean-block family overrides it.

#### Usage

    NormalBlockBase$predict(new_X)

#### Arguments

- `new_X`:

  new set of covariates.

#### Returns

A n\*p prediction matrix for new observations

------------------------------------------------------------------------

### `NormalBlockBase$latent_network()`

Extract the network encoded by Omega, as a matrix rather than a plot
(see `$plot_network()`): between clusters for the variance-block family,
between variables for the mean-block family.

#### Usage

    NormalBlockBase$latent_network(type = c("partial_cor", "support", "precision"))

#### Arguments

- `type`:

  edge value in the network: "support" (binary edges), "precision"
  (entries of the precision matrix, on the scale of the fit) or
  "partial_cor" (partial correlations)

#### Returns

a sparse square matrix, q x q (variance-block) or p x p (mean-block),
with named rows and columns

------------------------------------------------------------------------

### `NormalBlockBase$plot_loglik()`

plots the evolution of the objective (log-likelihood or ELBO) across the
(V)EM iterations of the last call to
[`optimize()`](https://rdrr.io/r/stats/optimize.html).

#### Usage

    NormalBlockBase$plot_loglik(show_increment = TRUE)

#### Arguments

- `show_increment`:

  whether to add a second panel with the (log10) absolute increment
  between iterations and the convergence `threshold` (distinguishes true
  convergence from a flat-looking objective trace).

#### Returns

a
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graph

------------------------------------------------------------------------

### `NormalBlockBase$plot_network()`

Plot the network encoded by Omega (see `$latent_network()` to extract it
as a matrix).

#### Usage

    NormalBlockBase$plot_network(
      type = c("partial_cor", "support"),
      output = c("igraph", "corrplot"),
      edge.color = c("#F8766D", "#00BFC4"),
      remove.isolated = FALSE,
      node.labels = NULL,
      layout = igraph::layout_in_circle,
      plot = TRUE
    )

#### Arguments

- `type`:

  edge value in the network: "partial_cor" (partial correlations) or
  "support" (binary edges).

- `output`:

  Output type. Either `igraph` (for the network) or `corrplot` (for the
  adjacency matrix)

- `edge.color`:

  Length 2 color vector. Color for positive/negative edges. Default is
  `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.

- `remove.isolated`:

  if `TRUE`, isolated node are remove before plotting. Only relevant for
  igraph output.

- `node.labels`:

  vector of character. The labels of the nodes. By default, the cluster
  names (variance-block) or the variable names (mean-block), as in
  `$latent_network()`.

- `layout`:

  an igraph layout function, applied to the graph without its edge
  weights, or a matrix of node coordinates. Only relevant for igraph
  output.

- `plot`:

  logical. Should the final network be displayed or only sent back to
  the user. Default is `TRUE`.

------------------------------------------------------------------------

### `NormalBlockBase$plot()`

plots the evolution of the objective during model optimization (see
`plot_loglik()`)

#### Usage

    NormalBlockBase$plot()

------------------------------------------------------------------------

### `NormalBlockBase$print()`

User friendly print method

#### Usage

    NormalBlockBase$print(model = paste0("A ", self$who_am_I, ".\n"))

#### Arguments

- `model`:

  First line of the print output

------------------------------------------------------------------------

### `NormalBlockBase$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockBase$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
# An internal abstract base class, never instantiated directly -- see
# normal_block() for how concrete models are created and fitted.
```
