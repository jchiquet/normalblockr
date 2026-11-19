# Root Base Class for Normal-Block Models

R6 abstract class shared by the variance-block (\[NormalBlockVarBase\])
and mean-block (\[NormalBlockMeanBase\]) model families.

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
  \`-sum(tau \* log(tau))\`: 0 when the clustering is known. Only the
  discrete partition is completed, as in the ICL of Biernacki, Celeux
  and Govaert (2000); the latent vectors W of the variance-block family
  are integrated out, not completed.

- `penalty_terms`:

  data frame with one row per group of parameters: its name (\`group\`),
  its number of free parameters (\`nu\`) and the number of independent
  observations informing it (\`N\`), so that the group costs \`nu \*
  log(N)\` in the BIC. As in the ICL of Daudin, Picard and Robin (2008)
  for SBMs and Keribin et al. (2015) for latent block models, the
  cluster proportions are charged \`log(p)\` (the variables are the
  units being clustered), parameters shared by every cell of a block
  \`log(n \* p)\`, and the remaining ones \`log(n)\`.

- `nb_param`:

  number of free parameters in the model, \`sum(penalty_terms\$nu)\`

- `penalty`:

  BIC penalty, \`sum(nu \* log(N))\` over \`penalty_terms\`

- `BIC`:

  Bayesian information criterion, \`deviance + penalty\` (lower is
  better); built on the ELBO when the clustering is unknown and the
  likelihood is not available in closed form.

- `ICL`:

  integrated completed likelihood criterion, \`BIC + 2 \* entropy\`
  (lower is better); equal to the BIC when the clustering is known.

- `EBIC`:

  extended BIC of Foygel and Drton (2010) with gamma = 1/2, \`BIC + 2 \*
  n_edges \* log(dim)\`, where \`dim\` is the dimension of the network:
  q for the variance-block family (network between clusters), p for the
  mean-block family (network between variables).

- `criteria`:

  a vector with loglik, BIC and number of parameters

- `sparsity`:

  (overall sparsity parameter)

- `sparsity_term`:

  (sparsity_term term in log-likelihood due to sparsity)

- `memberships`:

  cluster memberships

- `clustering`:

  given as the list of elements contained in each cluster

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

Create a new \[\`NormalBlockBase\`\] object.

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
  \`FALSE\`, the (costly) zero-inflation probability fit
  (\`kappa\`/\`B0\`) is skipped entirely, since it would otherwise never
  be used downstream.

#### Returns

A new \[\`NormalBlockBase\`\] object

------------------------------------------------------------------------

### `NormalBlockBase$update()`

Update a \[\`NormalBlockBase\`\] object

All possible parameters of the child classes

#### Usage

    NormalBlockBase$update(
      B = NA,
      dm1 = NA,
      C = NA,
      Omega = NA,
      gamma = NA,
      mu = NA,
      kappa = NA,
      alpha = NA,
      M = NA,
      S = NA,
      Psi = NA,
      Phi = NA,
      Lambda = NA,
      ll_list = NA,
      warm_started = NA,
      clustering_init = NA
    )

#### Arguments

- `B`:

  regression matrix \[all\]

- `dm1`:

  diagonal vector of inverse variance matrix (variables level) \[NBVar\]

- `C`:

  the matrix of groups memberships (posterior probabilities) \[all\]

- `Omega`:

  inverse variance matrix (cluster-level for Normal Block models,
  variable-level for Normal Mean Block models) \[all\]

- `gamma`:

  variance of posterior distribution of W \[NBVar - known\]

- `mu`:

  mean for posterior distribution of W \[NBVar - known\]

- `kappa`:

  vector of zero-inflation probabilities \[ZINBVar, ZINBMean\]

- `alpha`:

  vector of groups probabilities \[NBVar\]

- `M`:

  variational mean for posterior distribution of W \[NBVar - unknown\]

- `S`:

  variational diagonal of variances for posterior distribution of W
  \[NBVar - unknown\]

- `Psi`:

  variational expectation of C'Omega C, intermediary term in
  calculations \[NBMean - unknown\]

- `Phi`:

  variational correction term used in the Psi/ELBO computations
  \[NBMean - unknown\]

- `Lambda`:

  variational correction term used in the Sigma-hat update \[NBMean -
  unknown\]

- `ll_list`:

  list of log-lik (elbo) values

- `warm_started`:

  whether \`optim_initialize()\` should treat the model as already
  initialized (reuse B/Omega/dm1/C/alpha/M/S as they stand) rather than
  recomputing a fresh heuristic initialization. Set by
  \[warm_start_from()\] and by \[split()\]/\[merge()\].

- `clustering_init`:

  initial clustering

#### Returns

Update the current \[\`normal\`\] object

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

  whether to warn when the (V)EM stops at the \`niter\` cap without
  reaching \`threshold\` (see \`private\$warn_if_not_converged()\`). Set
  to \`FALSE\` for deliberately-truncated trial fits (cheap candidate
  scoring in \`candidates_split()\`/\`candidates_merge()\`, the
  sparsity-path warm-start probe in \[NormalBlockCollectionSparsity\])
  where stopping at the cap is expected and not a sign of trouble.

#### Returns

optimizes the model and updates its parameters

------------------------------------------------------------------------

### `NormalBlockBase$best_of_inits()`

Try several clustering-initialization heuristics and keep the best-ELBO
converged fit (see \[NB_control()\]'s \`clustering_init\` for the
heuristics and the rationale behind the default). Every candidate is
first screened with a short \`trial_niter\` run (same idea as
\`candidates_split()\`/\`candidates_merge()\`), and only the
\`max_training\` best-screened ones are fully retrained with
\`control\`.

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
  family's own preferred order (\`private\$default_inits\`)

- `trial_niter`:

  number of (V)EM iterations used to cheaply screen every candidate in
  \`inits\` before fully retraining the best few

- `max_training`:

  how many of the screened candidates (best \`loglik\` after
  \`trial_niter\` iterations) get fully retrained with \`control\`

- `control`:

  \`optimize()\` control list (\`niter\`/\`threshold\`) used for the
  final full retraining of the \`max_training\` best candidates

#### Returns

a new, already-optimized model. Does not mutate the current object;
reassign the result (\`model \<- model\$best_of_inits()\`).

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

  passed to \`split()\` as its \`rule\`: \`"kmeans"\`, \`"noise"\`, or
  \`NULL\` for the family's default

------------------------------------------------------------------------

### `NormalBlockBase$candidates_merge()`

generate and select a set of candidate models by merging the clusters of
the current model

#### Usage

    NormalBlockBase$candidates_merge(max_candidates = 30, trial_niter = 2)

#### Arguments

- `max_candidates`:

  merge candidates are, unlike split's, quadratic in q (\`choose(q,
  q-2)\` pairs): beyond \`max_candidates\` pairs, only the most
  promising ones are actually built and trial-optimized, ranked by the
  family's own \`private\$merge_score()\`. Set to \`Inf\` to always try
  every pair.

- `trial_niter`:

  see \[candidates_split()\]

------------------------------------------------------------------------

### `NormalBlockBase$heldout_loglik()`

Log-likelihood of new observations at the fitted parameters and, for an
unknown clustering, at the MAP partition (see \[normal_block_cv()\]). It
is exact, computed for variance-block models with q x q matrices only
(Woodbury identity), and given on Y's original scale.

#### Usage

    NormalBlockBase$heldout_loglik(newdata, per_observation = FALSE)

#### Arguments

- `newdata`:

  a \[NormalBlockData\] object with the same variables and covariates
  (and zero-inflation covariates) as the training data.

- `per_observation`:

  whether to return one value per observation rather than their sum.

#### Returns

the held-out log-likelihood, or a vector of \`newdata\$n\` values

------------------------------------------------------------------------

### `NormalBlockBase$predict()`

Predicts observations Y for new covariates X, in Y's original units
(like \`\$fitted\`, so that predicting on the training X reproduces it).

#### Usage

    NormalBlockBase$predict(new_X)

#### Arguments

- `new_X`:

  new set of covariates.

#### Returns

A n\*p prediction matrix for new observations

------------------------------------------------------------------------

### `NormalBlockBase$latent_network()`

Extract interaction network in the latent space, as a matrix rather than
a plot – see \`\$plot_network()\` to plot it instead.

#### Usage

    NormalBlockBase$latent_network(type = c("partial_cor", "support", "precision"))

#### Arguments

- `type`:

  edge value in the network. Can be "support" (binary edges),
  "precision" (coefficient of the precision matrix) or "partial_cor"
  (partial correlation between species)

#### Returns

a square matrix of size \`self\$q\`

------------------------------------------------------------------------

### `NormalBlockBase$plot_loglik()`

plots the evolution of the objective (log-likelihood or ELBO) across the
(V)EM iterations of the last call to \`optimize()\`.

#### Usage

    NormalBlockBase$plot_loglik(show_increment = TRUE)

#### Arguments

- `show_increment`:

  whether to add a second panel with the (log10) absolute increment
  between iterations and the convergence \`threshold\` (distinguishes
  true convergence from a flat-looking objective trace).

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### `NormalBlockBase$plot_network()`

plot the latent network. To extract the network as a matrix instead of
plotting it, use \`\$latent_network()\`.

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

  edge value in the network. Either "precision" (coefficient of the
  precision matrix) or "partial_cor" (partial correlation between
  species).

- `output`:

  Output type. Either \`igraph\` (for the network) or \`corrplot\` (for
  the adjacency matrix)

- `edge.color`:

  Length 2 color vector. Color for positive/negative edges. Default is
  \`c("#F8766D", "#00BFC4")\`. Only relevant for igraph output.

- `remove.isolated`:

  if \`TRUE\`, isolated node are remove before plotting. Only relevant
  for igraph output.

- `node.labels`:

  vector of character. The labels of the nodes. The default will use the
  column names ot the response matrix.

- `layout`:

  an igraph layout function, applied to the graph without its edge
  weights, or a matrix of node coordinates. Only relevant for igraph
  output.

- `plot`:

  logical. Should the final network be displayed or only sent back to
  the user. Default is \`TRUE\`.

------------------------------------------------------------------------

### `NormalBlockBase$plot()`

plots the evolution of the objective during model optimization (see
\`plot_loglik()\`)

#### Usage

    NormalBlockBase$plot()

------------------------------------------------------------------------

### `NormalBlockBase$print()`

User friendly print method

#### Usage

    NormalBlockBase$print(model = paste("A", self$who_am_I, ".\n"))

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
