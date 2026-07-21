# R6 abstract class for a generic sparse Normal Block model

R6 abstract class for a generic sparse Normal Block model

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

  number of zi variables (dimensions in X0)

- `q`:

  number of blocks

- `n_edges`:

  number of edges of the network (non null coefficient of the sparse
  precision matrix Omegaq)

- `model_par`:

  a list with the matrices of the model parameters: B (covariates), dm1
  (species variance), Omegaq (groups precision matrix)). On the internal
  fitting scale (\`self\$data\$Y\`, possibly column-rescaled by
  \`NormalBlockData(scale = TRUE)\`) – use
  \`\$B_original\`/\`\$dm1_original\` for the same quantities converted
  back to Y's original units.

- `B_original`:

  regression coefficients (d x p), converted back to Y's original units
  (undoing \`NormalBlockData(scale = TRUE)\`'s column-wise rescaling, if
  any). Use \`model_par\$B\` instead for the coefficients on the
  internal fitting scale.

- `dm1_original`:

  inverse residual variance per variable (1 / Var(Y_j)), converted back
  to Y's original units. Use \`model_par\$dm1\` instead for the internal
  fitting scale. With \`noise_covariance = "spherical"\`,
  \`model_par\$dm1\` is a single value repeated p times (one shared
  variance on the fitting scale); once converted back per-variable, the
  p values returned here generally differ from one another whenever Y's
  columns were rescaled by different factors – correctly so, since a
  single shared \*scaled\* variance does not correspond to a single
  shared variance in the original, heterogeneous-scale units.

- `nb_param`:

  number of parameters in the model

- `objective`:

  evolution of the objective function during (V)EM algorithm

- `loglik`:

  (or its variational lower bound)

- `deviance`:

  (or its variational lower bound)

- `BIC`:

  (or its variational lower bound)

- `entropy`:

  Entropy of the conditional distribution when applicable

- `ICL`:

  variational lower bound of the ICL

- `EBIC`:

  variational lower bound of the EBIC

- `criteria`:

  a vector with loglik, BIC and number of parameters

- `sparsity`:

  (overall sparsity parameter)

- `sparsity_weights`:

  (weights associated to each pair of groups)

- `sparsity_term`:

  (sparsity_term term in log-likelihood due to sparsity)

- `get_res_covariance`:

  whether the residual covariance is diagonal or spherical

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

- [`NormalBlockBase$best_of_inits()`](#method-NormalBlockBase-best_of_inits)

- [`NormalBlockBase$optimize()`](#method-NormalBlockBase-optimize)

- [`NormalBlockBase$warm_start_from()`](#method-NormalBlockBase-warm_start_from)

- [`NormalBlockBase$split()`](#method-NormalBlockBase-split)

- [`NormalBlockBase$candidates_split()`](#method-NormalBlockBase-candidates_split)

- [`NormalBlockBase$candidates_merge()`](#method-NormalBlockBase-candidates_merge)

- [`NormalBlockBase$merge()`](#method-NormalBlockBase-merge)

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

    NormalBlockBase$new(
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
      Omegaq = NA,
      gamma = NA,
      mu = NA,
      kappa = NA,
      alpha = NA,
      M = NA,
      S = NA,
      ll_list = NA,
      warm_started = NA,
      clustering_init = NA
    )

#### Arguments

- `B`:

  regression matrix

- `dm1`:

  diagonal vector of inverse variance matrix (variables level)

- `C`:

  the matrix of groups memberships (posterior probabilities)

- `Omegaq`:

  groups inverse variance matrix

- `gamma`:

  variance of posterior distribution of W

- `mu`:

  mean for posterior distribution of W

- `kappa`:

  vector of zero-inflation probabilities

- `alpha`:

  vector of groups probabilities

- `M`:

  variational mean for posterior distribution of W

- `S`:

  variational diagonal of variances for posterior distribution of W

- `ll_list`:

  list of log-lik (elbo) values

- `warm_started`:

  whether \`EM_initialize()\` should treat the model as already
  initialized (reuse B/Omegaq/dm1/C/alpha/M/S as they stand) rather than
  recomputing a fresh heuristic initialization – set by
  \[warm_start_from()\] and by \[split()\]/\[merge()\].

- `clustering_init`:

  name of a clustering heuristic to switch to, re-derived at the next
  \`optimize()\` call instead of reusing the current state (see
  \`NB_control(clustering_init = )\`). Used by \[best_of_inits()\].

#### Returns

Update the current \[\`normal\`\] object

------------------------------------------------------------------------

### `NormalBlockBase$best_of_inits()`

Try several clustering-initialization heuristics and keep the best-ELBO
converged fit (see \`NB_control(clustering_init = )\` and
\`inst/methods_initialization_and_refine.md\` for the rationale). Every
candidate is first screened with a short \`trial_niter\` run (same idea
as \`candidates_split()\`/\`candidates_merge()\`), and only the
\`max_training\` best-screened ones are fully retrained with
\`control\`.

#### Usage

    NormalBlockBase$best_of_inits(
      inits = c("ward2", "kmeans", "spectral"),
      trial_niter = 10,
      max_training = 2,
      control = list(niter = 500, threshold = 1e-04)
    )

#### Arguments

- `inits`:

  vector of clustering-heuristic names to try

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

a new, already-optimized \[\`NormalBlockBase\`\] object. Does not mutate
the current object; reassign the result (\`model \<-
model\$best_of_inits()\`).

------------------------------------------------------------------------

### `NormalBlockBase$optimize()`

calls optimization (EM or heuristic) and updates relevant fields

#### Usage

    NormalBlockBase$optimize(
      control = list(niter = 500, threshold = 1e-04),
      warn = TRUE
    )

#### Arguments

- `control`:

  a list for controlling the optimization proces

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

### `NormalBlockBase$warm_start_from()`

Seed this model's starting parameters from another, already-optimized
model with the same q, instead of the heuristic clustering-derived
values set at construction time. Used by
\[NormalBlockCollectionSparsity\] to warm-start each penalty in a
sparsity path from the previous (adjacent) one's converged solution –
adjacent penalties along a sorted path usually have similar optima, so
this typically needs far fewer EM iterations than starting cold each
time (the same rationale as warm-starting in glmnet/glassoFast's own
regularization paths). \`B0\`/\`kappa\` (zero-inflation) are
deliberately left untouched: they depend only on the data, not on
sparsity/blocks, so they are already set correctly and independently on
every model.

#### Usage

    NormalBlockBase$warm_start_from(other)

#### Arguments

- `other`:

  a \[NormalBlockBase\] object, already optimized

#### Returns

Update the current object in place with \`other\`'s parameters

------------------------------------------------------------------------

### `NormalBlockBase$split()`

Create a clone of the current \[\`NormalBlockBase\`\] object after
splitting cluster \`cl\` We split the cluster according to the species
variances

#### Usage

    NormalBlockBase$split(index, in_place = FALSE)

#### Arguments

- `index`:

  index (integer) of the cluster to split

- `in_place`:

  should the split applied to the object itself, or should a copy be
  sent? default FALSE (send a copy)

#### Returns

A new \[\`NormalBlockBase\`\] object

------------------------------------------------------------------------

### `NormalBlockBase$candidates_split()`

generate and select a set of candidate models by splitting the clusters
of the current model

#### Usage

    NormalBlockBase$candidates_split(trial_niter = 5)

#### Arguments

- `trial_niter`:

  number of EM iterations used to cheaply score each candidate before
  \[SelectionNClusters\] fully re-optimizes the best few
  (\`train_best_candidates()\`'s \`max_training\`) – kept short on
  purpose.

------------------------------------------------------------------------

### `NormalBlockBase$candidates_merge()`

generate and select a set of candidate models by merging the clusters of
the current model

#### Usage

    NormalBlockBase$candidates_merge(max_candidates = 30, trial_niter = 2)

#### Arguments

- `max_candidates`:

  merge candidates are, unlike split's, quadratic in q (\`choose(q,
  q-2)\` pairs) – beyond \`max_candidates\` pairs, only the most
  promising ones (largest \`\|Omegaq\[i, j\]\|\`, i.e. the most strongly
  related cluster pairs in the current fit) are actually built and
  trial-optimized, since merging two nearly independent blocks is rarely
  competitive anyway. Set to \`Inf\` to always try every pair.

- `trial_niter`:

  see \[candidates_split()\]

------------------------------------------------------------------------

### `NormalBlockBase$merge()`

Create a clone of the current \[\`NormalBlockBase\`\] object after
merging clusters \`cl1\` and \`cl2\`

#### Usage

    NormalBlockBase$merge(indices, in_place = FALSE)

#### Arguments

- `indices`:

  indices (couple of integer) of the clusters to merge

- `in_place`:

  should the split applied to the object itself, or should a copy be
  sent? default FALSE (send a copy)

#### Returns

A new \[\`NormalBlockBase\`\] object

------------------------------------------------------------------------

### `NormalBlockBase$predict()`

Predicts observations Y for new covariates X.

#### Usage

    NormalBlockBase$predict(new_X)

#### Arguments

- `new_X`:

  new set of covariates.

#### Returns

A n\*p prediction matrix for new observations

------------------------------------------------------------------------

### `NormalBlockBase$latent_network()`

Extract interaction network in the latent space

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

  whether to add, below the objective trace, a second panel with the
  (log10) absolute increment between consecutive iterations and the
  convergence \`threshold\` used to stop optimize() (dashed line). That
  second panel is what actually tells convergence apart from merely
  running out of iterations: the objective trace alone tends to look
  flat well before the increment has actually crossed the threshold,
  especially as the number of blocks grows (see inst/CSDA_analyses).

#### Returns

a \[\`ggplot2::ggplot\`\] graph

------------------------------------------------------------------------

### `NormalBlockBase$plot_network()`

plot the latent network.

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

  an optional igraph layout. Only relevant for igraph output.

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
