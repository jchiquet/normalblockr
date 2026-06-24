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
  (species variance), Omegaq (groups precision matrix))

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

- [`NormalBlockBase$optimize()`](#method-NormalBlockBase-optimize)

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

    NormalBlockBase$new(data, q, sparsity = 0, control = NB_control())

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
      ll_list = NA
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

#### Returns

Update the current \[\`normal\`\] object

------------------------------------------------------------------------

### `NormalBlockBase$optimize()`

calls optimization (EM or heuristic) and updates relevant fields

#### Usage

    NormalBlockBase$optimize(control = list(niter = 100, threshold = 1e-04))

#### Arguments

- `control`:

  a list for controlling the optimization proces

#### Returns

optimizes the model and updates its parameters

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

    NormalBlockBase$candidates_split()

------------------------------------------------------------------------

### `NormalBlockBase$candidates_merge()`

generate and select a set of candidate models by merging the clusters of
the current model

#### Usage

    NormalBlockBase$candidates_merge()

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

plots log-likelihood values during model optimization

#### Usage

    NormalBlockBase$plot_loglik(type = "b", log = "xy", neg = TRUE)

#### Arguments

- `type`:

  char for line type (see plot.default)

- `log`:

  char for logarithmic axes (see plot.default)

- `neg`:

  boolean plot negative log-likelihood (useful when log="y")

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

plot together latent network and log-likelihood values during model
optimization

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
