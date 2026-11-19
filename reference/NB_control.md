# NB_control

Control the model settings and various optimization parameters

## Usage

``` r
NB_control(
  niter = 500,
  threshold = 1e-04,
  fixed_point_niter = 5,
  sparsity_weights = NULL,
  sparsity_penalties = NULL,
  n_sparsity_penalties = 30,
  min_ratio = 0.01,
  fixed_tau = FALSE,
  clustering_init = NULL,
  verbose = TRUE,
  heuristic = FALSE,
  noise_covariance = NULL,
  refine = TRUE,
  acceleration = NULL
)
```

## Arguments

- niter:

  number of iterations in model optimization

- threshold:

  loglikelihood / elbo threshold under which optimization stops

- fixed_point_niter:

  number of sweeps of the tau update for Normal-Block-Mean with unknown
  clusters. Each sweep visits the rows of tau sequentially and maximizes
  the ELBO exactly in each, so it can never decrease the ELBO.

- sparsity_weights:

  weights with which the penalty should be applied in case sparsity is
  required, non-0 values on the diagonal mean diagonal shall be
  penalized too (default is non-penalized diagonal and 1s off-diagonal)

- sparsity_penalties:

  list of penalties the user wants to test, other parameters are only
  used if penalties is not specified

- n_sparsity_penalties:

  number of penalties to test.

- min_ratio:

  ratio for sparsity between max penalty (0 edge penalty) and min
  penalty to test

- fixed_tau:

  whether tau should be fixed at clustering_init during optimization
  useful for calls to fixed_q models in stability_selection

- clustering_init:

  how to obtain the initial clustering of the q unknown blocks: a
  heuristic name ("ward2", "kmeans", "sbm" or "spectral"), an actual
  clustering (a vector of labels or a p x q indicator matrix, or a list
  of either per q for a collection), or "best_of_inits" to try several
  heuristics per model and keep the best-ELBO fit (see
  \[NormalBlockVarBase\]'s \`best_of_inits()\`; not supported with
  \`sparsity = TRUE\`). Default \`NULL\`, i.e. "kmeans" for both model
  families; "ward2" can work better on zero-inflated abundance data.

- verbose:

  telling if information should be printed during optimization

- heuristic:

  whether to use the heuristic approach (moment-based, no (V)EM
  recursion) instead of the full (V)EM. Default is FALSE. In heuristic
  mode, no likelihood/ELBO is computed, so \`entropy\`, \`loglik\`,
  \`BIC\`, \`ICL\` and \`EBIC\` are all \`NA\` on the resulting model.

- noise_covariance:

  shape of the residual covariance. Variance-block models accept
  "diagonal" (variable-specific) or "spherical" (common); mean-block
  models, whose Sigma is the full p x p residual covariance, also accept
  "full". Default \`NULL\`, resolved per model family at fit time to
  "diagonal" – except for a mean-block model with \`sparsity \> 0\`,
  which implies "full" since a penalty on a diagonal precision matrix
  would have nothing to act on (explicitly asking for both is an error).
  A full Sigma has many parameters and needs n \> p; use it when the
  residual associations between variables are of interest.

- refine:

  for collections over a range of cluster counts (both families):
  whether \`optimize()\` automatically calls \`refine()\` afterwards.
  Default \`TRUE\`: the split and merge search repairs the poor local
  optima that independent fits can leave, at a few times the cost of the
  collection. With \`FALSE\`, \`collection\$refine()\` can be called
  later.

- acceleration:

  whether the (V)EM is accelerated by SQUAREM (Varadhan & Roland, 2008),
  an extrapolation built on two successive (V)EM steps and accepted only
  if it does not lower the objective. Default \`NULL\`: on, except for
  mean-block models with a diagonal or spherical Sigma, whose EM is
  exact and fast. \`TRUE\` or \`FALSE\` forces it. Never used when
  \`sparsity \> 0\`.

## Value

A named list of parameters to pass to \[normal_block()\]'s \`control\`
argument.
