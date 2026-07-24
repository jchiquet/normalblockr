# NB_control

Control the model settings and various optimization parameters

## Usage

``` r
NB_control(
  niter = 500,
  threshold = 1e-04,
  sparsity_weights = NULL,
  sparsity_penalties = NULL,
  n_sparsity_penalties = 30,
  min_ratio = 0.01,
  fixed_tau = FALSE,
  clustering_init = "ward2",
  verbose = TRUE,
  heuristic = FALSE,
  noise_covariance = c("diagonal", "spherical"),
  refine = FALSE
)
```

## Arguments

- niter:

  number of iterations in model optimization

- threshold:

  loglikelihood / elbo threshold under which optimization stops

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
  heuristic name ("ward2", the default, "kmeans", "sbm" or "spectral"),
  an actual clustering (a vector of labels or a p x q indicator matrix,
  or a list of either per q for a collection), or "best_of_inits" to try
  several heuristics per model and keep the best-ELBO fit (see
  \[NormalBlockVarBase\]'s \`best_of_inits()\`; not supported with
  \`sparsity = TRUE\`). See
  \`inst/methods_initialization_and_refine.md\` for the heuristics'
  rationale, why no single one dominates, and how this interacts with
  \`refine\` (below).

- verbose:

  telling if information should be printed during optimization

- heuristic:

  whether to use the heuristic approach (moment-based, no (V)EM
  recursion) instead of the full (V)EM. Default is FALSE. In heuristic
  mode, no likelihood/ELBO is computed, so \`entropy\`, \`loglik\`,
  \`BIC\`, \`ICL\` and \`EBIC\` are all \`NA\` on the resulting model.

- noise_covariance:

  variance can be variable specific ("diagonal", the default) or common
  ("spherical")

- refine:

  for \[NormalBlockVarCollectionClusters\] only: whether \`optimize()\`
  should automatically call \`refine()\` afterwards. Default \`FALSE\`
  since it adds real cost; call \`collection\$refine()\` directly at any
  point afterwards for the same effect without setting this.

## Value

A named list of parameters to pass to \[normal_block()\]'s \`control\`
argument.
