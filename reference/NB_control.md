# NB_control

Control the model settings and various optimization parameters

## Usage

``` r
NB_control(
  niter = 100,
  threshold = 1e-04,
  sparsity_weights = NULL,
  sparsity_penalties = NULL,
  n_sparsity_penalties = 30,
  min_ratio = 0.01,
  fixed_tau = FALSE,
  clustering_init = NULL,
  verbose = TRUE,
  heuristic = FALSE,
  noise_covariance = c("diagonal", "spherical"),
  clustering_approx = c("ward2", "kmeans", "sbm")
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

  proposal of initial value for clustering, when q is unknown, can be a
  list with one clustering for each q value

- verbose:

  telling if information should be printed during optimization

- heuristic:

  weather to use heuristic approach or not. Default is FALSE

- noise_covariance:

  variance can be variable specific ("diagonal", the default) or common
  ("spherical")

- clustering_approx:

  to use for clustering with heuristic inference method
