# Control Parameters for Normal-Block Models

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

  maximal number of (V)EM iterations.

- threshold:

  the (V)EM stops when the objective (log-likelihood or ELBO) increases
  by less than this. A negative value runs exactly `niter` iterations,
  without the warning issued when `niter` is reached.

- fixed_point_niter:

  number of sweeps of the tau update for Normal-Block-Mean models with
  unknown clusters. Each sweep visits the rows of tau sequentially and
  maximizes the ELBO exactly in each, so it can never decrease the ELBO.

- sparsity_weights:

  a symmetric matrix of weights applied to the penalty, entry by entry:
  q x q (variance-block) or p x p (mean-block). By default, 1 off the
  diagonal and 0 on it (the diagonal is not penalized).

- sparsity_penalties:

  the penalties of a sparsity path (`sparsity = TRUE`). By default,
  `n_sparsity_penalties` values from the smallest penalty that removes
  every edge down to `min_ratio` times it, on a log scale.

- n_sparsity_penalties:

  number of penalties of the default path.

- min_ratio:

  ratio between the smallest and the largest penalty of the default
  path.

- fixed_tau:

  whether the posterior cluster probabilities stay at their initial
  value (`clustering_init`) during the optimization; used to refit a
  sparsity path with its clustering held fixed (stability selection).

- clustering_init:

  how to obtain the initial clustering of the q unknown blocks: a
  heuristic name ("ward2", "kmeans", "sbm" or "spectral"), an actual
  clustering (a vector of labels or a p x q indicator matrix, or a list
  of either per q for a collection), or "best_of_inits" to try several
  heuristics per model and keep the best-ELBO fit (see
  [NormalBlockBase](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)'s
  `best_of_inits()`; not supported with `sparsity = TRUE`). Default
  `NULL`, i.e. "kmeans" for both model families; "ward2" can work better
  on zero-inflated abundance data.

- verbose:

  whether to print the progress of the optimization.

- heuristic:

  whether to use the heuristic approach (moment-based, no (V)EM
  recursion) instead of the full (V)EM. Default is FALSE. In heuristic
  mode, no likelihood/ELBO is computed, so `entropy`, `loglik`, `BIC`,
  `ICL` and `EBIC` are all `NA` on the resulting model.

- noise_covariance:

  shape of the residual covariance. Variance-block models accept
  "diagonal" (variable-specific) or "spherical" (common); mean-block
  models, whose Sigma is the full p x p residual covariance, also accept
  "full". Default `NULL`, resolved per model family at fit time to
  "diagonal" – except for a mean-block model with `sparsity > 0`, which
  implies "full" since a penalty on a diagonal precision matrix would
  have nothing to act on (explicitly asking for both is an error). A
  full Sigma has many parameters and needs n \> p; use it when the
  residual associations between variables are of interest.

- refine:

  for collections over a range of cluster counts (both families):
  whether [`optimize()`](https://rdrr.io/r/stats/optimize.html)
  automatically calls `refine()` afterwards. Default `TRUE`: the split
  and merge search repairs the poor local optima that independent fits
  can leave, at a few times the cost of the collection. With `FALSE`,
  `collection$refine()` can be called later.

- acceleration:

  whether the (V)EM is accelerated by SQUAREM (Varadhan & Roland, 2008),
  an extrapolation built on two successive (V)EM steps and accepted only
  if it does not lower the objective. Default `NULL`: on, except for
  mean-block models with a diagonal or spherical Sigma, whose EM is
  exact and fast. `TRUE` or `FALSE` forces it. Never used when
  `sparsity > 0`.

## Value

A named list of parameters to pass to
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)'s
`control` argument.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
## silent fits, without the split and merge search
ctrl <- NB_control(verbose = FALSE, refine = FALSE)
collection <- normal_block(data, blocks = 2:4, control = ctrl)
summary(collection)
#> A collection of diagonal normal-block-var models over q 
#> ===========================================================================
#>     q ranging from 2 to 4 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -562.340 1124.681 1295.894 1295.894 1297.280    11
#>        48 3       3        0 -467.614  935.229 1121.173 1121.173 1127.765     7
#>        53 4       6        0 -467.618  935.235 1139.824 1139.824 1156.459    34
```
