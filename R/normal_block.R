#' Fit Normal-Block Models
#'
#' Fits a Normal-Block model, or a collection of them, by (variational) EM.
#' The arguments select the model: the family (`model`), whether the
#' clustering is known (`blocks`), whether the network is penalized
#' (`sparsity`) and whether Y is zero-inflated (`zero_inflation`).
#' @param data a [NormalBlockData] object, holding the responses (Y, n x p)
#' and the covariates (X, n x d).
#' @param blocks the clustering of the variables: a p x q indicator matrix
#' for a known clustering (e.g. `model.matrix(~ 0 + factor(labels))`), a
#' single integer q for an unknown clustering into q clusters, or a vector of
#' distinct integers for a collection over these values of q.
#' @param sparsity the graphical-lasso penalty on Omega: `0` (the default)
#' for none, a positive number for a single penalty, or `TRUE` for a
#' collection along a path of penalties (see [NB_control()]).
#' @param zero_inflation whether to model an excess of exact zeros in Y.
#' @param control a list of settings, built with [NB_control()].
#' @param model the model family: `"var"` (the default) clusters the
#' variables by how they covary, and infers a network between clusters (see
#' [NormalBlockVarBase]); `"mean"` clusters them by how they respond to the
#' covariates, `mu_i = C B' X_i` (see [NormalBlockMeanBase]). Every
#' combination of the other arguments is available for both families, except
#' a penalty together with zero-inflation for the mean-block family, whose
#' zero-inflated models have a diagonal or spherical Sigma.
#' @return a fitted model (an R6 object, see the classes listed in the
#' reference index), or a collection of fitted models.
#' @examples
#' ## Variance-block models over a range of q
#' ex_data <- generate_normal_block_var_data(n = 100, p = 30, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' my_normal_block <- normal_block(data, blocks = 1:5)
#' plot(my_normal_block)
#' best <- my_normal_block$get_best_model("ICL")
#' plot(ex_data$Y, fitted(best)); abline(0, 1)
#' ## the same with zero-inflation
#' ex_data_zi <- generate_normal_block_var_data(n = 50, p = 30, d = 1, q = 3,
#'                                              kappa = rep(0.5, 30))
#' zidata <- NormalBlockData$new(ex_data_zi$Y, ex_data_zi$X)
#' my_zi_block <- normal_block(zidata, blocks = 3, zero_inflation = TRUE)
#' ## a mean-block model (clustering in the mean rather than the covariance)
#' ex_mean <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' mean_data <- NormalBlockData$new(ex_mean$Y, ex_mean$X)
#' my_mean_block <- normal_block(mean_data, blocks = 3, model = "mean")
#' coef(my_mean_block)
#' @export
normal_block <- function(data,
                         blocks,
                         sparsity = 0,
                         zero_inflation = FALSE,
                         control = NB_control(),
                         model = c("var", "mean")) {
  ## Recovering the requested model from the function arguments
  model <- match.arg(model)
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  if (is.numeric(blocks) && !is.matrix(blocks) && length(blocks) > 1 && anyDuplicated(blocks))
    stop("`blocks` looks like a vector of cluster labels; a known clustering must be passed as a ", "p x q indicator matrix (e.g. model.matrix(~ 0 + factor(labels))). A numeric vector is ", "read as the range of cluster counts to explore, so its values must be distinct.", call. = FALSE)
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if (!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))
  if (is.list(control$clustering_init)) stopifnot(length(control$clustering_init) == length(blocks))
  stopifnot(
    "clustering_init = 'best_of_inits' is not supported together with sparsity = TRUE (a sparsity path warm-starts a single clustering across all penalties)" =
      !(uses_best_of_inits(control) && isTRUE(sparsity))
  )

  fit <- get_model(data, blocks, sparsity = sparsity,
                   zero_inflation = zero_inflation,
                   model = model,
                   control = control)

  ## Estimation/optimization
  if (control$verbose) cat("Fitting a", fit$who_am_I, "\n")

  ## a known clustering has no initialization to choose
  if (uses_best_of_inits(control) && isNB(fit) && !is.matrix(blocks)) {
    fit <- fit$best_of_inits(control = control)
  } else {
    fit$optimize(control)
  }

  ## Finishing
  if (control$verbose) cat("\nDONE\n")

  fit
}

#' Control Parameters for Normal-Block Models
#'
#' Control the model settings and various optimization parameters
#'
#' @param niter maximal number of (V)EM iterations.
#' @param fixed_point_niter number of sweeps of the tau update for
#' Normal-Block-Mean models with unknown clusters. Each sweep visits the rows
#' of tau sequentially and maximizes the ELBO exactly in each, so it can
#' never decrease the ELBO.
#' @param threshold the (V)EM stops when the objective (log-likelihood or
#' ELBO) increases by less than this. A negative value runs exactly `niter`
#' iterations, without the warning issued when `niter` is reached.
#' @param sparsity_weights a symmetric matrix of weights applied to the
#' penalty, entry by entry: q x q (variance-block) or p x p (mean-block). By
#' default, 1 off the diagonal and 0 on it (the diagonal is not penalized).
#' @param sparsity_penalties the penalties of a sparsity path (`sparsity =
#' TRUE`). By default, `n_sparsity_penalties` values from the smallest
#' penalty that removes every edge down to `min_ratio` times it, on a log
#' scale.
#' @param n_sparsity_penalties number of penalties of the default path.
#' @param min_ratio ratio between the smallest and the largest penalty of the
#' default path.
#' @param fixed_tau whether the posterior cluster probabilities stay at their
#' initial value (`clustering_init`) during the optimization; used to refit
#' a sparsity path with its clustering held fixed (stability selection).
#' @param clustering_init how to obtain the initial clustering of the q
#' unknown blocks: a heuristic name ("ward2", "kmeans", "sbm" or
#' "spectral"), an actual clustering (a vector of labels or a p x q indicator
#' matrix, or a list of either per q for a collection), or "best_of_inits" to
#' try several heuristics per model and keep the best-ELBO fit (see
#' [NormalBlockBase]'s `best_of_inits()`; not supported with `sparsity =
#' TRUE`). Default `NULL`, i.e. "kmeans" for both model families; "ward2"
#' can work better on zero-inflated abundance data.
#' @param verbose whether to print the progress of the optimization.
#' @param noise_covariance shape of the residual covariance. Variance-block
#' models accept "diagonal" (variable-specific) or "spherical" (common);
#' mean-block models, whose Sigma is the full p x p residual covariance, also
#' accept "full". Default `NULL`, resolved per model family at fit time to
#' "diagonal" -- except for a mean-block model with `sparsity > 0`, which
#' implies "full" since a penalty on a diagonal precision matrix would have
#' nothing to act on (explicitly asking for both is an error). A full Sigma
#' has many parameters and needs n > p; use it when the residual associations
#' between variables are of interest.
#' @param heuristic whether to use the heuristic approach (moment-based, no (V)EM
#' recursion) instead of the full (V)EM. Default is FALSE. In heuristic mode, no
#' likelihood/ELBO is computed, so `entropy`, `loglik`, `BIC`, `ICL` and `EBIC`
#' are all `NA` on the resulting model.
#' @param refine for collections over a range of cluster counts (both
#' families): whether `optimize()` automatically calls `refine()` afterwards.
#' Default `TRUE`: the split and merge search repairs the poor local optima
#' that independent fits can leave, at a few times the cost of the
#' collection. With `FALSE`, `collection$refine()` can be called later.
#' @param acceleration whether the (V)EM is accelerated by SQUAREM
#' (Varadhan & Roland, 2008), an extrapolation built on two successive (V)EM
#' steps and accepted only if it does not lower the objective. Default `NULL`:
#' on, except for mean-block models with a diagonal or spherical Sigma, whose
#' EM is exact and fast. `TRUE` or `FALSE` forces it. Never used when
#' `sparsity > 0`.
#' @return A named list of parameters to pass to [normal_block()]'s `control`
#' argument.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' ## silent fits, without the split and merge search
#' ctrl <- NB_control(verbose = FALSE, refine = FALSE)
#' collection <- normal_block(data, blocks = 2:4, control = ctrl)
#' summary(collection)
#' @export
NB_control <- function(
    niter                = 500,
    threshold            = 1e-4,
    fixed_point_niter    = 5,
    sparsity_weights     = NULL,
    sparsity_penalties   = NULL,
    n_sparsity_penalties = 30,
    min_ratio            = 0.01,
    fixed_tau            = FALSE,
    clustering_init      = NULL,
    verbose              = TRUE,
    heuristic            = FALSE,
    noise_covariance     = NULL,
    refine               = TRUE,
    acceleration         = NULL) {

  if (!is.null(sparsity_weights))
    stopifnot(all(is.matrix(sparsity_weights), isSymmetric(sparsity_weights)))
  if (!is.null(noise_covariance))
    stopifnot("noise_covariance must be one of 'diagonal', 'spherical' or (mean-block models only) 'full'" =
                length(noise_covariance) == 1 &&
                noise_covariance %in% c("diagonal", "spherical", "full"))
  if (is.character(clustering_init) && length(clustering_init) == 1) {
    stopifnot("clustering_init, when given as a single string, must name a known heuristic ('kmeans', 'ward2', 'sbm' or 'spectral'), or be 'best_of_inits' -- otherwise pass an actual clustering (a vector of labels, a p x q indicator matrix, or a list of either for a collection over several q values)" =
                clustering_init %in% c("kmeans", "ward2", "sbm", "spectral", "best_of_inits"))
  }

  structure(list(niter                = niter                ,
                 threshold            = threshold            ,
                 fixed_point_niter    = fixed_point_niter    ,
                 sparsity_weights     = sparsity_weights     ,
                 sparsity_penalties   = sparsity_penalties   ,
                 n_sparsity_penalties = n_sparsity_penalties ,
                 min_ratio            = min_ratio            ,
                 fixed_tau            = fixed_tau            ,
                 clustering_init      = clustering_init      ,
                 verbose              = verbose              ,
                 heuristic            = heuristic            ,
                 noise_covariance     = noise_covariance     ,
                 refine               = refine               ,
                 acceleration         = acceleration         ))
}


#' Create a Normal-Block Model Object
#'
#' Creates, without fitting it, the model (or collection of models) that
#' [normal_block()] would fit for the same arguments.
#' @inheritParams normal_block
#' @return an R6 object: a model or a collection of models, not optimized yet.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' ## the model is created but not fitted, unlike with normal_block()
#' model <- get_model(data, blocks = 3, control = NB_control(verbose = FALSE))
#' class(model)
#' @export
get_model <- function(data,
                      blocks,
                      sparsity = 0,
                      zero_inflation = FALSE,
                      control = NB_control(),
                      model = c("var", "mean")) {
  model <- match.arg(model)
  ## family default for scale = NULL (see NormalBlockData): resolved once
  ## here, so that every model of a collection shares the same data object
  data <- data$for_family(model)

  changing_sparsity <- isTRUE(sparsity)
  unknown_q_list    <- !is.matrix(blocks) && length(blocks) > 1
  family <- if (model == "var") "Var" else "Mean"

  if (changing_sparsity || unknown_q_list) {
    kind <- if (changing_sparsity && unknown_q_list) "CollectionClustersSparsity"
            else if (changing_sparsity) "CollectionSparsity" else "CollectionClusters"
    collection <- get(paste0("NormalBlock", family, kind))
    if (changing_sparsity) return(collection$new(data, blocks, zero_inflation, control = control))
    return(collection$new(data, blocks, zero_inflation, sparsity, control = control))
  }
  class_name <- paste0(if (zero_inflation) "ZI", "NormalBlock", family,
                       if (is.matrix(blocks)) "Known" else "Unknown", "Clusters")
  get(class_name)$new(data, blocks, sparsity, control = control)
}
