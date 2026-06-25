#' Normal-block model
#'
#' Fit a normal-block model with a variational or heuristic algorithm
#' @param data NormalBlockData object, contains the matrix of responses (Y, n x p) and the design matrix (X, n x d), must be created with NormalBlockData$new.
#' @param blocks either a integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * q matrix (for indicating block membership when its known)
#' @param sparsity either TRUE to run the optimization for different sparsity penalty values
#' OR float to run model with a single sparsity penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param control a list-like structure for detailed control on parameters should be
#' generated with NB_control().
#' @return an R6 object with one of the model classes (or a collection of model objects).
#' @examples
#' ## Normal Data
#' ex_data <- generate_normal_block_data(n=50, p=50, d=1, q=3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' my_normal_block <- normal_block(data, blocks = 1:6)
#' \dontrun{
#' my_normal_block$plot(c("deviance", "BIC", "ICL"))
#' Y_hat <- my_normal_block$get_best_model()$fitted
#' plot(data$Y, Y_hat, log = "xy"); abline(0,1)
#' }
#' ## Normal Data with Zero Inflation
#' ex_data_zi <- generate_normal_block_data(n=50, p=50, d=1, q=3, kappa = rep(0.5,50))
#' zidata <- NormalBlockData$new(ex_data_zi$Y, ex_data_zi$X)
#' my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)
#'
#' @export
normal_block <- function(data,
                         blocks,
                         sparsity = 0,
                         zero_inflation = FALSE,
                         control = NB_control()) {
  ## Recovering the requested model from the function arguments
  stopifnot(is.numeric(blocks) | is.matrix(blocks))
  stopifnot(is.null(control$sparsity_weights) | is.matrix(control$sparsity_weights))
  if (!is.null(control$sparsity_weights)) stopifnot(isSymmetric(control$sparsity_weights))
  if (is.list(control$clustering_init)) stopifnot(length(control$clustering_init) == length(blocks))

  model <- get_model(data, blocks, sparsity = sparsity,
                     zero_inflation = zero_inflation,
                     control = control)

  ## Estimation/optimization
  if (control$verbose) cat("Fitting a", model$who_am_I, "\n")

  model$optimize(control)

  ## Finishing
  if (control$verbose) cat("\nDONE\n")

  model
}

#' NB_control
#'
#' Control the model settings and various optimization parameters
#'
#' @param niter number of iterations in model optimization
#' @param threshold loglikelihood / elbo threshold under which optimization stops
#' @param sparsity_weights weights with which the penalty should be applied in case
#' sparsity is required, non-0 values on the diagonal mean diagonal shall be
#' penalized too (default is non-penalized diagonal and 1s off-diagonal)
#' @param sparsity_penalties list of penalties the user wants to test, other parameters
#' are only used if penalties is not specified
#' @param n_sparsity_penalties number of penalties to test.
#' @param min_ratio ratio for sparsity between max penalty (0 edge penalty) and min penalty to test
#' @param fixed_tau whether tau should be fixed at clustering_init during optimization
#' useful for calls to fixed_q models in stability_selection
#' @param clustering_init how to obtain the initial clustering of the q unknown
#' blocks. Either the name of a clustering heuristic -- one of "ward2"
#' (default), "kmeans", "sbm" or "spectral" (a cheap proxy for "sbm": k-means
#' on the row-normalized eigenvectors of cov(R)) -- or an actual clustering to
#' use directly, as a vector of labels or a p x q indicator matrix. When q is
#' unknown (a collection over several q values), can also be a list with one
#' such heuristic name/clustering per q value. No single heuristic dominates
#' on every dataset (see inst/clustering_initialization_benchmark); "ward2" is
#' the default for giving the best balance of BIC rank and how rarely its
#' deviance path violates the model's "non-increasing in q" guarantee -- a
#' reliability signal "kmeans" lacks despite a marginally better raw rank.
#' With the heuristic name "sbm" on a collection (and no per-q list of
#' explicit clusterings), a single SBM exploration runs over the whole range
#' of q and is reused for every model instead of repeating one per q; any q
#' it doesn't reach falls back to a cheap "ward2" clustering.
#' @param verbose telling if information should be printed during optimization
#' @param noise_covariance variance can be variable specific ("diagonal", the default) or common ("spherical")
#' @param heuristic whether to use the heuristic approach (moment-based, no (V)EM
#' recursion) instead of the full (V)EM. Default is FALSE. In heuristic mode, no
#' likelihood/ELBO is computed, so `entropy`, `loglik`, `BIC`, `ICL` and `EBIC`
#' are all `NA` on the resulting model.
#' @param refine for [NormalBlockCollectionClusters] only: whether
#' `optimize()` should automatically call `refine()` afterwards (see its
#' documentation for the rationale and an empirical before/after comparison).
#' Default `FALSE` -- it adds real cost (roughly 3x the time of fitting the
#' collection alone), so it is opt-in; call `collection$refine()` directly at
#' any point afterwards for the same effect without setting this.
#' @export
NB_control <- function(
    niter                = 100,
    threshold            = 1e-4,
    sparsity_weights     = NULL,
    sparsity_penalties   = NULL,
    n_sparsity_penalties = 30,
    min_ratio            = 0.01,
    fixed_tau            = FALSE,
    clustering_init      = "ward2",
    verbose              = TRUE,
    heuristic            = FALSE,
    noise_covariance     = c("diagonal", "spherical"),
    refine               = FALSE) {

  if (!is.null(sparsity_weights))
    stopifnot(all(is.matrix(sparsity_weights), isSymmetric(sparsity_weights)))
  if (is.character(clustering_init) && length(clustering_init) == 1) {
    stopifnot("clustering_init, when given as a single string, must name a known heuristic ('kmeans', 'ward2', 'sbm' or 'spectral') -- otherwise pass an actual clustering (a vector of labels, a p x q indicator matrix, or a list of either for a collection over several q values)" =
                clustering_init %in% c("kmeans", "ward2", "sbm", "spectral"))
  }

  structure(list(niter                = niter                ,
                 threshold            = threshold            ,
                 sparsity_weights     = sparsity_weights     ,
                 sparsity_penalties   = sparsity_penalties   ,
                 n_sparsity_penalties = n_sparsity_penalties ,
                 min_ratio            = min_ratio            ,
                 fixed_tau            = fixed_tau            ,
                 clustering_init      = clustering_init      ,
                 verbose              = verbose              ,
                 heuristic            = heuristic            ,
                 noise_covariance     = match.arg(noise_covariance),
                 refine               = refine               ))
}


#' Creates appropriate new normal block model depending on the parametrization
#' @param blocks either an integer (number of blocks), a vector of integer (list of possible number of block)
#'  or a p * q matrix (for indicating block membership when its known)
#' @param sparsity boolean to say whether the model should have a changing penalty
#' OR float to run model with a single penalty value
#' @param zero_inflation boolean to indicate if Y is zero-inflated and adjust fitted model as a consequence
#' @param control a list-like structure for detailed control on parameters should be
#' generated with normal_block_control() for collections of sparse models
#' @param data contains the matrix of responses (Y) and the design matrix (X).
#' @export
get_model <- function(data,
                      blocks,
                      sparsity = 0,
                      zero_inflation = FALSE,
                      control = NB_control()) {

  changing_sparsity <- isTRUE(sparsity)
  unknown_q_list    <- !is.matrix(blocks) && length(blocks) > 1
  is_collection     <- changing_sparsity || unknown_q_list

  if (is_collection) {
    class_name <- if (changing_sparsity && unknown_q_list) {
      "NormalBlockCollectionClustersSparsity"
    } else if (changing_sparsity) {
      "NormalBlockCollectionSparsity"
    } else {
      "NormalBlockCollectionClusters"
    }
  } else {
    class_name <- if (is.matrix(blocks)) "NormalBlockKnownClusters" else "NormalBlockUnknownClusters"
    if (zero_inflation) class_name <- paste0("ZI", class_name)
  }

  myClass <- eval(str2lang(class_name))
  if (is_collection) {
    if (changing_sparsity) {
      model <- myClass$new(data, blocks, zero_inflation, control = control)
    } else {
      model <- myClass$new(data, blocks, zero_inflation, sparsity, control = control)
    }
  } else {
    model <- myClass$new(data, blocks, sparsity, control = control)
  }
  model
}
