## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionSparsity #########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block-Var Models over a Sparsity Path
#'
#' R6 class for a collection of normal-block models with a fixed clustering
#' (blocks) and different sparsity levels. The default penalty grid is read
#' off the latent block covariance of a short unpenalized fit.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = ex$parameters$C, sparsity = TRUE,
#'                        control = NB_control(verbose = FALSE, n_sparsity_penalties = 5))
#' models$get_best_model("BIC")$sparsity
#' @export
NormalBlockVarCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockVarCollectionSparsity",
  inherit   = NormalBlockCollectionSparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    family = "var",
    ## the latent block covariance of a short unpenalized fit
    penalty_scale = function(data, blocks, zero_inflation, control) {
      init_model <- get_model(data, blocks, 0, zero_inflation, control)
      init_model$optimize(control = list(niter = 5, threshold = 1e-4, verbose = FALSE), warn = FALSE)
      Sigmaq   <- chol2inv(chol(init_model$model_par$Omega))
      weights  <- init_model$sparsity_weights
      diag_pen <- max(diag(weights)) > 0
      abs((Sigmaq / weights)[upper.tri(Sigmaq, diag = diag_pen)])
    },
    ## StARS hooks (see NormalBlockCollectionSparsity): the clustering is held
    ## fixed on the subsamples -- the known one, or the variational tau of the
    ## least penalized model with fixed_tau = TRUE
    subsample_path = function(data, control) {
      if (!is.matrix(private$blocks_)) {
        control$clustering_init <- self$get_model(min(private$sparsity_))$var_par$tau
        control$fixed_tau <- TRUE
      }
      NormalBlockVarCollectionSparsity$new(data, private$blocks_, isTRUE(self$control$zero_inflation), control)
    },
    network_nodes = function() cluster_names(self$q)
  )
)
