## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionSparsity #########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block-Var Models over a Sparsity Path
#'
#' R6 class for a collection of normal-block models with a fixed clustering
#' (blocks) and different sparsity levels.
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
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarCollectionSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param blocks either a clustering matrix (known, fixed clustering) or a single integer (number of blocks to infer)
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockVarCollectionSparsity`] object
    initialize = function(mydata, blocks, zero_inflation = FALSE,
                          control = NB_control()) {

      ## Store user-defined fields
      self$data    <- mydata
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      private$blocks_ <- blocks

      ## Check block format
      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)

      ## extract a collection of sparsifying penalties
      if (!is.null(control$sparsity_penalties)){
        stopifnot("All penalties must be strictly positive" =
                    (min(control$sparsity_penalties) > 0))
        sparsity <- control$sparsity_penalties
      } else {
        init_model <- get_model(mydata, blocks, 0, zero_inflation, control)
        init_model$optimize(control = list(niter=5, threshold=1e-4, verbose=FALSE), warn = FALSE)
        Sigmaq   <- chol2inv(chol(init_model$model_par$Omega))
        diag_pen <- max(diag(init_model$sparsity_weights)) > 0
        weights  <- init_model$sparsity_weights
        weights  <- abs((Sigmaq / weights)[upper.tri(Sigmaq, diag = diag_pen)])
        if (length(weights) > 0) {
          max_pen  <- max(weights)
          sparsity <- 10^seq(log10(max_pen), log10(max_pen * control$min_ratio), len = control$n_sparsity_penalties)
        } else {
          sparsity <- rep(0, len = control$n_sparsity_penalties)
        }
      }

      private$sparsity_ <- sort(sparsity, decreasing = TRUE)
      ## Instantiation of the models in the collection
      self$models <- map(private$sparsity_, function(lambda)
        model <- get_model(mydata, blocks, lambda, zero_inflation, control)
      )
      private$progress_field <- "sparsity"
      private$progress_label <- "penalty"
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
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
    network_nodes = function() paste0("group_", seq_len(self$q))
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field sparsity_details list of information about model's penalties
    sparsity_details = function()
      list("n_penalties" = length(self$sparsity), "min_ratio" = min(self$sparsity)/max(self$sparsity),
           "min_sparsity" = min(self$sparsity), "max_sparsity" = max(self$sparsity)),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){
      nc <- if (is.null(self$control$noise_covariance)) "diagonal" else self$control$noise_covariance
      zi <- if (isTRUE(self$control$zero_inflation)) "zero-inflated " else ""
      paste0("Collection of ", zi, nc, " normal-block-var models with ",
             ifelse(is.matrix(private$blocks_), "fixed blocks", "fixed q"),
        ", with different sparsity penalties.")}
  )
)
