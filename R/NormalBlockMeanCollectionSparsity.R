## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionSparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Mean-Block Models over a Sparsity Path
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase]) with
#' a fixed clustering (or a fixed number of blocks) and different sparsity
#' levels applied to the p x p precision matrix of the variables. Mirrors
#' [NormalBlockVarCollectionSparsity], StARS included: the network relates
#' the variables, so its edges do not depend on cluster labels, and the
#' subsample fits simply hold the clustering fixed (the known one, or the MAP
#' partition of the least penalized model).
#' @examples
#' ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 3, sparsity = TRUE, model = "mean",
#'                        control = NB_control(n_sparsity_penalties = 5))
#' models$plot(c("BIC", "EBIC"))
#' @export
NormalBlockMeanCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockMeanCollectionSparsity",
  inherit   = NormalBlockCollectionSparsity,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanCollectionSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param blocks either a clustering matrix (known, fixed clustering) or a single integer (number of blocks to infer)
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockMeanCollectionSparsity`] object
    initialize = function(mydata, blocks, control = NB_control()) {
      self$data    <- mydata
      self$control <- control
      private$blocks_ <- blocks

      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)

      if (!is.null(control$sparsity_penalties)) {
        stopifnot("All penalties must be strictly positive" =
                    (min(control$sparsity_penalties) > 0))
        sparsity <- control$sparsity_penalties
      } else {
        ## The variance-block path reads its scale off a short unpenalized fit.
        Sigma    <- stats::cov(ols_residuals(mydata))
        weights  <- matrix(1, mydata$p, mydata$p)
        diag(weights) <- 0
        if (!is.null(control$sparsity_weights)) weights <- control$sparsity_weights
        diag_pen <- max(diag(weights)) > 0
        scale    <- abs((Sigma / weights)[upper.tri(Sigma, diag = diag_pen)])
        scale    <- scale[is.finite(scale)]
        if (length(scale) > 0) {
          max_pen  <- max(scale)
          sparsity <- 10^seq(log10(max_pen), log10(max_pen * control$min_ratio),
                             len = control$n_sparsity_penalties)
        } else {
          sparsity <- rep(0, len = control$n_sparsity_penalties)
        }
      }

      private$sparsity_ <- sort(sparsity, decreasing = TRUE)
      self$models <- map(private$sparsity_, function(lambda)
        get_model(mydata, blocks, sparsity = lambda, model = "mean", control = control)
      )
      private$progress_field <- "sparsity"
      private$progress_label <- "penalty"
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## StARS hooks (see NormalBlockCollectionSparsity): the subsample paths
    ## hold the clustering fixed, as a known one -- the given blocks, or the
    ## MAP partition of the least penalized model (empty clusters dropped)
    subsample_path = function(data, control) {
      blocks <- private$blocks_
      if (!is.matrix(blocks)) {
        cl <- self$get_model(min(private$sparsity_))$clustering
        blocks <- as_indicator(match(cl, sort(unique(cl))))
      }
      NormalBlockMeanCollectionSparsity$new(data, blocks, control)
    },
    network_nodes = function() {
      nm <- colnames(self$data$Y)
      if (is.null(nm)) paste0("Y", seq_len(self$data$p)) else nm
    }
  ),
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function() "normal-block-mean model with sparsity path"
  )
)
