## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionSparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Normal-Block-Mean Models over a Sparsity Path
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase]) with
#' a fixed clustering (or a fixed number of blocks) and different sparsity
#' levels applied to the p x p precision matrix of the variables. Mirrors
#' [NormalBlockVarCollectionSparsity], StARS included; the default penalty
#' grid is read off the covariance of the least-squares residuals.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 3, sparsity = TRUE, model = "mean",
#'                        control = NB_control(n_sparsity_penalties = 5))
#' plot(models, criteria = c("BIC", "EBIC"))
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
    #' @param zero_inflation must be `FALSE`: the path penalizes a full Sigma,
    #' which zero-inflated mean-block models do not carry
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockMeanCollectionSparsity`] object
    initialize = function(mydata, blocks, zero_inflation = FALSE, control = NB_control()) {
      stopifnot(
        "zero-inflation is not available for a mean-block sparsity path: it needs a full Sigma, which zero-inflated mean-block models do not carry" =
          !zero_inflation
      )
      super$initialize(mydata, blocks, zero_inflation, control)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    family = "mean",
    ## the covariance of the least-squares residuals
    penalty_scale = function(data, blocks, zero_inflation, control) {
      Sigma   <- data$ols_fit()$Sigma
      weights <- matrix(1, data$p, data$p)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) weights <- control$sparsity_weights
      diag_pen <- max(diag(weights)) > 0
      abs((Sigma / weights)[upper.tri(Sigma, diag = diag_pen)])
    },
    ## StARS hooks (see NormalBlockCollectionSparsity): the subsample paths
    ## hold the clustering fixed, as a known one -- the given blocks, or the
    ## MAP partition of the least penalized model (empty clusters dropped)
    subsample_path = function(data, control) {
      blocks <- private$blocks_
      if (!is.matrix(blocks)) {
        cl <- self$get_model(min(private$sparsity_))$clustering
        blocks <- as_indicator(match(cl, sort(unique(cl))))
      }
      NormalBlockMeanCollectionSparsity$new(data, blocks, control = control)
    },
    network_nodes = function() {
      nm <- colnames(self$data$Y)
      if (is.null(nm)) paste0("Y", seq_len(self$data$p)) else nm
    }
  )
)
