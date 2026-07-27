## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanBase ##########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a generic sparse Normal Block model
NormalBlockMeanBase <- R6::R6Class(
  classname = "NormalBlockMeanBase",
  inherit = NormalBlockBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NormalBlockMeanBase`] object.
    #' @param data object of NormalMeanBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_Mean_control
    #' @return A new [`NormalBlockMeanBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, self$data$p, self$data$p)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    Psi               = NA,
    Phi               = NA,
    Lambda            = NA,

    heuristic_cluster_B_from_variable_B = function(B_variable, C){
      B <- B_variable %*% C / rep(colSums(C), each = nrow(B_variable))
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalMeanBlockData(scale = TRUE)`) -- use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, Omega = private$Omega),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      as.integer(self$q * self$d + self$p + self$n_edges)
    },
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a p x p matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$p))
        private$weights <- value
      }
    }
)
)
