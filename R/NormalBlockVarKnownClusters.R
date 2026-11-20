## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarKnownClusters #############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Normal-Block-Var Model with Known Clustering
#'
#' R6 class for a normal-block model with known clustering.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = ex$parameters$C, control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
NormalBlockVarKnownClusters <- R6::R6Class(
  classname = "NormalBlockVarKnownClusters",
  inherit   = NormalBlockVarBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarKnownClusters`] object.
    #' @param data object of NormalBlockVarData class, with responses and design matrix
    #' @param C p x q clustering matrix, C_jk = 1 if variable j belongs to cluster k
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockVarKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      check_known_clustering(C, data$p)
      super$initialize(data, ncol(C), sparsity, control)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    known_clusters = TRUE,

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      dm1       <- private$dm1_from_residuals(reg_res$R)
      Sigmaq    <- private$heuristic_Sigmaq_from_Sigma(reg_res$Sigma)
      Omega    <- private$get_Omega(Sigmaq)
      list(B = reg_res$B, dm1 = dm1, Omega = Omega)
    },

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, dm1 = private$dm1, Omega = private$Omega,
             gamma = private$gamma, mu = private$mu)
      } else c(private$get_heuristic_parameters(),  list(
        gamma = diag(1, self$q, self$q),
        mu    = matrix(0, self$n, self$q)
        )
      )
    },

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockVarKnownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- NormalBlockVarKnownClusters_fit(
        Y = self$data$Y, X = self$data$X, C = private$C,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance,
        accelerate = !isFALSE(control$acceleration),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, dm1 = res$dm1, Omega = res$Omega,
           gamma = res$gamma, mu = res$mu, ll_list = res$objective)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function() list(gamma = private$gamma, mu = private$mu)
  )
)
