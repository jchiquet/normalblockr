## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarUnknownClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Normal-Block-Var Model with Unknown Clustering
#'
#' R6 class for a normal-block model with a fixed number of clusters (but
#' unknown clustering).
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
NormalBlockVarUnknownClusters <- R6::R6Class(
  classname = "NormalBlockVarUnknownClusters",
  inherit = NormalBlockVarBase,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarUnknownClusters`] object.
    #' @param data contains the matrix of responses (Y) and the design matrix (X).
    #' @param q required number of groups
    #' @param sparsity sparsity penalty to add on blocks precision matrix for sparsity
    #' @param control structured list for specific parameters
    #' @return A new [`NormalBlockVarUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      if (anyNA(private$C)) # if no initial clustering provided
        private$C <- private$heuristic_clustering(reg_res$R)
      private$C <- clip_probabilities(private$C)
      Sigmaq    <- private$heuristic_Sigmaq_from_Sigma(reg_res$Sigma)
      Omega    <- private$get_Omega(Sigmaq)
      dm1       <- private$dm1_from_residuals(reg_res$R)
      list(B = reg_res$B, Omega = Omega, dm1 = dm1,
           C = private$C, alpha = colMeans(private$C))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, Omega = private$Omega, dm1 = private$dm1,
             C = private$C, alpha = private$alpha, M = private$M, S = private$S)
      } else c(private$get_heuristic_parameters(),  list(
            M = matrix(rep(0, self$n * self$q), nrow = self$n),
            S = rep(0.1, self$q)
          )
      )
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockVarUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- NormalBlockVarUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance, fixed_tau = isTRUE(control$fixed_tau),
        accelerate = !isFALSE(control$acceleration),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, dm1 = res$dm1, alpha = res$alpha,
           C = res$C, M = res$M, S = res$S, ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par the parameters of [NormalBlockVarBase] (B, B0, dm1, Omega), plus alpha (cluster proportions)
    model_par  = function() {
      parameters       <- super$model_par
      parameters$alpha <- private$alpha
      parameters},
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$C)
  )
)
