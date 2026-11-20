## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockVarUnknownClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Normal-Block-Var Model with Unknown Clustering
#'
#' R6 class for a zero-inflated normal-block model with a fixed number of
#' clusters (but unknown clustering).
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3, kappa = rep(0.3, 20))
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = 3, zero_inflation = TRUE,
#'                       control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
ZINormalBlockVarUnknownClusters <- R6::R6Class(
  classname = "ZINormalBlockVarUnknownClusters",
  inherit   = NormalBlockVarBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockVarUnknownClusters`] object.
    #' @param data object of NormalBlockVarData class, with responses and design matrix
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param q required number of groups
    #' @param control structured list of more specific parameters
    #' @return A new [`ZINormalBlockVarUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control, zero_inflation = TRUE)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS -------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## split(rule = "kmeans"): residuals of the non-zero-inflated cells only
    split_target = function() (self$data$Y - self$data$X %*% private$B) * self$data$zeros_bar,

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function() {
      zi_diag <- self$data$zi_ols_fit()
      if (anyNA(private$C))
        private$C <- private$heuristic_clustering(zi_diag$R)
      private$C <- clip_probabilities(private$C)
      Sigmaq <- private$heuristic_Sigmaq_from_Sigma(cov(zi_diag$R))
      Omega <- private$get_Omega(Sigmaq)
      list(B = zi_diag$B, dm1 = private$zi_dm1(), Omega = Omega,
           alpha = colMeans(private$C), C = private$C)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, dm1 = private$dm1, Omega = private$Omega,
             alpha = private$alpha, C = private$C, M = private$M, S = private$S)
      } else c(private$get_heuristic_parameters(),  list(
          M = matrix(rep(0, self$n * self$q), nrow = self$n),
          S = matrix(rep(0.1, self$n * self$q), nrow = self$n)
        )
      )
    },

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockVarUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockVarUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance, fixed_tau = isTRUE(control$fixed_tau),
        accelerate = !isFALSE(control$acceleration),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, dm1 = res$dm1, Omega = res$Omega, alpha = res$alpha,
           C = res$C, M = res$M, S = res$S, ll_list = res$objective)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field var_par a list with variational parameters
    var_par  = function() {list(M = private$M, S = private$S, tau = private$C)},
    #' @field model_par the parameters of [NormalBlockVarBase] (B, B0, dm1, Omega), plus kappa (zero-inflation probabilities) and alpha (cluster proportions)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par$alpha <- private$alpha
      par
    }
  )
)
