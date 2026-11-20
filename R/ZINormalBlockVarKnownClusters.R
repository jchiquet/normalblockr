## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockVarKnownClusters #############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Normal-Block-Var Model with Known Clustering
#'
#' R6 class for a zero-inflated normal-block model with a known clustering.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3, kappa = rep(0.3, 20))
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' model <- normal_block(data, blocks = ex$parameters$C, zero_inflation = TRUE,
#'                       control = NB_control(verbose = FALSE))
#' model$clustering
#' @export
ZINormalBlockVarKnownClusters <- R6::R6Class(
  classname = "ZINormalBlockVarKnownClusters",
  inherit   = NormalBlockVarBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockVarKnownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param C p x q clustering matrix, C_jk = 1 if variable j belongs to cluster k
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINormalBlockVarKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      check_known_clustering(C, data$p)
      super$initialize(data, ncol(C), sparsity, control = control, zero_inflation = TRUE)
      private$C  <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    known_clusters = TRUE,
    ## split(rule = "kmeans"): residuals of the non-zero-inflated cells only
    split_target = function() (self$data$Y - self$data$X %*% private$B) * self$data$zeros_bar,

    get_heuristic_parameters = function(){
      zi_diag <- self$data$zi_ols_fit()
      Sigmaq  <- private$heuristic_Sigmaq_from_Sigma(cov(zi_diag$R))
      Omega  <- private$get_Omega(Sigmaq)
      list(B = zi_diag$B, dm1 = private$zi_dm1(), Omega = Omega)
    },

    optim_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, dm1 = private$dm1, Omega = private$Omega,
             gamma = private$gamma, mu = private$mu)
      } else c(private$get_heuristic_parameters(),  list(
        gamma = rep(list(diag(1, self$q, self$q)), self$n),
        mu    = matrix(0, self$n, self$q)
        )
      )
    },

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockVarKnownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockVarKnownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean, C = private$C,
        B0 = init$B, dm1_0 = init$dm1, Omega0 = init$Omega,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance,
        accelerate = !isFALSE(control$acceleration),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, dm1 = res$dm1, Omega = res$Omega,
           gamma = purrr::array_branch(res$gamma, 3), mu = res$mu, ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function() list(gamma = private$gamma, mu = private$mu),
    #' @field model_par the parameters of [NormalBlockVarBase] (B, B0, dm1, Omega), plus kappa (zero-inflation probabilities)
    model_par  = function() {
      par       <- super$model_par
      par$kappa <- private$kappa
      par
    }
  )
)
