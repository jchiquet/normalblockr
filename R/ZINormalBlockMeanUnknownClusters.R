## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockMeanUnknownClusters #############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Normal-Block-Mean Model with Unknown Clustering
#'
#' R6 class for a zero-inflated Normal-Block-Mean model with a fixed number of
#' clusters (but unknown clustering), inferred by variational EM. Sigma is
#' diagonal or spherical here. See [NormalBlockMeanBase] for why a full one
#' is out of reach under a mask.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' Y <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
#' data <- NormalBlockData$new(Y, ex$X)
#' model <- normal_block(data, blocks = 3, model = "mean", zero_inflation = TRUE)
#' model$clustering
#' @export
ZINormalBlockMeanUnknownClusters <- R6::R6Class(
  classname = "ZINormalBlockMeanUnknownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockMeanUnknownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of clusters
    #' @param sparsity unused here, kept for signature symmetry (must be 0)
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINormalBlockMeanUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control, zero_inflation = TRUE)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## split(rule = "kmeans"): masked least-squares profiles, as the initialization
    split_target = function() self$data$X %*% self$data$zi_ols_fit()$B,

    ## Runs the VEM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockMeanUnknownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockMeanUnknownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
        noise_covariance = private$res_covariance,
        fixed_tau = isTRUE(control$fixed_tau),
        accelerate = private$use_acceleration(control),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, C = res$C, alpha = res$alpha,
           ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field var_par a list with the variational parameter: tau (posterior group probabilities)
    var_par = function() list(tau = private$C),
    #' @field model_par the parameters of [NormalBlockMeanBase] (B, Omega),
    #' plus B0 (zero-inflation coefficients), kappa (zero-inflation
    #' probabilities) and alpha (cluster proportions)
    model_par = function() c(super$model_par, list(B0 = private$B0, kappa = private$kappa, alpha = private$alpha))
  )
)
