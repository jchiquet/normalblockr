## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS ZINormalBlockMeanKnownClusters ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Zero-Inflated Normal-Block-Mean Model with Known Clustering
#'
#' R6 class for a zero-inflated Normal-Block-Mean model with a known
#' clustering. Sigma is diagonal or spherical here. See
#' [NormalBlockMeanBase] for why a full one is out of reach under a mask.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' Y <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
#' data <- NormalBlockData$new(Y, ex$X)
#' model <- normal_block(data, blocks = ex$parameters$C, model = "mean",
#'                       zero_inflation = TRUE)
#' model$clustering
#' @export
ZINormalBlockMeanKnownClusters <- R6::R6Class(
  classname = "ZINormalBlockMeanKnownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`ZINormalBlockMeanKnownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param C p x q clustering matrix, C_jk = 1 if variable j belongs to cluster k
    #' @param sparsity unused here, kept for signature symmetry (must be 0)
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`ZINormalBlockMeanKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      check_known_clustering(C, data$p)
      super$initialize(data, ncol(C), sparsity, control, zero_inflation = TRUE)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    known_clusters = TRUE,
    ## split(rule = "kmeans"): masked least-squares profiles, as the initialization
    split_target = function() self$data$X %*% self$data$zi_ols_fit()$B,

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## ZINormalBlockMeanKnownClusters_fit).
    EM_optimize = function(control) {
      init <- private$optim_initialize()
      res  <- ZINormalBlockMeanKnownClusters_fit(
        Y = self$data$Y, X = self$data$X,
        zeros_bar = self$data$zeros_bar, zi_cond_mean = private$ZI_cond_mean,
        C = private$C, B0 = init$B, Omega0 = init$Omega,
        noise_covariance = private$res_covariance,
        accelerate = private$use_acceleration(control),
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, Omega = res$Omega, ll_list = res$objective)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par the parameters of [NormalBlockMeanBase] (B, Omega),
    #' plus B0 (zero-inflation coefficients) and kappa (zero-inflation
    #' probabilities)
    model_par = function() c(super$model_par, list(B0 = private$B0, kappa = private$kappa))
  )
)
