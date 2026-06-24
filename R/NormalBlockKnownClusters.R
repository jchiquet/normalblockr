## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockKnownClusters #############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a normal-block model with known clustering.
#' @export
NormalBlockKnownClusters <- R6::R6Class(
  classname = "NormalBlockKnownClusters",
  inherit   = NormalBlockBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockKnownClusters`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param C clustering matrix C_jk = 1 if species j belongs to cluster k
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockKnownClusters`] object
    initialize = function(data, C, sparsity = 0, control = NB_control()) {
      stopifnot("C must be a matrix" = is.matrix(C))
      stopifnot("There cannot be empty clusters" = min(colSums(C)) > 0)
      super$initialize(data, ncol(C), sparsity, control)
      private$C <- C
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      ddiag <- colMeans(reg_res$R^2)
      dm1   <- switch(private$res_covariance,
                      "diagonal"  = 1 / as.vector(ddiag),
                      "spherical" = rep(1/mean(ddiag), self$p))
      Sigmaq    <- private$heuristic_Sigmaq_from_Sigma(reg_res$Sigma)
      Omegaq    <- private$get_Omegaq(Sigmaq)
      list(B = reg_res$B, dm1 = dm1, Omegaq = Omegaq)
    },

    EM_initialize = function() {
      if (private$warm_started) {
        list(B = private$B, dm1 = private$dm1, Omegaq = private$Omegaq,
             gamma = private$gamma, mu = private$mu)
      } else c(private$get_heuristic_parameters(),  list(
        gamma = diag(1, self$q, self$q),
        mu    = matrix(0, self$n, self$q)
        )
      )
    },

    ## Runs the EM recursion via the Rcpp/Armadillo core (src/exports.cpp,
    ## NormalBlockKnownClusters_fit); see inst/normal_block_models.qmd §1/§2.
    EM_optimize = function(control) {
      init <- private$EM_initialize()
      res  <- NormalBlockKnownClusters_fit(
        Y = self$data$Y, X = self$data$X, C = private$C,
        B0 = init$B, dm1_0 = init$dm1, Omegaq0 = init$Omegaq,
        sparsity = self$sparsity, sparsity_weights = self$sparsity_weights,
        noise_covariance = private$res_covariance,
        niter = control$niter, threshold = control$threshold
      )
      private$niter <- res$niter
      list(B = res$B, dm1 = res$dm1, Omegaq = res$Omegaq,
           gamma = res$gamma, mu = res$mu, ll_list = res$objective)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function() list(gamma = private$gamma, mu = private$mu),
    #' @field entropy Entropy of the conditional distribution
    entropy    = function() {
      if (!private$approx){
        res <- .5 * self$n * self$q * log(2 * pi * exp(1)) +
          .5 * self$n * as.numeric(determinant(private$gamma)$modulus)
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        self$data$X %*% private$B
      } else {
        self$data$X %*% private$B + private$mu %*% t(private$C)
      }
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
      {paste(private$res_covariance, "normal-block model with fixed blocks")}
  )
)

