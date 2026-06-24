## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_blocks ##############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a normal-block model with known clustering.
#' @export
NB_fixed_blocks <- R6::R6Class(
  classname = "NB_fixed_blocks",
  inherit   = NB,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NB_fixed_blocks`] object.
    #' @param data object of NB_data class, with responses and design matrix
    #' @param C clustering matrix C_jk = 1 if species j belongs to cluster k
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NB_fixed_blocks`] object
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

    compute_loglik  = function(B, Omegaq, dm1 = NA, gamma = NA, mu = NA) {
      log_det_Omegaq <- as.numeric(determinant(Omegaq, logarithm = TRUE)$modulus)
      log_det_Gamma  <- as.numeric(determinant(gamma , logarithm = TRUE)$modulus)

      J <- -.5 * (self$n * self$p * log(2 * pi * exp(1)) - self$n * sum(log(dm1)))
      J <- J + .5 * self$n * (log_det_Omegaq + log_det_Gamma)
      if (self$sparsity > 0) {
        ## when not sparse, this terms equal -n q /2 by definition of Omegaq_hat
        J <- J + self$n*self$q / 2 - .5 * sum(diag(Omegaq %*% (self$n * gamma + crossprod(mu))))
        J <- J - self$sparsity * sum(abs(self$sparsity_weights * Omegaq))
      }
      J
    },

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
      c(private$get_heuristic_parameters(),  list(
        gamma = diag(1, self$q, self$q),
        mu    = matrix(0, self$n, self$q)
        )
      )
    },

    EM_step = function(B, dm1, Omegaq, gamma, mu) {
      ## E step
      gamma <- solve(Omegaq + diag(colSums(dm1 * private$C), self$q, self$q))
      mu    <- (self$data$Y - self$data$X %*% B) %*% (dm1 * private$C) %*% gamma

      ## M step
      YmmuCT <- self$data$Y - mu %*% t(private$C)
      B      <- self$data$XtXm1 %*% crossprod(self$data$X, YmmuCT)
      ddiag  <- colMeans((YmmuCT - self$data$X %*% B)^2) + private$C %*% diag(gamma)
      dm1  <- switch(private$res_covariance,
        "diagonal"  = 1 / as.vector(ddiag),
        "spherical" = rep(1/mean(ddiag), self$p))
      Omegaq <- private$get_Omegaq(crossprod(mu)/self$n + gamma)
      list(B = B, Omegaq = Omegaq, dm1 = dm1, gamma = gamma, mu = mu)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field posterior_par a list with the parameters of posterior distribution W | Y
    posterior_par = function(value) list(gamma = private$gamma, mu = private$mu),
    #' @field entropy Entropy of the conditional distribution
    entropy    = function(value) {
      if (!private$approx){
        res <- .5 * self$n * self$q * log(2 * pi * exp(1)) +
          .5 * self$n * as.numeric(determinant(private$gamma)$modulus)
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model
    fitted = function(value){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + private$mu %*% t(private$C)
      }
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function(value)
      {paste(private$res_covariance, "normal-block model with fixed blocks")}
  )
)

