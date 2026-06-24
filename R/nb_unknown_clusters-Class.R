## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NB_fixed_q ###################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a normal-block model with fixed number of clusters (but unknown clustering).
#' @export
NB_fixed_q <- R6::R6Class(
  classname = "NB_fixed_q",
  inherit = NB,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS --------------------------------------
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field fixed_tau whether tau should be fixed at clustering_init during optimization, useful for stability selection
    fixed_tau = NULL,

    #' @description Create a new [`NB_fixed_q`] object.
    #' @param data contains the matrix of responses (Y) and the design matrix (X).
    #' @param q required number of groups
    #' @param sparsity sparsity penalty to add on blocks precision matrix for sparsity
    #' @param control structured list for specific parameters
    #' @return A new [`NB_fixed_q`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
      self$fixed_tau <- control$fixed_tau
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(

    compute_loglik  = function(B, dm1, Omegaq, alpha, C, M, S) {
      log_det_Omegaq <- as.numeric(determinant(Omegaq, logarithm = TRUE)$modulus)

      J <- -.5 * self$n * self$p * log(2 * pi * exp(1)) + .5 * self$n * sum(log(dm1))
      J <- J + .5 * self$n * log_det_Omegaq + .5 * self$n * sum(log(S))
      J <- J + sum(C %*% log(alpha)) - sum(xlogx(C))

      if (private$sparsity_ > 0) {
        ## when not sparse, this terms equal -n q /2 by definition of Omegaq_hat and simplifies
        J <- J + self$n*self$q / 2 - .5 * sum(diag(Omegaq %*% (crossprod(M) + self$n * diag(S, self$q, self$q))))
        J <- J - private$sparsity_ * sum(abs(private$sparsity__weights * Omegaq))
      }
      J
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference -----------------------

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      if (anyNA(private$C)) # if no initial clustering provided
        private$C <- private$heuristic_clustering(reg_res$R)
      private$C <- check_one_boundary(check_zero_boundary(private$C))
      Sigmaq    <- private$heuristic_Sigmaq_from_Sigma(reg_res$Sigma)
      Omegaq    <- private$get_Omegaq(Sigmaq)
      ddiag <- colMeans(reg_res$R^2)
      dm1   <- switch(private$res_covariance,
                      "diagonal"  = 1 / as.vector(ddiag),
                      "spherical" = rep(1/mean(ddiag), self$p))
      list(B = reg_res$B, Omegaq = Omegaq, dm1 = dm1,
           C = private$C, alpha = colMeans(private$C))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated inference ------------------------

    EM_initialize = function() {
      c(private$get_heuristic_parameters(),  list(
            M = matrix(rep(0, self$n * self$q), nrow = self$n),
            S = rep(0.1, self$q)
          )
      )
    },

    EM_step = function(B, Omegaq, dm1, alpha, C, M, S) {

      ## Auxiliary variables
      R     <- self$data$Y - self$data$X %*% B
      Gamma <- solve(Omegaq + diag(colSums(dm1 * C), self$q, self$q))

      # E step
      M <- R %*% (dm1 * C) %*% Gamma
      S <- diag(Gamma)

      if (self$q > 1 & !self$fixed_tau) {
        eta <- dm1 * crossprod(R, M) - .5 * outer(dm1,  colSums(M^2) + self$n * S)
        eta <- eta + outer(rep(1, self$p), log(alpha)) - 1
        C   <- t(check_zero_boundary(check_one_boundary(apply(eta, 1, softmax))))
      }

      # M step
      MCT   <- tcrossprod(M, C)
      B     <- self$data$XtXm1 %*% crossprod(self$data$X, self$data$Y - MCT)
      ddiag <- colMeans(R^2 - 2 * R * MCT + tcrossprod(M^2 + outer(rep(1, self$n), S), C))
      dm1   <- switch(private$res_covariance,
                     "diagonal"  = 1 / as.vector(ddiag),
                     "spherical" = rep(1/mean(ddiag), self$p))
      alpha  <- colMeans(C)
      Omegaq <- private$get_Omegaq(crossprod(M)/self$n +  diag(S, self$q, self$q))

      list(B = B, Omegaq = Omegaq, dm1 = dm1, alpha = alpha, C = C, M = M, S = S)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omegaq (groups precision matrix))
    model_par  = function() {
      parameters       <- super$model_par
      parameters$alpha <- private$alpha
      parameters},
    #' @field nb_param number of parameters in the model
    nb_param = function() {as.integer(super$nb_param + self$q - 1)}, # adding alpha
    #' @field var_par a list with the matrices of the variational parameters: M (means), S (variances), tau (posterior group probabilities)
    var_par    = function() list(M = private$M,  S = private$S, tau = private$C),
    #' @field entropy Entropy of the conditional distribution
    entropy    = function() {
      if (!private$approx){
        res <- .5 * self$n * self$q * log(2 * pi * exp(1)) + .5 * self$n * sum(log(private$S))
        res <- res - sum(xlogx(private$C))
      } else {res <- NA}
      res
    },
    #' @field fitted Y values predicted by the model
    fitted = function(){
      if (private$approx) {
        res <- self$data$X %*% private$B
      } else {
        res <- self$data$X %*% private$B + tcrossprod(private$M, private$C)
      }
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(value) {
      paste(private$res_covariance, "normal-block model with", self$q, "unknown blocks")
    }
  )
)
