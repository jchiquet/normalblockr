## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanUnknownClusters ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a normal-block model with known clustering.
#' @export
NormalBlockMeanUnknownClusters <- R6::R6Class(
  classname = "NormalBlockMeanUnknownClusters",
  inherit   = NormalBlockMeanBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanUnknownClusters`] object.
    #' @param data object of NormalMeanBlockData class, with responses and design matrix
    #' @param q number of clusters
    #' @param sparsity to apply on variance matrix when calling GLASSO
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockMeanUnknownClusters`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      super$initialize(data, q, sparsity, control)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    optim_initialize = function() {
      # warm restart to integrate later
      return(private$get_heuristic_parameters())
    },

    get_heuristic_parameters = function(){
      reg_res   <- private$multivariate_normal_inference()
      Omega     <- private$get_Omega(reg_res$Sigma)

      kmeans_init <- kmeans(t(reg_res$R), centers = self$q, nstart = 20)
      tau         <- as_indicator(kmeans_init$cluster)
      tau         <- check_one_boundary(check_zero_boundary(tau))
      tau         <- tau / rowSums(tau)
      B           <- private$heuristic_cluster_B_from_variable_B(reg_res$B,
                                                                 tau)
      list(B = B, Omega = Omega, tau = tau)
    },


    B_estimator = function(Omega = private$Omega,
                           tau = private$C,
                           Psi = private$Psi){
      return(self$data$XtXm1 %*% self$data$XtY %*% Omega %*% tau %*% solve(Psi))
    },

    Lambda_estimator = function(B = private$B,
                                tau = private$C) {

    M      <- crossprod(B, self$data$XtX) %*% B
    diag_M <- diag(M)

    diag_val1 <- as.vector(tau %*% diag_M)
    term1 <- diag(diag_val1, nrow = self$data$p, ncol = self$data$p)

    diag_val2 <- diag(tau %*% M %*% t(tau))
    term2 <- diag(diag_val2, nrow = self$data$p, ncol = self$data$p)

    return((term1 - term2) / self$data$n)
    },

    Sigma_estimator = function(Omega = private$Omega, B = private$B,
                               tau = private$C, Lambda = private$Lambda){

      R_bar <- self$data$Y - self$data$X %*% B %*% t(tau)

      Sigma <- (crossprod(R_bar) / self$data$n) + Lambda
      return(Sigma)
    },

    Psi_estimator = function(Omega = private$Omega, tau = private$C) {
      diag_Omega <- diag(Omega)
      term1 <- crossprod(tau, Omega) %*% tau

      diag_product <- as.vector(crossprod(tau, diag_Omega))
      term2 <- diag(diag_product, nrow = self$q, ncol = self$q)

      term3 <- crossprod(tau, diag(diag_Omega, self$data$p)) %*% tau

      Psi <- term1 + term2 - term3

      Psi <- Psi + 1e-8 * diag(self$q)

      return(Psi)
    },

    Phi_estimator = function(Omega = private$Omega, tau = private$C,
                             Psi = private$Psi) {
      tau_Omega_tau <- crossprod(tau, Omega) %*% tau
      Phi <- Psi - tau_Omega_tau
      return(Phi)
    },

    alpha_estimator = function(tau = private$C) {
      return(colMeans(tau))
    },

    tau_estimator = function(Omega = private$Omega,
                             B     = private$B,
                             alpha = private$alpha,
                             tau = private$C) {
      M      <- crossprod(B, self$data$XtX) %*% B
      diag_M <- diag(M)

      YtXB <- crossprod(self$data$Y, self$data$X %*% B)

      term1 <- - matrix(1, nrow = self$data$p, ncol = self$q)
      term2 <- Omega %*% (YtXB - tau %*% M)
      term2 <- pmin(pmax(term2, -50), 50)
      term3 <- -0.5 * matrix(diag(Omega), nrow = self$data$p, ncol = 1) %*% matrix(1, nrow = 1, ncol = self$q) %*% diag(diag_M)
      term4 <- diag(diag(Omega)) %*% tau %*% M

      exposant <- term1 + term2 + term3 + term4

      log_alpha <- log(pmax(alpha, 1e-300))
      log_numerator <- matrix(log_alpha, nrow = self$data$p,
                               ncol = self$q, byrow = TRUE) + exposant
      log_numerator <- log_numerator - apply(log_numerator, 1, max)
      tau <- exp(log_numerator)
      tau <- tau / rowSums(tau)
      tau <- check_one_boundary(check_zero_boundary(tau))
      tau <- tau / rowSums(tau)
      return(tau)
    },


    compute_loglik = function(B = private$B, Omega = private$Omega,
                              tau = private$C, alpha = private$alpha,
                              Phi = private$Phi){

      R_bar <- self$data$Y - tcrossprod(self$data$X %*% B, tau)
      M     <- crossprod(B, self$data$XtX) %*% B
      l <- - (self$n * self$p / 2) * log(2 * pi) +
           (self$n / 2) * as.numeric(determinant(Omega, logarithm = TRUE)$modulus) +
           sum(tau * outer(rep(1, self$p), log(alpha))) - sum(tau * log(tau)) -
           0.5 * (sum(R_bar * (R_bar %*% Omega)) + sum(diag(Phi %*% M)))

      return(as.numeric(l))
    },

    EM_optimize = function(control){
      init_params <- private$optim_initialize()
      B           <- init_params$B
      Omega       <- init_params$Omega
      tau         <- init_params$tau
      alpha       <- private$alpha_estimator(tau)
      Psi         <- private$Psi_estimator(Omega, tau)
      Phi         <- private$Phi_estimator(Omega, tau, Psi)
      Lambda      <- private$Lambda_estimator(B, tau)
      ll_prev     <- private$compute_loglik(B, Omega, tau, alpha, Phi)
      ll_list     <- c(ll_prev)

      for(i in 1:control$niter){
        # M-step
        Psi    <- private$Psi_estimator(Omega, tau)
        B      <- private$B_estimator(Omega, tau, Psi)
        Lambda <- private$Lambda_estimator(B, tau)
        Sigma  <- private$Sigma_estimator(Omega, B, tau, Lambda)
        Omega  <- private$get_Omega(Sigma)
        alpha  <- private$alpha_estimator(tau)

        # VE-step
        for (k in 1:control$fixed_point_niter) {
          tau <- private$tau_estimator(Omega, B, alpha, tau)
        }


        # elbo updating
        ll_current <- private$compute_loglik(B, Omega, tau, alpha, Phi)
        ll_list     <- c(ll_list, ll_current)
        if(abs(ll_current - ll_prev) < control$threshold){
          break
        }else{ll_prev <- ll_current}
      }
      private$niter <- i
      list(B = B, Omega = Omega, C = tau, alpha = alpha, Psi = Psi, Phi = Phi,
           Lambda = Lambda, ll_list = ll_list)
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field fitted Y values predicted by the model, in Y's original units
    fitted = function(){
      res <- if (private$approx) {
        self$data$X %*% private$B
      } else {
        self$data$X %*% private$B + private$mu %*% t(private$C)
      }
      private$rescale_to_original(res)
    },
    #' @field var_par a list with the variational parameter: tau (posterior group probabilities)
    var_par    = function() list(tau = private$C),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function()
    {paste("normal-block-mean model with unknown blocks")}
  )
)

