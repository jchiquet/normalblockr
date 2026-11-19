## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanBase ##########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for Normal-Block-Mean Models
#'
#' R6 abstract class for the Normal-Block models where the clustering
#' structures the mean (mu_i = C B' X_i).
#' @examples
#' # An internal abstract base class, never instantiated directly -- use
#' # NormalBlockMeanKnownClusters / NormalBlockMeanUnknownClusters.
#' @keywords internal
NormalBlockMeanBase <- R6::R6Class(
  classname = "NormalBlockMeanBase",
  inherit = NormalBlockBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NormalBlockMeanBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_Mean_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' @return A new [`NormalBlockMeanBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control(),
                          zero_inflation = FALSE) {
      ## family default for scale = NULL (see NormalBlockData): unscaled
      data <- data$for_family("mean")
      ## family default for the initial clustering
      if (is.null(control$clustering_init)) control$clustering_init <- "kmeans"
      ## "diagonal" by default. Asking for sparsity means asking for off-diagonal
      ## structure, so it implies a full Sigma unless one was named explicitly.
      if (is.null(control$noise_covariance))
        control$noise_covariance <- if (isTRUE(sparsity > 0) && !zero_inflation) "full" else "diagonal"
      private$res_covariance <- control$noise_covariance
      ## A full Sigma ties the variables together within each row, and the
      ## zero-inflation mask leaves a different set of them observed in every
      ## row: each row would then need its own submatrix inverse (a
      ## missing-data EM rather than the reweighting the diagonal shape allows).
      stopifnot(
        "zero-inflated mean-block models only support noise_covariance = 'diagonal' or 'spherical'" =
          !(zero_inflation && control$noise_covariance == "full"),
        "sparsity > 0 is not available for zero-inflated mean-block models: it penalizes the off-diagonal terms of a full Sigma, which they do not carry" =
          !(zero_inflation && isTRUE(sparsity > 0))
      )
      stopifnot(
        "sparsity > 0 needs noise_covariance = 'full': a diagonal or spherical Sigma has no off-diagonal coefficient for the graphical lasso to penalize" =
          !(isTRUE(sparsity > 0) && control$noise_covariance != "full")
      )
      ## A full Sigma is p x p and estimated from n residuals, singular as
      ## soon as n <= p. The graphical lasso regularizes it.
      if (control$noise_covariance == "full" && !isTRUE(sparsity > 0) && data$n <= data$p)
        stop("mean-block models estimate a full p x p covariance from n observations, ",
             "so they need n > p (here n = ", data$n, ", p = ", data$p,
             "). Use sparsity > 0 to regularize it through the graphical lasso, ",
             "or reduce the number of variables.", call. = FALSE)
      super$initialize(data, q, sparsity, control, zero_inflation)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, self$data$p, self$data$p)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
    },

    #' @description Predicts observations Y for new covariates X, in Y's
    #' original units. The mean-block mean is mu_i = C B' X_i, so the
    #' cluster-level predictor has to be mapped back to the variables
    #' through C.
    #' @param new_X new set of covariates.
    #' @return A n*p prediction matrix for new observations
    predict = function(new_X) {
      private$rescale_to_original(tcrossprod(new_X %*% private$B, private$C))
    },

    #' @description Seed this model's starting parameters from another,
    #' already-optimized model with the same q, instead of the heuristic
    #' clustering-derived values set at construction time. Used by
    #' [split()]/[merge()].
    #' @param other a [NormalBlockMeanBase] object, already optimized
    #' @return Update the current object in place with `other`'s parameters
    warm_start_from = function(other) {
      stopifnot("warm_start_from() requires both models to have the same q" = self$q == other$q)
      args <- other$model_par
      ## var_par only exists on the unknown-clusters side; a known clustering
      ## is fixed and must not be carried over
      if (!is.null(other$var_par)) {
        args$C     <- other$var_par$tau
        args$alpha <- colMeans(other$var_par$tau)
      }
      do.call(self$update, args)
      private$warm_started <- TRUE
      invisible(self)
    },

    #' @description Create a clone of the current [`NormalBlockMeanBase`]
    #' object after splitting cluster `index`. Unlike the variance-block
    #' family, Omega and the sparsity weights are p x p here and do not
    #' depend on q, so they carry over unchanged; only C (tau) and B (one
    #' column per cluster) are affected.
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split be applied to the object itself, or
    #' should a copy be sent? Default FALSE (send a copy)
    #' @param rule how the members of the cluster are divided: `"kmeans"`
    #' runs a 2-means on their fitted mean trajectories `X B_ols`; `"noise"`
    #' cuts them at the median of their noise variance `1 / diag(Omega)`.
    #' Default `NULL`, i.e. `"noise"` for this family.
    #' @return A new [`NormalBlockMeanBase`] object
    split = function(index, in_place = FALSE, rule = NULL) {
      ## split2 moves to the new cluster q + 1, split1 stays in `index`
      cl     <- self$clustering == index
      split2 <- private$split_members(index, rule)
      split1 <- cl & !split2

      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## the next M-step differentiates the duplicated column via C/tau
      new_B <- cbind(private$B, private$B[, index])

      if (in_place) {
        self$update(C = new_C, B = new_B, alpha = colMeans(new_C), warm_started = TRUE)
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, B = new_B, alpha = colMeans(new_C), warm_started = TRUE)
        return(invisible(new_NB))
      }
    },

    #' @description Create a clone of the current [`NormalBlockMeanBase`]
    #' object after merging clusters `indices`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the merge be applied to the object itself, or
    #' should a copy be sent? Default FALSE (send a copy)
    #' @return A new [`NormalBlockMeanBase`] object
    merge = function(indices, in_place = FALSE) {
      indices <- sort(indices)

      new_C <- private$C[, -indices[2], drop = FALSE]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]

      new_B <- private$B[, -indices[2], drop = FALSE]
      new_B[, indices[1]] <- .5 * (private$B[, indices[1]] + private$B[, indices[2]])

      if (in_place) {
        self$update(C = new_C, B = new_B, alpha = colMeans(new_C), warm_started = TRUE)
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, B = new_B, alpha = colMeans(new_C), warm_started = TRUE)
        return(new_NB)
      }
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## SQUAREM by default only for a full Sigma: with a diagonal or spherical
    ## one the EM is exact and fast, and the acceleration only adds cost.
    use_acceleration = function(control) {
      if (is.null(control$acceleration)) private$res_covariance == "full"
      else isTRUE(control$acceleration)
    },
    ## heldout_loglik() (NormalBlockBase): mean C B'X_i of the new rows, and
    ## their Gaussian log-likelihoods with precision Omega. Zero-inflated
    ## mean-block models have a diagonal Omega, so each non-zero cell
    ## contributes independently.
    heldout_mean = function(X, C) tcrossprod(X %*% private$B, C),
    heldout_rows = function(R, C, observed) {
      Omega <- private$Omega
      if (is.null(observed)) {
        cO <- chol(Omega)
        quad <- rowSums(tcrossprod(R, cO)^2)
        return(-0.5 * (ncol(R) * log(2 * pi) - 2 * sum(log(diag(cO))) + quad))
      }
      w <- matrix(diag(Omega), nrow(R), ncol(R), byrow = TRUE)
      rowSums(ifelse(observed, 0.5 * log(w / (2 * pi)) - 0.5 * w * R^2, 0))
    },
    ## heuristics tried by best_of_inits()
    default_inits     = c("kmeans", "ward2", "spectral"),
    ## split() rule used by refine() unless told otherwise (see split())
    default_split_rule = "noise",
    res_covariance    = NA, # shape of Sigma: "full", "diagonal" or "spherical"

    ## Each cluster profile beta_k is shared by every cell of its block, hence
    ## charged log(n * p) like the block parameters of a latent block model;
    ## so is a spherical variance. Per-variable variances and a full
    ## precision matrix (diagonal plus non-zero off-diagonal terms) are
    ## informed by the n observations.
    penalty_groups = function() {
      cov <- switch(private$res_covariance,
                    "full"      = penalty_group("Omega", self$p + self$n_edges, self$n),
                    "diagonal"  = penalty_group("noise", self$p, self$n),
                    "spherical" = penalty_group("noise", 1, self$n * self$p))
      rbind(penalty_group("B", self$q * self$d, self$n * self$p), cov)
    },
    ## the network relates the p variables (full Sigma only)
    network_dim = function() self$p,
    ## split(rule = "kmeans") clusters each variable's own fitted mean
    ## trajectory, as the initial clustering does
    split_target = function() self$data$X %*% self$data$ols_fit()$B,
    ## split(rule = "noise") cuts at the median of the noise variances
    noise_variance = function() 1 / diag(private$Omega),
    Psi               = NA,
    Phi               = NA,
    Lambda            = NA,

    ## Precision matrix for the requested shape of Sigma. Mirrors
    ## omega_from_residuals() in src/normal_block_mean_base.h
    omega_from_sigma = function(Sigma) {
      switch(private$res_covariance,
             "full"      = private$get_Omega(Sigma),
             "diagonal"  = diag(1 / diag(Sigma), nrow(Sigma)),
             "spherical" = diag(1 / mean(diag(Sigma)), nrow(Sigma)))
    },

    ## Warm-started counterpart of omega_from_sigma(), used only by the R
    ## reference recursion so that it mirrors the C++ core's M-step, which
    ## resumes each graphical lasso from the previous one (nb_omega::estimate,
    ## src/omega_estimation.h). Kept separate so omega_from_sigma() stays
    ## stateless for the heuristic and split()/merge() paths, where a carried
    ## -over state would be stale.
    omega_from_sigma_warm = function(Sigma, warm = NULL) {
      if (private$res_covariance != "full" || !isTRUE(private$sparsity_ > 0))
        return(list(Omega = private$omega_from_sigma(Sigma), warm = NULL))
      out <- graphical_lasso_fit(Sigma, private$sparsity_ * self$sparsity_weights,
                                 thr = NB_GLASSO_THRESHOLD,
                                 w_init = warm$w, wi_init = warm$wi)
      if (anyNA(out$wi)) return(list(Omega = chol2inv(chol(Sigma)), warm = NULL))
      list(Omega = Matrix::symmpart(out$wi), warm = list(w = out$w, wi = out$wi))
    },

    ## Omega is p x p here, so it carries no cluster-pair information: the
    ## closest mean profiles (columns of B) are the promising merges instead.
    merge_score = function(pairs) {
      -map_dbl(pairs, function(ij) sum((private$B[, ij[1]] - private$B[, ij[2]])^2))
    },

    ## Masked counterpart of multivariate_normal_inference(), used by the ZI
    ## subclasses to initialize.
    zi_mean_inference = function() {
      fit <- self$data$zi_ols_fit()
      list(B = fit$B, Omega = private$omega_from_sigma(diag(1 / fit$dm1, self$p)), R = fit$R)
    },

    heuristic_cluster_B_from_variable_B = function(B_variable, C){
      B <- B_variable %*% C / rep(colSums(C), each = nrow(B_variable))
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalBlockData(scale = TRUE)`): use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, Omega = private$Omega),
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a p x p matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$p))
        private$weights <- value
      }
    }
)
)
