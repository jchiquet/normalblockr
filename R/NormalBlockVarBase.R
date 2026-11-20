## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for Normal-Block-Var Models
#'
#' R6 abstract class for the sparse Normal-Block models, where the clustering
#' structures the latent covariance.
#' @examples
#' # An internal abstract base class, never instantiated directly. See
#' # normal_block() for how concrete models (NormalBlockVarKnownClusters,
#' # NormalBlockVarUnknownClusters, and their zero-inflated variants) are
#' # actually created and fitted.
#' @keywords internal
NormalBlockVarBase <- R6::R6Class(
  classname = "NormalBlockVarBase",
  inherit   =  NormalBlockBase,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(

    #' @description Create a new [`NormalBlockVarBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' When `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`)
    #' is skipped entirely, since it would otherwise never be used downstream.
    #' @return A new [`NormalBlockVarBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control(), zero_inflation = FALSE) {
      ## family defaults (see NormalBlockMeanBase for the other ones)
      data <- data$for_family("var")
      if (is.null(control$clustering_init)) control$clustering_init <- "kmeans"
      if (is.null(control$noise_covariance)) control$noise_covariance <- "diagonal"
      stopifnot("noise_covariance = 'full' only applies to mean-block models, whose Sigma is the full p x p residual covariance" =
                  control$noise_covariance %in% c("diagonal", "spherical"))
      super$initialize(data, q, sparsity, control, zero_inflation)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, q, q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
      ## either diagonal or spherical residuals covariance
      private$res_covariance <- control$noise_covariance
    },

    #' @description Checks that the blocks gather positively co-varying
    #' variables, as the model requires (two variables of a block have a
    #' positive correlation). Compares, for each block, the empirical
    #' correlations of the residuals `Y - XB` with those implied by the model.
    #' @return a data frame with one row per block: `block`, `size`, the mean
    #' and minimum of the empirical within-block correlations (`mean_cor`,
    #' `min_cor`), the proportion of negative ones (`prop_negative`), and the
    #' mean correlation implied by the model (`model_cor`). `NA` for blocks
    #' with fewer than two variables.
    within_block_correlations = function() {
      R <- self$data$Y - self$data$X %*% private$B
      if (is.matrix(private$kappa)) R[self$data$zeros == 1] <- NA
      emp_cor   <- suppressWarnings(stats::cor(R, use = "pairwise.complete.obs"))
      Sigma_kk  <- diag(chol2inv(chol(private$Omega)))
      noise_var <- 1 / private$dm1
      cl <- self$clustering
      do.call(rbind, lapply(seq_len(self$q), function(k) {
        idx <- which(cl == k)
        if (length(idx) < 2)
          return(data.frame(block = k, size = length(idx), mean_cor = NA_real_, min_cor = NA_real_,
                            prop_negative = NA_real_, model_cor = NA_real_))
        pairs   <- upper.tri(diag(length(idx)))
        emp     <- emp_cor[idx, idx][pairs]
        sd_tot  <- sqrt(Sigma_kk[k] + noise_var[idx])
        implied <- (Sigma_kk[k] / tcrossprod(sd_tot))[pairs]
        data.frame(block = k, size = length(idx),
                   mean_cor = mean(emp, na.rm = TRUE), min_cor = min(emp, na.rm = TRUE),
                   prop_negative = mean(emp < 0, na.rm = TRUE), model_cor = mean(implied))
      }))
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after splitting cluster `index`
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @param rule how the members of the cluster are divided: `"kmeans"`
    #' runs a 2-means on their columns of residuals `Y - XB`; `"noise"` cuts
    #' them at the median of their noise variance. Default `NULL`, i.e.
    #' `"kmeans"` for this family.
    #' @return A new [`NormalBlockVarBase`] object
    split = function(index, in_place = FALSE, rule = NULL) {
      ## update private fields related to group parameters
      ## C, Omega, M, S, sparsity_weights

      ## split2 moves to the new cluster q + 1, split1 stays in `index`
      cl     <- self$clustering == index
      split2 <- private$split_members(index, rule)
      split1 <- cl & !split2

      new_C <- private$split_memberships(index, split1, split2)

      ## New block starts as a copy of its parent's M/S column; split1/split2
      ## index variables (length p) and must not be reused on M/S, which
      ## are indexed by *individuals* (length n).
      new_M <- cbind(private$M, private$M[, index])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- cbind(private$S, private$S[, index])
      } else {
        new_S <- c(private$S, private$S[index])
      }

      ## Sparsity weights
      if (self$q == 1) {
        new_weights <- matrix(c(0,1,1,0), 2, 2)
      } else {
        weights_cl <-  private$weights[index, setdiff(1:self$q, index)]
        weights_cl <-  c(weights_cl, mean(weights_cl))
        new_weights <- cbind(rbind(private$weights, weights_cl, deparse.level = 0),
                             c(weights_cl, 0))
      }

      ## Precision matrix re-derived from new_M/new_S.
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      private$reseeded(in_place, new_C, Omega = new_Omega, M = new_M, S = new_S, weights = new_weights)
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockVarBase`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      new_C <- private$merge_memberships(indices)

      ## Variational means
      new_M <- private$M[, -indices[2], drop = FALSE]
      new_M[, indices[1]] <- .5 * (private$M[, indices[1]] + private$M[, indices[2]])

      ## Variational variances
      if (is.matrix(private$S)) {
        new_S <- private$S[, -indices[2], drop = FALSE]
        new_S[, indices[1]] <- .5 * (private$S[, indices[1]] + private$S[, indices[2]])
      } else {
        new_S <- private$S[-indices[2]]
        new_S[indices[1]] <- .5 * (private$S[indices[1]] + private$S[indices[2]])
      }

      ## Sparsity weights
      new_weights <-  private$weights[-indices[2], -indices[2], drop = FALSE]

      ## Precision matrix re-derived from new_M/new_S.
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      private$reseeded(in_place, new_C, Omega = new_Omega, M = new_M, S = new_S, weights = new_weights)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## heldout_loglik() (NormalBlockBase): mean of the new rows, and their
    ## exact log-likelihoods, marginal over W: Y_i ~ N(B'X_i, C Sigma C' + D).
    ## Woodbury and the determinant lemma, with A = Omega + C'D^{-1}C (q x q):
    ##   log|C Sigma C' + D| = log|D| - log|Omega| + log|A|
    ##   r'(C Sigma C' + D)^{-1} r = r'D^{-1}r - u'A^{-1}u, u = C'D^{-1}r.
    ## `observed` (zero-inflated models): the marginal of each row's non-zero
    ## cells, i.e. the same formulas on the observed sub-rows of C, D and r.
    heldout_mean = function(X, C) X %*% private$B,
    heldout_rows = function(R, C, observed) {
      Omega <- private$Omega; dm1 <- private$dm1
      ld_Omega <- as.numeric(determinant(Omega)$modulus)
      if (is.null(observed)) {
        cA <- chol(Omega + crossprod(C * dm1, C))
        RD <- R * matrix(dm1, nrow(R), ncol(R), byrow = TRUE)
        W  <- backsolve(cA, t(RD %*% C), transpose = TRUE)
        quad <- rowSums(R * RD) - colSums(W^2)
        logdet <- -sum(log(dm1)) - ld_Omega + 2 * sum(log(diag(cA)))
        return(-0.5 * (ncol(R) * log(2 * pi) + logdet + quad))
      }
      vapply(seq_len(nrow(R)), function(i) {
        o <- observed[i, ]; d <- dm1[o]; r <- R[i, o]; Co <- C[o, , drop = FALSE]
        cA <- chol(Omega + crossprod(Co * d, Co))
        w  <- backsolve(cA, crossprod(Co, d * r), transpose = TRUE)
        logdet <- -sum(log(d)) - ld_Omega + 2 * sum(log(diag(cA)))
        -0.5 * (sum(o) * log(2 * pi) + logdet + sum(d * r^2) - sum(w^2))
      }, numeric(1))
    },
    ## heuristics tried by best_of_inits()
    default_inits     = c("kmeans", "ward2", "spectral"),
    ## split() rule used by refine() unless told otherwise (see split())
    default_split_rule = "kmeans",
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    family_name       = "normal-block-var",

    ## Variable-specific B and d are informed by the n observations, and so is
    ## the q x q precision matrix of the n latent vectors W_i (diagonal plus
    ## the non-zero off-diagonal terms). A spherical noise variance is shared
    ## by all n * p cells.
    family_penalty_groups = function() {
      noise <- if (private$res_covariance == "diagonal") {
        penalty_group("noise", self$p, self$n)
      } else {
        penalty_group("noise", 1, self$n * self$p)
      }
      rbind(penalty_group("B", self$p * self$d, self$n),
            penalty_group("Omega", self$q + self$n_edges, self$n),
            noise)
    },

    ## Cluster pairs for merge() are ranked by the signed correlation between
    ## the two blocks' latent signals, highest first: a merged block can only
    ## describe positively co-varying variables, so merging two nearly
    ## independent or negatively correlated blocks is rarely competitive.
    merge_score = function(pairs) {
      block_cor <- cov2cor(chol2inv(chol(private$Omega)))
      map_dbl(pairs, function(ij) block_cor[ij[1], ij[2]])
    },
    ## split(rule = "kmeans") clusters the columns of the residuals, which
    ## carry the latent signal W shared by the members of a block
    split_target = function() self$data$Y - self$data$X %*% private$B,
    ## split(rule = "noise") cuts at the median of the noise variances
    noise_variance = function() 1 / private$dm1,

    ## Closed-form M-step estimate of Omega from M/S (Sigma_hat = M'M/n +
    ## diag(S)); used by split()/merge() to reseed Omega from the new
    ## clustering's own M/S rather than the parent's stale Omega.
    omega_from_M_S = function(M, S, weights) {
      s_vec <- if (is.matrix(S)) colMeans(S) else S
      Sigma_hat <- crossprod(M) / self$n + diag(s_vec, ncol(M))
      ensure_pd(private$get_Omega(Sigma_hat, weights))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------

    ## Per-variable inverse variance from a residual matrix ("diagonal": one
    ## dm1 per variable; "spherical": one shared value). ddiag is floored
    ## away from 0 to avoid dm1 = Inf on a near-constant residual.
    dm1_from_residuals = function(R) {
      ddiag <- pmax(colMeans(R^2), .Machine$double.eps)
      switch(private$res_covariance,
             "diagonal"  = 1 / as.vector(ddiag),
             "spherical" = rep(1 / mean(ddiag), self$p))
    },

    ## Masked counterpart for zero-inflated models: the per-variable inverse
    ## variances of the masked least-squares fit, pooled over the non-zero
    ## cells when the noise is spherical
    zi_dm1 = function() {
      fit <- self$data$zi_ols_fit()
      switch(private$res_covariance,
             "diagonal"  = fit$dm1,
             "spherical" = rep(zi_pooled_precision(fit$dm1, self$data$nY), self$p))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    heuristic_optimize = function(control){
      parameters <- private$get_heuristic_parameters()
      c(parameters, list(ll_list = NA))
    },

    heuristic_Sigmaq_from_Sigma = function(Sigma){
      Sigma_q <- crossprod(private$C, Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
      ## NA when a cluster is empty (colSums(C) == 0 there -> division by 0)
      if (anyNA(Sigma_q)) {
        diag(Sigma_q)[is.na(diag(Sigma_q))] <- mean(diag(Sigma_q)[!is.na(diag(Sigma_q))])
        Sigma_q[is.na(Sigma_q)] <- 0
      }
      Sigma_q
    }

  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field B_original regression coefficients (d x p), converted back to
    #' Y's original units (undoing `NormalBlockData(scale = TRUE)`'s
    #' column-wise rescaling, if any). Use `model_par$B` instead for the
    #' coefficients on the internal fitting scale.
    B_original = function() private$rescale_to_original(private$B, power = 1),
    #' @field fitted Y values predicted by the model, in Y's original units:
    #' `X B`, plus the posterior mean of the latent block signals for (V)EM
    #' fits, with the observed zeros of zero-inflated models set to 0
    fitted = function() {
      res <- self$data$X %*% private$B
      if (!private$approx)
        res <- res + tcrossprod(if (private$known_clusters) private$mu else private$M, private$C)
      if (private$zi) res <- res * self$data$zeros_bar
      private$rescale_to_original(res)
    },
    #' @field model_par a list with the model parameters: B (d x p regression
    #' coefficients), B0 (zero-inflation coefficients, if any), dm1 (inverse
    #' noise variances of the variables) and Omega (q x q precision matrix of
    #' the latent block signals), on the scale of the fit. Use
    #' `$B_original`/`$dm1_original` for B and dm1 on Y's original scale.
    model_par = function() list(B = private$B, B0 = private$B0,
                                dm1 = private$dm1, Omega = private$Omega),
    #' @field sparsity_weights (weights associated to each pair of groups)
    sparsity_weights = function(value) {
      if (missing(value)) {
        private$weights
      } else {
        stopifnot("must be a q x q matrix" =
                    all(is.matrix(value), nrow(value) == ncol(value), ncol(value) == self$q))
        private$weights <- value
      }
    },
    #' @field dm1_original inverse residual variance per variable
    #' (1 / Var(Y_j)), converted back to Y's original units. Use
    #' `model_par$dm1` instead for the internal fitting scale. With
    #' `noise_covariance = "spherical"`, `model_par$dm1` is a single value
    #' repeated p times; once
    #' converted back per-variable, the p values returned here generally
    #' differ from one another whenever Y's columns were rescaled by
    #' different factors.
    dm1_original = function() private$rescale_to_original(private$dm1, power = -2)
  )
)
