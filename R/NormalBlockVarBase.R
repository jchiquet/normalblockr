## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for Variance-Block Models
#'
#' R6 abstract class for the sparse Normal-Block models, where the clustering
#' structures the latent covariance.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
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
      ## NB_control()'s clustering_init default is NULL, resolved here to
      ## this family's own default (see NormalBlockMeanBase for the other
      ## family's, and why they differ).
      if (is.null(control$clustering_init)) control$clustering_init <- "ward2"
      super$initialize(data, q, sparsity, control)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, q, q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
      ## variant (either diagonal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance

      ## kappa/B0/ZI_cond_mean are only used by the ZI subclasses; skip the
      ## p logistic regressions entirely for plain models.
      if (zero_inflation) {
        if(self$data$npY < self$n * self$p){
          B0_list <- lapply(1:self$data$p,
                            f <- function(j){
                              df <- data.frame("zeros" = data$zeros[,j], self$data$X0)
                              model <- glm(zeros ~ 0 + ., family=binomial(link = "logit"), data=df)
                              return(model$coefficients)})
          private$B0 <- t(sapply(B0_list, unlist))
        }else{
          private$B0 <- matrix(rep(-Inf, self$data$p * self$data$d0), nrow = self$data$d0)
        }

        private$kappa <- apply(self$data$X0 %*% private$B0, MARGIN = c(1,2), FUN = sigmoid)
        private$ZI_cond_mean <-
          sum(xlogy(data$zeros, private$kappa)) +
          sum(xlogy(data$zeros_bar, 1 - private$kappa))
      }
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update a [`NormalBlockBase`] object
    #'
    #' All possible parameters of the child classes
    #' @param B regression matrix
    #' @param dm1 diagonal vector of inverse variance matrix (variables level)
    #' @param C the matrix of groups memberships (posterior probabilities)
    #' @param Omega groups inverse variance matrix
    #' @param gamma  variance of posterior distribution of W
    #' @param mu mean for posterior distribution of W
    #' @param kappa vector of zero-inflation probabilities
    #' @param alpha vector of groups probabilities
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list  list of log-lik (elbo) values
    #' @param warm_started whether `EM_initialize()` should treat the model as
    #' already initialized (reuse B/Omega/dm1/C/alpha/M/S as they stand)
    #' rather than recomputing a fresh heuristic initialization -- set by
    #' [warm_start_from()] and by [split()]/[merge()].
    #' @param clustering_init name of a clustering heuristic to switch to,
    #' re-derived at the next `optimize()` call instead of reusing the
    #' current state (see `NB_control(clustering_init = )`). Used by
    #' [best_of_inits()].
    #' @return Update the current [`normal`] object
    update = function(B = NA,
                      dm1 = NA,
                      C = NA,
                      Omega = NA,
                      gamma = NA,
                      mu = NA,
                      kappa = NA,
                      alpha = NA,
                      M = NA,
                      S = NA,
                      ll_list = NA,
                      warm_started = NA,
                      clustering_init = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(C))       private$C       <- C
      if (!anyNA(Omega))  private$Omega  <- Omega
      if (!anyNA(gamma))   private$gamma   <- gamma
      if (!anyNA(kappa))   private$kappa   <- kappa
      if (!anyNA(mu))      private$mu      <- mu
      if (!anyNA(alpha))   private$alpha   <- alpha
      if (!anyNA(M))       private$M       <- M
      if (!anyNA(S))       private$S       <- S
      if (!anyNA(ll_list)) private$ll_list <- ll_list
      if (!anyNA(warm_started)) private$warm_started <- warm_started
      if (!anyNA(clustering_init)) {
        stopifnot("clustering_init must be a single heuristic name" =
                    is.character(clustering_init) && length(clustering_init) == 1)
        private$clustering_approx <- clustering_init
        private$C <- matrix(NA, self$p, self$q)
        private$warm_started <- FALSE
      }
    },

    #' @description Try several clustering-initialization heuristics and keep
    #' the best-ELBO converged fit (see `NB_control(clustering_init = )` and
    #' `inst/methods_initialization_and_refine.md` for the rationale). Every
    #' candidate is first screened with a short `trial_niter` run (same idea
    #' as `candidates_split()`/`candidates_merge()`), and only the
    #' `max_training` best-screened ones are fully retrained with `control`.
    #' @param inits vector of clustering-heuristic names to try
    #' @param trial_niter number of (V)EM iterations used to cheaply screen
    #' every candidate in `inits` before fully retraining the best few
    #' @param max_training how many of the screened candidates (best `loglik`
    #' after `trial_niter` iterations) get fully retrained with `control`
    #' @param control `optimize()` control list (`niter`/`threshold`) used
    #' for the final full retraining of the `max_training` best candidates
    #' @return a new, already-optimized [`NormalBlockVarBase`] object. Does not
    #' mutate the current object; reassign the result
    #' (`model <- model$best_of_inits()`).
    best_of_inits = function(inits = c("ward2", "kmeans", "spectral"),
                             trial_niter = 10, max_training = 2,
                             control = list(niter = 500, threshold = 1e-4)) {
      stopifnot(
        "best_of_inits() requires (V)EM inference (not the heuristic-only mode, see NB_control(heuristic = ))" =
          !private$approx,
        "best_of_inits() only applies when the initial clustering is inferred by a heuristic (see NB_control(clustering_init = ))" =
          !is.na(private$clustering_approx)
      )
      candidates <- map(inits, function(init) {
        cand <- self$clone()
        cand$update(clustering_init = init)
        cand$optimize(list(niter = trial_niter, threshold = 1e-4), warn = FALSE)
        cand
      })
      ibest <- order(map_dbl(candidates, "loglik"), decreasing = TRUE)[1:min(max_training, length(candidates))]
      best_candidates <- candidates[ibest]
      map(best_candidates, function(cand) cand$optimize(control, warn = FALSE))
      best_candidates[[which.max(map_dbl(best_candidates, "loglik"))]]
    },

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization process
    #' @param warn whether to warn when the (V)EM stops at the `niter` cap
    #' without reaching `threshold` (see `private$warn_if_not_converged()`).
    #' Set to `FALSE` for deliberately-truncated trial fits (cheap candidate
    #' scoring in `candidates_split()`/`candidates_merge()`, the sparsity-path
    #' warm-start probe in [NormalBlockCollectionSparsity]) where stopping at
    #' the cap is expected and not a sign of trouble.
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 500, threshold = 1e-4), warn = TRUE) {
      private$niter_max  <- control$niter
      private$threshold  <- control$threshold
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
      if (warn) private$warn_if_not_converged()
      invisible(self)
    },

    #' @description Seed this model's starting parameters from another,
    #' already-optimized model with the same q, instead of a fresh heuristic
    #' clustering. Used by [NormalBlockVarCollectionSparsity] to warm-start
    #' each penalty in a sparsity path from the previous one's solution.
    #' `B0`/`kappa` (zero-inflation) are left untouched: they depend only on
    #' the data, already set correctly and independently on every model.
    #' @param other a [NormalBlockVarBase] object, already optimized
    #' @return Update the current object in place with `other`'s parameters
    warm_start_from = function(other) {
      stopifnot("warm_start_from() requires both models to have the same q" = self$q == other$q)
      args <- other$model_par
      args$B0 <- NULL
      if (!is.null(other$var_par)) {
        args$C     <- other$var_par$tau
        args$M     <- other$var_par$M
        args$S     <- other$var_par$S
        args$alpha <- colMeans(other$var_par$tau)
      } else if (!is.null(other$posterior_par)) {
        args$gamma <- other$posterior_par$gamma
        args$mu    <- other$posterior_par$mu
      }
      do.call(self$update, args)
      private$warm_started <- TRUE
      invisible(self)
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after splitting cluster `cl`
    #' We split the cluster according to the species variances
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockVarBase`] object
    split = function(index, in_place = FALSE) {
      ## update private fields related to group parameters
      ## C, Omega, M, S, sparsity_weights

      ## indices of individuals split within the cluster
      cl  <- self$clustering == index
      var <- 1/private$dm1; var_median <- median(var[cl])
      split1 <- (var > var_median) & cl ;  split2 <- (var <= var_median) & cl

      ## Cluster split
      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## New block starts as a copy of its parent's M/S column; split1/split2
      ## index *variables* (length p) and must not be reused on M/S, which
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

      ## Precision matrix: re-derived from new_M/new_S (already consistent
      ## with the split), not hand-edited from the parent's Omega -- see
      ## omega_from_M_S().
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## Mark as already initialized so EM_initialize() reuses this state
      ## instead of a fresh heuristic clustering.
      if (in_place) {
        self$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                      alpha = colMeans(new_C), warm_started = TRUE)
        new_NB$sparsity_weights <- new_weights
        return(invisible(new_NB))
      }
    },

    #' @description generate and select a set of candidate models
    #' by splitting the clusters of the current model
    #' @param trial_niter number of EM iterations used to cheaply score each
    #' candidate before [SelectionNClusters] fully re-optimizes the best few
    #' (`train_best_candidates()`'s `max_training`) -- kept short on purpose.
    candidates_split = function(trial_niter = 5) {
      # do not split groups with less than 2 guys
      candidates <- map((1:self$q)[self$cluster_sizes > 1], self$split)
      # keep candidates with at least 2 guys per cluster and non empty split.
      # Compared against the number of currently *live* (non-empty) clusters
      # rather than self$q: the (V)EM can converge to an empty cluster (a
      # known general failure mode, e.g. plot_network()'s cluster_sizes fix),
      # in which case self$q overcounts and every split would otherwise look
      # invalid, silently stalling explore_forward() before n_clusters_range[2].
      clustering_sizes <- map(candidates, "clustering") %>% map(table)
      min_sizes  <- clustering_sizes %>% map_dbl(min)
      n_clusters <- clustering_sizes %>% map_dbl(length)
      n_live <- length(unique(self$clustering))
      candidates <- candidates[min_sizes > 1 & n_clusters == n_live + 1]

      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = trial_niter, threshold = 1e-4), warn = FALSE)
      candidates
    },

    #' @description generate and select a set of candidate models
    #' by merging the clusters of the current model
    #' @param max_candidates merge candidates are, unlike split's, quadratic
    #' in q (`choose(q, q-2)` pairs) -- beyond `max_candidates` pairs, only
    #' the most promising ones (largest `|Omega[i, j]|`, i.e. the most
    #' strongly related cluster pairs in the current fit) are actually built
    #' and trial-optimized, since merging two nearly independent blocks is
    #' rarely competitive anyway. Set to `Inf` to always try every pair.
    #' @param trial_niter see [candidates_split()]
    candidates_merge = function(max_candidates = 30, trial_niter = 2) {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      pairs <- combn(self$q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates) {
        score <- map_dbl(pairs, function(ij) abs(private$Omega[ij[1], ij[2]]))
        pairs <- pairs[order(score, decreasing = TRUE)[1:max_candidates]]
      }
      candidates <- map(pairs, self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = trial_niter, threshold = 1e-4), warn = FALSE)
      candidates
    },

    #' @description Create a clone of the current [`NormalBlockVarBase`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockVarBase`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      ## drop = FALSE: merging q = 2 down to q = 1 would otherwise silently
      ## drop these to a plain vector/scalar.
      new_C <- private$C[, -indices[2], drop = FALSE]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]

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

      ## Precision matrix: re-derived from new_M/new_S (already consistent
      ## with the merge), not hand-edited from the parent's Omega -- see
      ## omega_from_M_S().
      new_Omega <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## See split()'s comment: mark as warm-started so EM_initialize()
      ## reuses this state instead of a fresh heuristic Sigmaq/Omega.
      if (in_place) {
        self$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omega = new_Omega, M = new_M, S = new_S,
                      alpha = colMeans(new_C), warm_started = TRUE)
        new_NB$sparsity_weights <- new_weights
        return(new_NB)
      }
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    kappa             = NA, # vector of zero-inflation probabilities
    B0                = NA, # vector of zero-inflation regression matrix
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    res_covariance    = NA, # shape of the residuals covariance (diagonal or spherical)
    ZI_cond_mean      = NA, # conditional mean of the ZI component (fixed),

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated (V)EM inference --------------

    ## Closed-form M-step estimate of Omega from M/S (Sigma_hat = M'M/n +
    ## diag(S)); used by split()/merge() to reseed Omega from the new
    ## clustering's own M/S/C rather than the parent's stale Omega.
    omega_from_M_S = function(M, S, weights) {
      s_vec <- if (is.matrix(S)) colMeans(S) else S
      Sigma_hat <- crossprod(M) / self$n + diag(s_vec, ncol(M))
      Omega <- if (private$sparsity_ == 0) {
        chol2inv(chol(Sigma_hat))
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma_hat, rho = private$sparsity_ * weights)
        if (anyNA(glasso_out$wi)) chol2inv(chol(Sigma_hat)) else Matrix::symmpart(glasso_out$wi)
      }
      ensure_pd(Omega)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## Per-variable inverse variance from a residual matrix ("diagonal": one
    ## dm1 per variable; "spherical": one shared value). ddiag is floored
    ## away from 0 to avoid dm1 = Inf on a near-constant residual.
    dm1_from_residuals = function(R) {
      ddiag <- pmax(colMeans(R^2), .Machine$double.eps)
      switch(private$res_covariance,
             "diagonal"  = 1 / as.vector(ddiag),
             "spherical" = rep(1 / mean(ddiag), self$p))
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of ZI Diagonal Normal distribution
    ## (the weighted-least-squares fit itself lives in zi_weighted_fit(),
    ## R/utils.R, shared with the collection-level zi_residuals())
    zi_diag_normal_inference = function(){
      fit <- zi_weighted_fit(self$data)
      list(B = fit$B, dm1 = fit$dm1, kappa = private$kappa, R = fit$R)
    },

    heuristic_optimize = function(control){
      parameters <- private$get_heuristic_parameters()
      c(parameters, list(ll_list = NA))
    },

    heuristic_Sigmaq_from_Sigma = function(Sigma){
      Sigma_q <- (t(private$C) %*% Sigma %*% private$C) / outer(colSums(private$C), colSums(private$C))
      ## NA when a cluster is empty (colSums(C) == 0 there -> division by 0)
      if (anyNA(Sigma_q)) {
        diag(Sigma_q)[is.na(diag(Sigma_q))] <- mean(diag(Sigma_q)[!is.na(diag(Sigma_q))])
        Sigma_q[is.na(Sigma_q)] <- 0
      }
      Sigma_q
    },

    ## Registry of clustering heuristics used to turn the OLS/ZI residuals R
    ## (n x p) into an initial clustering of the p variables into self$q
    ## groups (a vector of length p with values in 1:q), selectable via
    ## NB_control(clustering_init = ...). See
    ## inst/methods_initialization_and_refine.md and
    ## inst/clustering_initialization_benchmark for the rationale and the
    ## empirical comparison behind the "ward2" default.
    clustering_methods = list(
      kmeans   = function(R, q) kmeans(t(R), q, nstart = 30, iter.max = 50)$cluster,
      ## ward2_tree() (R/utils.R) also backs sbm_clustering_path()'s own
      ## fallback -- same computation, shared rather than duplicated.
      ward2    = function(R, q) cutree(ward2_tree(R), q),
      sbm      = function(R, q) {
        options <- list(verbosity = 0, exploreMin = q, exploreMax = q, plot = FALSE, nbCores = 1)
        mySBM <- sbm::estimateSimpleSBM(cov(R), "gaussian", estimOptions = options)
        mySBM$setModel(q)
        mySBM$memberships
      },
      spectral = function(R, q) {
        U <- eigen(cov(R), symmetric = TRUE)$vectors[, seq_len(q), drop = FALSE]
        U <- U / pmax(sqrt(rowSums(U^2)), 1e-10)
        kmeans(U, q, nstart = 30, iter.max = 50)$cluster
      }
    ),

    heuristic_clustering = function(R) {
      clustering <- private$clustering_methods[[private$clustering_approx]](R, self$q)
      if (length(unique(clustering)) < self$q) {
        ## ward2's hclust()/cutree() always yields exactly q groups (modulo
        ## exact tied merge heights), making it a robust fallback whenever
        ## the chosen heuristic collapses to fewer than q clusters.
        clustering <- private$clustering_methods$ward2(R, self$q)
      }
      C <- as_indicator(clustering)
      if (min(colSums(C)) < 1) warning("Initialization failed to place elements in each cluster")
      C
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
    #' @field d0 number of zi variables (dimensions in X0)
    d0 = function() self$data$d0,
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalBlockData(scale = TRUE)`) -- use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, B0 = private$B0,
                                dm1 = private$dm1, Omega = private$Omega),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(self$p * self$d + self$q + self$n_edges + nb_param_D)
    },
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
    #' repeated p times (one shared variance on the fitting scale); once
    #' converted back per-variable, the p values returned here generally
    #' differ from one another whenever Y's columns were rescaled by
    #' different factors -- correctly so, since a single shared *scaled*
    #' variance does not correspond to a single shared variance in the
    #' original, heterogeneous-scale units.
    dm1_original = function() private$rescale_to_original(private$dm1, power = -2),
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function() private$res_covariance
  )
)
