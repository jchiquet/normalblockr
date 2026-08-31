## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanBase ##########################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for Mean-Block Models
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
    #' @param data object of NormalMeanBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_Mean_control
    #' @return A new [`NormalBlockMeanBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control()) {
      ## NB_control()'s clustering_init default is NULL, resolved here to
      ## "kmeans": clustering is done on each variable's fitted mean
      ## trajectory X %*% B_j (get_heuristic_parameters(), in the concrete
      ## subclasses); benchmarked against "ward2" on simulated data, kmeans
      ## consistently lands in a better basin (see best_of_inits() and the
      ## commit history for the numbers), so it is the more sensible default
      ## here even though "ward2" remains the variance-block family's.
      if (is.null(control$clustering_init)) control$clustering_init <- "kmeans"
      super$initialize(data, q, sparsity, control)
      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, self$data$p, self$data$p)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights
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
      args$C     <- other$var_par$tau
      args$alpha <- colMeans(other$var_par$tau)
      do.call(self$update, args)
      private$warm_started <- TRUE
      invisible(self)
    },

    #' @description Create a clone of the current [`NormalBlockMeanBase`]
    #' object after splitting cluster `index`. Unlike the variance-block
    #' family, Omega and the sparsity weights are p x p here and do not
    #' depend on q, so they carry over unchanged; only C (tau) and B (one
    #' column per cluster) are affected. Variables are split by their
    #' current noise variance (1 / diag(Omega)) around its within-cluster
    #' median -- the same criterion [NormalBlockVarBase]'s `split()` uses
    #' via `dm1`, since `diag(Omega)` plays the same per-variable-precision
    #' role here.
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split be applied to the object itself, or
    #' should a copy be sent? Default FALSE (send a copy)
    #' @return A new [`NormalBlockMeanBase`] object
    split = function(index, in_place = FALSE) {
      cl  <- self$clustering == index
      var <- 1 / diag(private$Omega); var_median <- median(var[cl])
      split1 <- (var > var_median) & cl; split2 <- (var <= var_median) & cl

      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C <- new_C / rowSums(new_C)

      ## The new cluster starts as a copy of its parent's B column; the
      ## next M-step differentiates the two since it is driven by C/tau,
      ## which already differs between them (same rationale as
      ## NormalBlockVarBase's split() duplicating M's column).
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

    #' @description generate and select a set of candidate models by
    #' splitting the clusters of the current model
    #' @param trial_niter number of (V)EM iterations used to cheaply score
    #' each candidate before fully re-optimizing the best few
    candidates_split = function(trial_niter = 5) {
      candidates <- map((1:self$q)[self$cluster_sizes > 1], self$split)
      ## keep candidates with at least 2 variables per cluster and a
      ## genuine split (see NormalBlockVarBase's candidates_split() for why
      ## this is compared against the number of currently *live* clusters
      ## rather than self$q)
      clustering_sizes <- map(candidates, "clustering") %>% map(table)
      min_sizes  <- clustering_sizes %>% map_dbl(min)
      n_clusters <- clustering_sizes %>% map_dbl(length)
      n_live <- length(unique(self$clustering))
      candidates <- candidates[min_sizes > 1 & n_clusters == n_live + 1]

      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = trial_niter, threshold = 1e-4, fixed_point_niter = 5), warn = FALSE)
      candidates
    },

    #' @description generate and select a set of candidate models by
    #' merging the clusters of the current model
    #' @param max_candidates merge candidates are quadratic in q -- beyond
    #' `max_candidates` pairs, only the most promising ones (smallest
    #' distance between the two clusters' B columns, i.e. the most similar
    #' mean profiles) are actually built and trial-optimized. Set to `Inf`
    #' to always try every pair.
    #' @param trial_niter see [candidates_split()]
    candidates_merge = function(max_candidates = 30, trial_niter = 2) {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      pairs <- combn(self$q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates) {
        dist <- map_dbl(pairs, function(ij) sum((private$B[, ij[1]] - private$B[, ij[2]])^2))
        pairs <- pairs[order(dist)[1:max_candidates]]
      }
      candidates <- map(pairs, self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = trial_niter, threshold = 1e-4, fixed_point_niter = 5), warn = FALSE)
      candidates
    },

    #' @description Create a clone of the current [`NormalBlockMeanBase`]
    #' object after merging clusters `indices`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the merge be applied to the object itself, or
    #' should a copy be sent? Default FALSE (send a copy)
    #' @return A new [`NormalBlockMeanBase`] object
    merge = function(indices, in_place = FALSE) {
      indices <- sort(indices)

      ## drop = FALSE: merging q = 2 down to q = 1 would otherwise silently
      ## drop these to a plain vector/scalar.
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
    },

    #' @description Try several clustering-initialization heuristics and
    #' keep the best-ELBO converged fit (see `NB_control(clustering_init = )`).
    #' Every candidate is first screened with a short `trial_niter` run, and
    #' only the `max_training` best-screened ones are fully retrained with
    #' `control`.
    #' @param inits vector of clustering-heuristic names to try
    #' @param trial_niter number of (V)EM iterations used to cheaply screen
    #' every candidate in `inits` before fully retraining the best few
    #' @param max_training how many of the screened candidates (best
    #' `loglik` after `trial_niter` iterations) get fully retrained with
    #' `control`
    #' @param control `optimize()` control list used for the final full
    #' retraining of the `max_training` best candidates
    #' @return a new, already-optimized [`NormalBlockMeanBase`] object. Does
    #' not mutate the current object; reassign the result
    #' (`model <- model$best_of_inits()`).
    best_of_inits = function(inits = c("kmeans", "ward2", "spectral"),
                             trial_niter = 10, max_training = 2,
                             control = list(niter = 500, threshold = 1e-4, fixed_point_niter = 5)) {
      stopifnot(
        "best_of_inits() requires (V)EM inference (not the heuristic-only mode, see NB_control(heuristic = ))" =
          !private$approx,
        "best_of_inits() only applies when the initial clustering is inferred by a heuristic (see NB_control(clustering_init = ))" =
          !is.na(private$clustering_approx)
      )
      candidates <- map(inits, function(init) {
        cand <- self$clone()
        cand$update(clustering_init = init)
        cand$optimize(list(niter = trial_niter, threshold = 1e-4, fixed_point_niter = 5), warn = FALSE)
        cand
      })
      ibest <- order(map_dbl(candidates, "loglik"), decreasing = TRUE)[1:min(max_training, length(candidates))]
      best_candidates <- candidates[ibest]
      map(best_candidates, function(cand) cand$optimize(control, warn = FALSE))
      best_candidates[[which.max(map_dbl(best_candidates, "loglik"))]]
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    Psi               = NA,
    Phi               = NA,
    Lambda            = NA,

    heuristic_cluster_B_from_variable_B = function(B_variable, C){
      B <- B_variable %*% C / rep(colSums(C), each = nrow(B_variable))
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omega (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalMeanBlockData(scale = TRUE)`) -- use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, Omega = private$Omega),
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      as.integer(self$q * self$d + self$p + self$n_edges)
    },
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
