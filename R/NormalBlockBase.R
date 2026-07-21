## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a generic sparse Normal Block model
NormalBlockBase <- R6::R6Class(
  classname = "NormalBlockBase",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data  = NULL,

    #' @description Create a new [`NormalBlockBase`] object.
    #' @param data object of NormalBlockData class, with responses and design matrix
    #' @param q number of block/cluster
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @param zero_inflation whether the concrete subclass models zero-inflation;
    #' set by the ZI subclasses themselves, not meant to be set by the end user.
    #' When `FALSE`, the (costly) zero-inflation probability fit (`kappa`/`B0`)
    #' is skipped entirely, since it would otherwise never be used downstream.
    #' @return A new [`NormalBlockBase`] object
    initialize = function(data, q, sparsity = 0, control = NB_control(), zero_inflation = FALSE) {
      self$data <- data

      stopifnot("There cannot be more blocks than there are entities to cluster" = q <= ncol(self$data$Y))

      ## variant (either diagonal or spherical residuals covariance)
      private$res_covariance <- control$noise_covariance

      ## pointer to the chosen optimization function
      private$optimizer <- ifelse(control$heuristic,
                                  private$heuristic_optimize,
                                  private$EM_optimize)
      private$approx <- control$heuristic

      ## penalty mask
      private$sparsity_ <- sparsity
      weights <- matrix(1, q, q)
      diag(weights) <- 0
      if (!is.null(control$sparsity_weights)) {
        weights <- control$sparsity_weights
      }
      private$weights <- weights

      ## control$clustering_init is either the name of a clustering heuristic
      ## (a single string, deferred to heuristic_clustering(), looked up in
      ## private$clustering_methods) or an actual clustering to use directly
      ## (a vector of labels or a p x q indicator matrix).
      cl0 <- control$clustering_init
      if (is.character(cl0) && length(cl0) == 1) {
        private$clustering_approx <- cl0
        private$C <- matrix(NA, self$data$n, q)
      } else if (!is.null(cl0)) {
        if (!is.vector(cl0) & !is.matrix(cl0)) stop("Labels must be encoded in vector of labels or indicator matrix")
        if (is.vector(cl0)) {
          if (any(cl0 < 1 | cl0 > q))
            stop("Cluster labels must be between 1 and q")
          if (length(cl0) != self$p)
            stop("Cluster labels must match the number of Y's columns")
          if (length(unique(cl0)) != q)
            stop("The number of clusters in the initial clustering must be equal to q.")
          cl0 <- as_indicator(cl0)
        } else {
          if (nrow(cl0) != self$p)
            stop("Cluster-indicating matrix must have as many rows as Y has columns")
          if (ncol(cl0) != q)
            stop("Cluster-indicating matrix must have q columns")
          if (min(colSums(cl0)) < 1)
            stop("There cannot be empty clusters in the initial clustering matrix.")
        }
        private$C <- cl0
      } else {
        private$C <- matrix(NA, self$data$n, q)
      }

      ## Zero-inflation probabilities (kappa/B0) and the resulting fixed
      ## log-likelihood contribution (ZI_cond_mean) are only ever read by the
      ## ZI subclasses (see zi_diag_normal_inference(), and the ZI EM_optimize()/
      ## fitted/model_par methods) -- skip the p logistic regressions entirely
      ## for plain (non zero-inflated) models, where they would just be wasted
      ## work (private$kappa/B0/ZI_cond_mean stay at their NA default).
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
    #' @param Omegaq groups inverse variance matrix
    #' @param gamma  variance of posterior distribution of W
    #' @param mu mean for posterior distribution of W
    #' @param kappa vector of zero-inflation probabilities
    #' @param alpha vector of groups probabilities
    #' @param M variational mean for posterior distribution of W
    #' @param S variational diagonal of variances for posterior distribution of W
    #' @param ll_list  list of log-lik (elbo) values
    #' @param warm_started whether `EM_initialize()` should treat the model as
    #' already initialized (reuse B/Omegaq/dm1/C/alpha/M/S as they stand)
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
                      Omegaq = NA,
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
      if (!anyNA(Omegaq))  private$Omegaq  <- Omegaq
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
    #' @return a new, already-optimized [`NormalBlockBase`] object. Does not
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
    #' @param control a list for controlling the optimization proces
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
    #' already-optimized model with the same q, instead of the heuristic
    #' clustering-derived values set at construction time. Used by
    #' [NormalBlockCollectionSparsity] to warm-start each penalty in a sparsity
    #' path from the previous (adjacent) one's converged solution -- adjacent
    #' penalties along a sorted path usually have similar optima, so this
    #' typically needs far fewer EM iterations than starting cold each time
    #' (the same rationale as warm-starting in glmnet/glassoFast's own
    #' regularization paths). `B0`/`kappa` (zero-inflation) are deliberately
    #' left untouched: they depend only on the data, not on sparsity/blocks,
    #' so they are already set correctly and independently on every model.
    #' @param other a [NormalBlockBase] object, already optimized
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

    #' @description Create a clone of the current [`NormalBlockBase`] object after splitting cluster `cl`
    #' We split the cluster according to the species variances
    #' @param index index (integer) of the cluster to split
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockBase`] object
    split = function(index, in_place = FALSE) {
      ## update private fields related to group parameters
      ## C, Omegaq, M, S, sparsity_weights

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

      ## Variational means: the new block starts as a copy of its parent's
      ## column. split1/split2 (above) are computed over *variables*
      ## (length p, from self$clustering) and must not be reused to index
      ## M/S, which are indexed by *individuals* (length n): doing so used
      ## to silently scramble rows via recycling whenever p happened to
      ## divide n evenly, and errors outright otherwise -- e.g. the
      ## zero-inflated case, where S is also n x q rather than length q.
      ## EM_initialize()'s first E-step differentiates the duplicated
      ## columns immediately afterwards anyway, since it is driven by C/tau,
      ## which already differs between the two new blocks.
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
      ## with the split), not hand-edited from the parent's Omegaq -- see
      ## omega_from_M_S().
      new_Omegaq <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## Mark the result as already initialized (C/Omegaq/M/S/alpha freshly
      ## derived above, B/dm1 carried over unchanged from the parent -- valid
      ## since their dimension doesn't depend on q) so that EM_initialize()
      ## reuses this state instead of discarding it for a fresh heuristic
      ## Sigmaq/Omegaq estimate, which is both wasted work and a numerical
      ## stability risk (a candidate clustering from a split/merge is not
      ## constrained to be "nice" the way a clustering heuristic's output is).
      if (in_place) {
        self$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(invisible(self))
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S,
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
    #' the most promising ones (largest `|Omegaq[i, j]|`, i.e. the most
    #' strongly related cluster pairs in the current fit) are actually built
    #' and trial-optimized, since merging two nearly independent blocks is
    #' rarely competitive anyway. Set to `Inf` to always try every pair.
    #' @param trial_niter see [candidates_split()]
    candidates_merge = function(max_candidates = 30, trial_niter = 2) {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      pairs <- combn(self$q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates) {
        score <- map_dbl(pairs, function(ij) abs(private$Omegaq[ij[1], ij[2]]))
        pairs <- pairs[order(score, decreasing = TRUE)[1:max_candidates]]
      }
      candidates <- map(pairs, self$merge)
      for (i in seq_along(candidates))
        candidates[[i]]$optimize(list(niter = trial_niter, threshold = 1e-4), warn = FALSE)
      candidates
    },

    #' @description Create a clone of the current [`NormalBlockBase`] object after merging clusters `cl1` and `cl2`
    #' @param indices indices (couple of integer) of the clusters to merge
    #' @param in_place should the split applied to the object itself, or should a copy be sent?
    #' default FALSE (send a copy)
    #' @return A new [`NormalBlockBase`] object
    merge = function(indices, in_place=FALSE) {

      ## sorting by increasing group label
      indices <- sort(indices)

      ## Cluster merge. drop = FALSE throughout: merging from q = 2 down to
      ## q = 1 would otherwise silently drop these to plain vectors/a scalar
      ## (R's default behavior when only one row/column remains), breaking
      ## the matrix-indexed assignments right below and omega_from_M_S()'s
      ## ncol() call further down.
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
      ## with the merge), not hand-edited from the parent's Omegaq -- see
      ## omega_from_M_S().
      new_Omegaq <- private$omega_from_M_S(new_M, new_S, new_weights)

      ## See split()'s comment: mark as warm-started so EM_initialize()
      ## reuses this state instead of a fresh heuristic Sigmaq/Omegaq.
      if (in_place) {
        self$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S,
                    alpha = colMeans(new_C), warm_started = TRUE)
        self$sparsity_weights <- new_weights
        return(self)
      } else {
        new_NB <- self$clone()
        new_NB$update(C = new_C, Omegaq = new_Omegaq, M = new_M, S = new_S,
                      alpha = colMeans(new_C), warm_started = TRUE)
        new_NB$sparsity_weights <- new_weights
        return(new_NB)
      }
    },

    #' @description Predicts observations Y for new covariates X.
    #' @param new_X new set of covariates.
    #' @return A n*p prediction matrix for new observations
    predict = function(new_X){
      return(new_X %*% private$B)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract interaction network in the latent space
    #' @param type edge value in the network. Can be "support" (binary edges), "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species)
    #' @return a square matrix of size `self$q`
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omegaq != 0 & !diag(TRUE, ncol(private$Omegaq))),
        "precision"   = private$Omegaq,
        "partial_cor" = {
          tmp <- -private$Omegaq / tcrossprod(sqrt(diag(private$Omegaq))); diag(tmp) <- 1
          tmp
        }
      )
      ## Enforce sparse Matrix encoding to avoid downstream problems with igraph::graph_from_adjacency_matrix
      ## as it fails when given dsyMatrix objects
      Matrix::Matrix(net, sparse = TRUE)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param show_increment whether to add, below the objective trace, a second
    #' panel with the (log10) absolute increment between consecutive iterations
    #' and the convergence `threshold` used to stop optimize() (dashed line).
    #' That second panel is what actually tells convergence apart from merely
    #' running out of iterations: the objective trace alone tends to look flat
    #' well before the increment has actually crossed the threshold, especially
    #' as the number of blocks grows (see inst/CSDA_analyses).
    #' @description plots the evolution of the objective (log-likelihood or ELBO)
    #' across the (V)EM iterations of the last call to `optimize()`.
    #' @return a [`ggplot2::ggplot`] graph
    plot_loglik = function(show_increment = TRUE) {
      ll  <- private$ll_list
      obj <- self$objective
      if (length(obj) == 0 || all(is.na(obj))) {
        message("No objective trace to plot (heuristic inference does not compute a log-likelihood/ELBO).")
        return(invisible(NULL))
      }

      ## The objective (log-lik/ELBO) trace climbs fast then flattens out
      ## visually on a linear scale, well before it has actually converged
      ## (the same effect show_increment exists to catch). When every value
      ## is negative -- always true in practice in this model, since the
      ## leading -0.5*n*p*log(2*pi) term dominates -- plotting log10(-obj)
      ## instead spreads out the near-convergence iterations (where -obj is
      ## smallest) at the expense of the early, already-obvious big jumps,
      ## the same log-scale trick already used for the increment panel
      ## below. Left as the raw value in the (rare) edge case where that
      ## doesn't hold, e.g. a tiny/degenerate fit where the objective is
      ## not guaranteed negative.
      obj_is_neg <- all(obj < 0)
      obj_facet  <- if (obj_is_neg) "log10(-objective)" else "objective"
      obj_value  <- if (obj_is_neg) log10(-obj) else obj

      dplot <- tibble::tibble(iteration = seq_along(obj), value = obj_value, facet = obj_facet)
      last_increment <- abs(diff(ll))[length(obj)]
      converged <- !is.na(private$threshold) && last_increment < private$threshold
      subtitle  <- if (is.na(private$niter_max)) {
        sprintf("%d iterations", length(obj))
      } else if (converged) {
        sprintf("converged after %d iterations (threshold = %.1e)", length(obj), private$threshold)
      } else {
        sprintf("stopped at the %d-iteration cap -- last increment (%.3g) still above threshold (%.1e)",
                private$niter_max, last_increment, private$threshold)
      }

      if (show_increment) {
        dplot <- dplyr::bind_rows(
          dplot,
          tibble::tibble(iteration = seq_along(obj), value = log10(abs(diff(ll))), facet = "log10(|increment|)")
        )
      }
      dplot$facet <- factor(dplot$facet, levels = c(obj_facet, "log10(|increment|)"))

      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = iteration, y = value)) +
        ggplot2::geom_line() + ggplot2::geom_point(size = .8) +
        ggplot2::ggtitle(label = "(V)EM optimization", subtitle = subtitle) +
        ggplot2::xlab("iteration") + ggplot2::ylab(NULL) +
        ggplot2::theme_bw()

      if (show_increment) {
        p <- p + ggplot2::facet_wrap(~ facet, ncol = 1, scales = "free_y", strip.position = "left")
        if (!is.na(private$threshold))
          p <- p + ggplot2::geom_hline(
            data = data.frame(facet = factor("log10(|increment|)", levels = levels(dplot$facet)),
                               yintercept = log10(private$threshold)),
            ggplot2::aes(yintercept = yintercept),
            linetype = "dashed", colour = "red", alpha = .6
          )
      }
      p
    },

    #' @description plot the latent network.
    #' @param type edge value in the network. Either "precision" (coefficient of the precision matrix) or "partial_cor" (partial correlation between species).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. The default will use the column names ot the response matrix.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an optional igraph layout. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = igraph::layout_in_circle,
                            plot = TRUE) {
      if(anyNA(private$Omegaq)) stop("NA in the precision matrix")

      type   <- match.arg(type)
      output <- match.arg(output)

      net <- self$latent_network(type)

      if (output == "igraph") {
        G <-  igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        if (!is.null(node.labels)) {
          igraph::V(G)$label <- node.labels
        } else {
          igraph::V(G)$label <- unlist(lapply(1:ncol(net), f <- function(x) paste0("Cluster_", x)))
        }
        ## Nice nodes
        V.deg <- igraph::degree(G)/sum(igraph::degree(G))
        igraph::V(G)$label.cex <- V.deg / max(V.deg) + .5
        igraph::V(G)$size <- tabulate(self$clustering, nbins = self$q) * 100 / self$p
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        igraph::V(G)$frame.color <- "transparent" # NA triggers an igraph plot() warning; same (invisible) rendering
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 15*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- igraph::delete.vertices(G, which(igraph::degree(G) == 0))
        }
        if (plot) plot(G, layout = layout)
      }
      if (output == "corrplot") {
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          G <- net
          diag(net) <- 0
          corrplot::corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
        } else  {
          G <- net
        }
      }
      invisible(G)
    },

    #' @description plots the evolution of the objective during model optimization
    #' (see `plot_loglik()`)
    plot = function(){
      self$plot_loglik()
    },


    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## S3 methods ----------------------------
    #' @description User friendly print method
    #' @param model First line of the print output
    print = function(model = paste("A", self$who_am_I, ".\n")) {
      cat(model)
      cat("===========================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("===========================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $posterior_par / $var_par, $clustering \n")
      cat("    $loglik, $BIC, $ICL, $objective, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), coef(), sigma(), fitted(), predict() \n")
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    B                 = NA, # regression matrix
    dm1               = NA, # diagonal vector of inverse variance matrix (variables level)
    C                 = NA, # the matrix of posterior probabilities (tau) or group affectation
    Omegaq            = NA, # precision matrix for clusters
    kappa             = NA, # vector of zero-inflation probabilities
    B0                = NA, # vector of zero-inflation regression matrix
    alpha             = NA, # vector of groups probabilities
    gamma             = NA, # variance of  posterior distribution of W
    mu                = NA, # mean for posterior distribution of W
    M                 = NA, # variational mean for posterior distribution of W
    S                 = NA, # variational diagonal of variances for posterior distribution of W
    optimizer         = NA, # a link to the function that perform the optimization
    ll_list           = NA, # list of log-likelihoods or ELBOs
    sparsity_         = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    res_covariance    = NA, # shape of the residuals covariance (diagonal or spherical)
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # name of the clustering heuristic, key into clustering_methods
    ZI_cond_mean      = NA, # conditional mean of the ZI component (fixed)
    niter             = NA, # number of EM iterations required by the inference, if applicable
    niter_max         = NA, # niter cap passed to the last optimize() call (for plot_loglik()'s
                             # and warn_if_not_converged()'s "did it actually converge or just
                             # hit the cap?" diagnostic)
    threshold         = NA, # convergence threshold passed to the last optimize() call (idem)
    warm_started      = FALSE, # set by warm_start_from() and by split()/merge(): tells
                                # EM_initialize() to reuse the current B/dm1/Omegaq (and
                                # C/M/S/alpha or gamma/mu) instead of (re-)deriving them
                                # from the heuristic clustering.
                                # Deliberately NOT inferred from "are these fields non-NA",
                                # which would also be true for split()/merge() clones (those
                                # keep their own, different, already-tested initialization path).

    ## Warns when the (V)EM stopped because it hit the niter cap rather than
    ## because it actually converged below threshold -- silently otherwise
    ## (heuristic fits have no ll_list/niter to check). This is increasingly
    ## likely as the number of blocks q grows: more latent structure means a
    ## higher fraction of missing information, hence a slower EM convergence
    ## rate (see inst/CSDA_analyses) -- the fix for genuinely slow cases is to
    ## raise `niter` in NB_control(), not to loosen `threshold`.
    warn_if_not_converged = function() {
      if (is.na(private$niter) || private$niter < private$niter_max) return(invisible())
      last_increment <- abs(diff(private$ll_list))[private$niter]
      if (last_increment >= private$threshold)
        warning(sprintf(
          "%s: (V)EM stopped at the niter cap (%d) without reaching the convergence threshold (last increment = %.3g, threshold = %.1e). Consider raising `niter` in NB_control(), especially with many blocks -- see plot_loglik() to check.",
          self$who_am_I, private$niter_max, last_increment, private$threshold), call. = FALSE)
      invisible()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for integrated (V)EM inference --------------
    ## Each concrete subclass overrides EM_optimize() to call its
    ## Rcpp/Armadillo core (see src/exports.cpp and
    ## inst/normal_block_models.qmd); EM_initialize() (heuristic
    ## initialization) stays in R and supplies the starting values.
    EM_initialize = function() {},

    get_Omegaq = function(Sigma) {
      if (private$sparsity_ == 0) {
        Omega <- solve(Sigma)
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma, rho = private$sparsity_ * self$sparsity_weights)
        if (anyNA(glasso_out$wi)) {
          warning(
            "GLasso fails, the penalty is probably too small and the system badly conditionned \n reciprocal condition number =",
            rcond(Sigma), "\n We send back the original matrix and its inverse (unpenalized)."
          )
          Omega <- solve(Sigma)
        } else {
          Omega <- Matrix::symmpart(glasso_out$wi)
        }
      }
      Omega
    },

    ## Closed-form M-step estimate of Omegaq directly from M/S (Sigma_hat =
    ## M'M/n + diag(S)), the same formula the (V)EM's own M-step uses (see
    ## Sigma_hat_ in src/normal_block_unknown_clusters.h). Used by split()/
    ## merge() to seed the new model's Omegaq from its freshly built M/S/C
    ## (already consistent with the new clustering) instead of hand-editing
    ## the parent's Omegaq, which reflects the *old* clustering. `weights`
    ## is passed explicitly (rather than reading self$sparsity_weights, as
    ## get_Omegaq() does) because split()/merge() call this before the new,
    ## resized sparsity_weights have been assigned anywhere.
    omega_from_M_S = function(M, S, weights) {
      s_vec <- if (is.matrix(S)) colMeans(S) else S
      Sigma_hat <- crossprod(M) / self$n + diag(s_vec, ncol(M))
      Omega <- if (private$sparsity_ == 0) {
        solve(Sigma_hat)
      } else {
        glasso_out <- glassoFast::glassoFast(Sigma_hat, rho = private$sparsity_ * weights)
        if (anyNA(glasso_out$wi)) solve(Sigma_hat) else Matrix::symmpart(glasso_out$wi)
      }
      ensure_pd(Omega)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Methods for heuristic inference----------------------
    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ## Converts a per-variable quantity from the internal fitting scale
    ## (data$Y, rescaled column-wise by NormalBlockData(scale = TRUE)) back to
    ## Y's original units. Y_scaled = Y_original / Y_scale, so an additive
    ## quantity like B or a fitted value picks up one factor of Y_scale
    ## (power = 1: X_original = Y_scale * X_scaled), while dm1 = 1/variance
    ## picks up Y_scale^-2 (Var(Y_original) = Y_scale^2 * Var(Y_scaled), so
    ## dm1_original = dm1_scaled / Y_scale^2, power = -2). Used by the public
    ## B_original/dm1_original bindings and by fitted() in the 4 leaf classes
    ## -- model_par$B/model_par$dm1 stay on the internal scale unchanged, since
    ## warm_start_from() copies them directly into another model's private
    ## state and must not have them silently rescaled.
    rescale_to_original = function(M, power = 1) {
      factor <- self$data$Y_scale^power
      if (is.matrix(M)) M * matrix(factor, nrow(M), ncol(M), byrow = TRUE) else M * factor
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## MLE of MV Normal distribution
    multivariate_normal_inference = function(){
      B     <- self$data$XtXm1 %*% self$data$XtY
      R     <- ols_residuals(self$data)
      Sigma <- cov(R)
      list(B = B, R = R, Sigma = Sigma)
    },

    ## Per-variable inverse variance from a residual matrix, respecting
    ## res_covariance ("diagonal": one dm1 per variable; "spherical": a
    ## single shared value repeated p times). Shared by
    ## NormalBlockKnownClusters/NormalBlockUnknownClusters's
    ## get_heuristic_parameters() (previously duplicated verbatim in both).
    ## ddiag is floored away from exact zero: a variable with an
    ## (near-)constant residual otherwise produces dm1 = Inf (e.g. a
    ## variable observed at very few points relative to X's degrees of
    ## freedom -- see zi_weighted_fit() in R/utils.R for the same guard on
    ## the zero-inflation side).
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
      ## NA happens when a cluster ends up empty (colSums(C) == 0 for that
      ## column): every entry of outer(colSums(C), colSums(C)) involving it
      ## is then a division by 0 (see the "Initialization failed to place
      ## elements in each cluster" warning a few lines below this method's
      ## caller -- the heuristic clustering this feeds from does not reject
      ## empty clusters outright, unlike the known-clusters constructors).
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
    #' @field inference_method inference procedure used (heuristic or integrated with EM)
    inference_method = function() ifelse(private$approx, "heuristic", "integrated"),
    #' @field n number of samples
    n = function() self$data$n,
    #' @field p number of responses per sample
    p = function() self$data$p,
    #' @field d number of variables (dimensions in X)
    d = function() self$data$d,
    #' @field d0 number of zi variables (dimensions in X0)
    d0 = function() self$data$d0,
    #' @field q number of blocks
    q = function() as.integer(ncol(private$C)),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix Omegaq)
    n_edges  = function() sum(private$Omegaq[upper.tri(private$Omegaq, diag = FALSE)] != 0),
    #' @field model_par a list with the matrices of the model parameters: B (covariates), dm1 (species variance), Omegaq (groups precision matrix)). On the internal fitting scale (`self$data$Y`, possibly column-rescaled by `NormalBlockData(scale = TRUE)`) -- use `$B_original`/`$dm1_original` for the same quantities converted back to Y's original units.
    model_par = function() list(B = private$B, B0 = private$B0,
                                     dm1 = private$dm1, Omegaq = private$Omegaq),
    #' @field B_original regression coefficients (d x p), converted back to
    #' Y's original units (undoing `NormalBlockData(scale = TRUE)`'s
    #' column-wise rescaling, if any). Use `model_par$B` instead for the
    #' coefficients on the internal fitting scale.
    B_original = function() private$rescale_to_original(private$B, power = 1),
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
    #' @field nb_param number of parameters in the model
    nb_param = function() {
      nb_param_D <- ifelse(private$res_covariance == "diagonal", self$p, 1)
      as.integer(self$p * self$d + self$q + self$n_edges + nb_param_D)
    },
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1],
    #' @field loglik (or its variational lower bound)
    loglik = function() if (private$approx) NA else private$ll_list[[length(private$ll_list)]] + self$sparsity_term,
    #' @field deviance (or its variational lower bound)
    deviance = function() -2 * self$loglik,
    #' @field BIC (or its variational lower bound)
    #' @field entropy Entropy of the conditional distribution when applicable
    entropy    = function() 0,
    BIC = function() self$deviance + log(self$n) * self$nb_param,
    #' @field ICL variational lower bound of the ICL
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field EBIC variational lower bound of the EBIC
    EBIC   = function() self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(self$q), 0),
    #' @field criteria a vector with loglik, BIC and number of parameters
    criteria   = function() {
      data.frame(nb_param = self$nb_param, q = self$q, n_edges = self$n_edges, sparsity = self$sparsity,
                 loglik = self$loglik, deviance = self$deviance, BIC = self$BIC, ICL = self$ICL, EBIC = self$EBIC,
                 niter = private$niter )
    },
    #' @field sparsity (overall sparsity parameter)
    sparsity = function(value) {
      if (missing(value)) {
        private$sparsity_
      } else {
        stopifnot("must be a positive scale" = value >= 0)
        private$sparsity_ <- value
      }
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
    #' @field sparsity_term (sparsity_term term in log-likelihood due to sparsity)
    sparsity_term = function() self$sparsity * sum(abs(self$sparsity_weights * private$Omegaq)),
    #' @field get_res_covariance whether the residual covariance is diagonal or spherical
    get_res_covariance = function() private$res_covariance,
    #' @field memberships cluster memberships
    memberships = function() private$C,
    #' @field clustering given as the list of elements contained in each cluster
    clustering = function() {
      cl <- get_clusters(private$C)
      names(cl) <- colnames(self$data$Y)
      cl
    },
    #' @field cluster_sizes given as a vector of cluster sizes
    cluster_sizes = function() tabulate(self$clustering, nbins = self$q),
    #' @field elements_per_cluster given as the list of elements contained in each cluster
    elements_per_cluster = function() {
      if (is.null(names(self$clustering)))
        base::split(1:self$p, self$clustering)
      else
        base::split(names(self$clustering), self$clustering)
    }
  )
)
