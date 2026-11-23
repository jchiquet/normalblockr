## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockBase ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Root Base Class for Normal-Block Models
#'
#' R6 abstract class shared by the variance-block ([NormalBlockVarBase]) and
#' mean-block ([NormalBlockMeanBase]) model families.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how concrete models are created and fitted.
#' @keywords internal
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
    initialize = function(data, q, sparsity = 0, control, zero_inflation = FALSE) {
      self$data <- data
      stopifnot("There cannot be more blocks than there are entities to cluster" = q <= ncol(self$data$Y))

      ## pointer to the chosen optimization function
      private$optimizer <- if (control$heuristic) private$heuristic_optimize else private$EM_optimize
      stopifnot("this model family does not implement NB_control(heuristic = TRUE)" =
                  is.function(private$optimizer))
      private$approx <- control$heuristic
      private$zi <- zero_inflation

      ## control$clustering_init is either the name of a clustering heuristic
      ## or an actual clustering to use directly
      cl0 <- control$clustering_init
      if (is.character(cl0) && length(cl0) == 1) {
        private$clustering_approx <- cl0
        private$C <- matrix(NA, self$p, q)
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
          ## an indicator matrix needs a member per cluster; a matrix of
          ## posterior probabilities (e.g. the tau a stability selection
          ## starts from) only needs some mass in each cluster
          hard <- all(cl0 == 0 | cl0 == 1)
          if (min(colSums(cl0)) < (if (hard) 1 else .Machine$double.eps))
            stop("There cannot be empty clusters in the initial clustering matrix.")
        }
        private$C <- cl0
      } else {
        private$C <- matrix(NA, self$p, q)
      }

      if (zero_inflation) private$fit_zi_component()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Setters    ------------------------
    #' @description
    #' Update the parameters of a model in place. Each argument applies to
    #' the families or variants named in parentheses; `NA` (the default)
    #' leaves the parameter unchanged.
    #' @param B regression coefficients (all models)
    #' @param dm1 inverse noise variances of the variables (variance-block)
    #' @param C cluster memberships: indicators, or posterior probabilities
    #' (all models)
    #' @param Omega precision matrix: between clusters (variance-block) or
    #' between variables (mean-block)
    #' @param gamma posterior variance of W (variance-block, known clusters)
    #' @param mu posterior mean of W (variance-block, known clusters)
    #' @param alpha cluster proportions (unknown clusters)
    #' @param M variational means of W (variance-block, unknown clusters)
    #' @param S variational variances of W (variance-block, unknown clusters)
    #' @param ll_list  list of log-lik (elbo) values
    #' @param warm_started whether `optim_initialize()` should treat the model as
    #' already initialized (reuse B/Omega/dm1/C/alpha/M/S as they stand)
    #' rather than recomputing a fresh heuristic initialization. Set by
    #' `$warm_start_from()` and by `$split()`/`$merge()`.
    #' @param clustering_init initial clustering
    #' @return Update the current object in place
    update = function(B = NA,
                      dm1 = NA,
                      C = NA,
                      Omega = NA,
                      gamma = NA,
                      mu = NA,
                      alpha = NA,
                      M = NA,
                      S = NA,
                      ll_list = NA,
                      warm_started = NA,
                      clustering_init = NA) {
      if (!anyNA(B))       private$B       <- B
      if (!anyNA(dm1))     private$dm1     <- dm1
      if (!anyNA(C))       private$C       <- C
      if (!anyNA(Omega))   private$Omega   <- Omega
      if (!anyNA(gamma))   private$gamma   <- gamma
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

    #' @description calls optimization (EM or heuristic) and updates relevant fields
    #' @param control a list for controlling the optimization process
    #' @param warn whether to warn when the (V)EM stops at the `niter` cap
    #' without reaching `threshold` (see `private$warn_if_not_converged()`).
    #' Set to `FALSE` for deliberately-truncated trial fits (cheap candidate
    #' scoring in `candidates_split()`/`candidates_merge()`, the sparsity-path
    #' warm-start probe in [NormalBlockCollectionSparsity]) where stopping at
    #' the cap is expected and not a sign of trouble.
    #' @return optimizes the model and updates its parameters
    optimize = function(control = list(niter = 500, threshold = 1e-4,
                                       fixed_point_niter = 5), warn = TRUE) {
      private$niter_max  <- control$niter
      private$threshold  <- control$threshold
      optim_out <- private$optimizer(control)
      do.call(self$update, optim_out)
      if (warn) private$warn_if_not_converged()
      invisible(self)
    },

    #' @description Try several clustering-initialization heuristics and keep
    #' the best-ELBO converged fit (see [NB_control()]'s `clustering_init`
    #' for the heuristics and the rationale behind the default). Every
    #' candidate is first screened with a short `trial_niter` run (same idea
    #' as `candidates_split()`/`candidates_merge()`), and only the
    #' `max_training` best-screened ones are fully retrained with `control`.
    #' @param inits vector of clustering-heuristic names to try; defaults to
    #' the model family's own preferred order (`private$default_inits`)
    #' @param trial_niter number of (V)EM iterations used to cheaply screen
    #' every candidate in `inits` before fully retraining the best few
    #' @param max_training how many of the screened candidates (best `loglik`
    #' after `trial_niter` iterations) get fully retrained with `control`
    #' @param control `optimize()` control list (`niter`/`threshold`) used
    #' for the final full retraining of the `max_training` best candidates
    #' @return a new, already-optimized model. Does not mutate the current
    #' object; reassign the result (`model <- model$best_of_inits()`).
    best_of_inits = function(inits = private$default_inits,
                             trial_niter = 10, max_training = 2,
                             control = list(niter = 500, threshold = 1e-4, fixed_point_niter = 5)) {
      stopifnot(
        "best_of_inits() only applies to an unknown clustering: a known one has no initialization to choose" =
          !private$known_clusters,
        "best_of_inits() requires (V)EM inference (not the heuristic-only mode, see NB_control(heuristic = ))" =
          !private$approx,
        "best_of_inits() only applies when the initial clustering is inferred by a heuristic (see NB_control(clustering_init = ))" =
          !is.na(private$clustering_approx)
      )
      candidates <- map(inits, function(init) {
        cand <- self$clone()
        cand$update(clustering_init = init)
        cand$optimize(private$trial_control(trial_niter), warn = FALSE)
        cand
      })
      ibest <- order(map_dbl(candidates, "loglik"), decreasing = TRUE)[1:min(max_training, length(candidates))]
      best_candidates <- candidates[ibest]
      map(best_candidates, function(cand) cand$optimize(control, warn = FALSE))
      best_candidates[[which.max(map_dbl(best_candidates, "loglik"))]]
    },

    #' @description Seed this model's starting parameters from another,
    #' already-optimized model of the same family and q, instead of a fresh
    #' heuristic clustering. Used by the sparsity paths to warm-start each
    #' penalty from the previous one's solution.
    #' @param other a model of the same class, already optimized
    #' @return Update the current object in place with `other`'s parameters
    warm_start_from = function(other) {
      stopifnot("warm_start_from() requires both models to have the same q" = self$q == other$q)
      args <- other$model_par
      if (!is.null(other$var_par)) {
        args$C     <- other$var_par$tau
        args$M     <- other$var_par$M
        args$S     <- other$var_par$S
        args$alpha <- colMeans(other$var_par$tau)
      } else if (!is.null(other$posterior_par)) {
        args$gamma <- other$posterior_par$gamma
        args$mu    <- other$posterior_par$mu
      }
      ## the zero-inflation component (B0, kappa) depends on the data only
      args <- Filter(Negate(is.null), args[intersect(names(args), names(formals(self$update)))])
      do.call(self$update, args)
      private$warm_started <- TRUE
      invisible(self)
    },

    #' @description generate and select a set of candidate models
    #' by splitting the clusters of the current model
    #' @param trial_niter number of (V)EM iterations used to cheaply score each
    #' candidate before the caller fully re-optimizes the best few -- kept
    #' short on purpose.
    #' @param split_rule passed to `split()` as its `rule`: `"kmeans"`,
    #' `"noise"`, or `NULL` for the family's default
    candidates_split = function(trial_niter = 5, split_rule = NULL) {
      # do not split groups with less than 2 guys
      candidates <- map((1:self$q)[self$cluster_sizes > 1],
                        function(k) tryCatch(self$split(k, rule = split_rule), error = function(e) NULL))
      candidates <- Filter(Negate(is.null), candidates)
      # keep candidates with at least 2 guys per cluster and non empty split.
      # Compared against the number of currently *live* (non-empty) clusters
      # rather than self$q: the (V)EM can converge to an empty cluster,
      # in which case self$q overcounts and every split would otherwise look
      # invalid, silently stalling explore_forward() before n_clusters_range[2].
      clustering_sizes <- map(map(candidates, "clustering"), table)
      min_sizes  <- map_dbl(clustering_sizes, min)
      n_clusters <- map_dbl(clustering_sizes, length)
      n_live <- length(unique(self$clustering))
      candidates <- candidates[min_sizes > 1 & n_clusters == n_live + 1]
      private$trial_fits(candidates, trial_niter)
    },

    #' @description generate and select a set of candidate models
    #' by merging the clusters of the current model
    #' @param max_candidates merge candidates are, unlike split's, quadratic
    #' in q (`choose(q, q-2)` pairs): beyond `max_candidates` pairs, only
    #' the most promising ones are actually built and trial-optimized, ranked
    #' by the family's own `private$merge_score()`. Set to `Inf` to always
    #' try every pair.
    #' @param trial_niter see `$candidates_split()`
    candidates_merge = function(max_candidates = 30, trial_niter = 2) {
      stopifnot("need at least two clusters to merge them" = self$q > 1)
      pairs <- combn(self$q, 2, simplify = FALSE)
      if (length(pairs) > max_candidates)
        pairs <- pairs[order(private$merge_score(pairs), decreasing = TRUE)[1:max_candidates]]
      candidates <- map(pairs, function(ij) tryCatch(self$merge(ij), error = function(e) NULL))
      private$trial_fits(Filter(Negate(is.null), candidates), trial_niter)
    },

    #' @description Log-likelihood of new observations at the fitted
    #' parameters and, for an unknown clustering, at the MAP partition (see
    #' [normal_block_cv()]). It is exact, computed for variance-block models
    #' with q x q matrices only (Woodbury identity), and given on Y's original
    #' scale.
    #' @param newdata a [NormalBlockData] object with the same variables and
    #' covariates (and zero-inflation covariates) as the training data.
    #' @param per_observation whether to return one value per observation
    #' rather than their sum.
    #' @return the held-out log-likelihood, or a vector of `newdata$n` values
    heldout_loglik = function(newdata, per_observation = FALSE) {
      stopifnot("newdata must be a NormalBlockData object" = inherits(newdata, "NormalBlockData"),
                "newdata must have the same variables as the training data" = newdata$p == self$p,
                "newdata must have the same covariates as the training data" = ncol(newdata$X) == ncol(self$data$X),
                "a heuristic fit (NB_control(heuristic = TRUE)) has no likelihood" = !isTRUE(private$approx))
      s <- self$data$Y_scale
      ## newdata's raw Y, put on the scale the model was fitted on
      Y <- newdata$Y_original / matrix(s, newdata$n, self$p, byrow = TRUE)
      C <- matrix(0, self$p, self$q)
      C[cbind(seq_len(self$p), self$clustering)] <- 1
      R <- Y - private$heldout_mean(newdata$X, C)
      ll <- if (is.matrix(private$kappa)) {
        stopifnot("newdata must have the same zero-inflation covariates" = ncol(newdata$X0) == nrow(private$B0))
        kappa <- if (all(private$kappa == 0)) matrix(0, newdata$n, self$p)
                 else stats::plogis(newdata$X0 %*% private$B0)
        observed <- newdata$zeros == 0
        private$heldout_rows(R, C, observed) +
          rowSums(ifelse(observed, log1p(-kappa), log(kappa)))
      } else {
        private$heldout_rows(R, C, NULL)
      }
      ll <- ll - sum(log(s))
      if (per_observation) ll else sum(ll)
    },

    #' @description Predicts the mean of Y for new covariates X, in Y's
    #' original units. For the variance-block family, `X B`: the latent block
    #' signals average out, so that on the training X it differs from
    #' [fitted()], which adds their posterior mean. The mean-block family
    #' overrides it.
    #' @param new_X new set of covariates.
    #' @return A n*p prediction matrix for new observations
    predict = function(new_X){
      private$rescale_to_original(new_X %*% private$B)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Extractors ------------------------
    #' @description Extract the network encoded by Omega, as a matrix rather
    #' than a plot (see `$plot_network()`): between clusters for the
    #' variance-block family, between variables for the mean-block family.
    #' @param type edge value in the network: "support" (binary edges),
    #' "precision" (entries of the precision matrix, on the scale of the fit)
    #' or "partial_cor" (partial correlations)
    #' @return a sparse square matrix, q x q (variance-block) or p x p
    #' (mean-block), with named rows and columns
    latent_network = function(type = c("partial_cor", "support", "precision")) {
      net <- switch(
        match.arg(type),
        "support"     = 1 * (private$Omega != 0 & !diag(TRUE, ncol(private$Omega))),
        "precision"   = private$Omega,
        "partial_cor" = {
          tmp <- -private$Omega / tcrossprod(sqrt(diag(private$Omega))); diag(tmp) <- 1
          tmp
        }
      )
      nodes <- network_nodes_names(self)
      if (length(nodes) == ncol(net)) dimnames(net) <- list(nodes, nodes)
      ## Sparse encoding: igraph::graph_from_adjacency_matrix() fails on dsyMatrix
      Matrix::Matrix(net, sparse = TRUE)
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Graphical methods------------------
    #' @param show_increment whether to add a second panel with the (log10)
    #' absolute increment between iterations and the convergence `threshold`
    #' (distinguishes true convergence from a flat-looking objective trace).
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

      ## log10(-obj) spreads out near-convergence iterations on the plot,
      ## same trick as the increment panel below; raw value if not all-negative.
      obj_is_neg <- all(obj < 0)
      obj_facet  <- if (obj_is_neg) "log10(-objective)" else "objective"
      obj_value  <- if (obj_is_neg) log10(-obj) else obj

      dplot <- data.frame(iteration = seq_along(obj), value = obj_value, facet = obj_facet)
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
        dplot <- rbind(
          dplot,
          data.frame(iteration = seq_along(obj), value = log10(abs(diff(ll))), facet = "log10(|increment|)")
        )
      }
      dplot$facet <- factor(dplot$facet, levels = c(obj_facet, "log10(|increment|)"))

      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = .data$iteration, y = .data$value)) +
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
            ggplot2::aes(yintercept = .data$yintercept),
            linetype = "dashed", colour = "red", alpha = .6
          )
      }
      p
    },

    #' @description Plot the network encoded by Omega (see
    #' `$latent_network()` to extract it as a matrix).
    #' @param type edge value in the network: "partial_cor" (partial
    #' correlations) or "support" (binary edges).
    #' @param output Output type. Either `igraph` (for the network) or `corrplot` (for the adjacency matrix)
    #' @param edge.color Length 2 color vector. Color for positive/negative edges. Default is `c("#F8766D", "#00BFC4")`. Only relevant for igraph output.
    #' @param node.labels vector of character. The labels of the nodes. By
    #' default, the cluster names (variance-block) or the variable names
    #' (mean-block), as in `$latent_network()`.
    #' @param remove.isolated if `TRUE`, isolated node are remove before plotting. Only relevant for igraph output.
    #' @param layout an igraph layout function, applied to the graph without its edge weights, or a matrix of node coordinates. Only relevant for igraph output.
    #' @param plot logical. Should the final network be displayed or only sent back to the user. Default is `TRUE`.
    plot_network = function(type            = c("partial_cor", "support"),
                            output          = c("igraph", "corrplot"),
                            edge.color      = c("#F8766D", "#00BFC4"),
                            remove.isolated = FALSE,
                            node.labels     = NULL,
                            layout          = igraph::layout_in_circle,
                            plot = TRUE) {
      if(anyNA(private$Omega)) stop("NA in the precision matrix")

      type   <- match.arg(type)
      output <- match.arg(output)
      net <- self$latent_network(type)

      if (output == "igraph") {
        G <-  igraph::graph_from_adjacency_matrix(net, mode = "undirected", weighted = TRUE, diag = FALSE)

        igraph::V(G)$label <- if (!is.null(node.labels)) node.labels else rownames(net)
        ## Nice nodes: a cluster's area grows with its number of variables
        sizes <- if (inherits(self, "NormalBlockVarBase")) tabulate(self$clustering, nbins = ncol(net))
                 else rep(1, ncol(net))
        igraph::V(G)$size <- 6 + 14 * sqrt(sizes / max(sizes))
        igraph::V(G)$color <- "grey90"
        igraph::V(G)$frame.color <- "grey40"
        igraph::V(G)$label.cex <- 0.8
        igraph::V(G)$label.color <- rgb(0, 0, .2, .8)
        ## Nice edges
        igraph::E(G)$color <- ifelse(igraph::E(G)$weight > 0, edge.color[1], edge.color[2])
        if (type == "support")
          igraph::E(G)$width <- abs(igraph::E(G)$weight)
        else
          igraph::E(G)$width <- 10*abs(igraph::E(G)$weight)

        if (remove.isolated) {
          G <- igraph::delete_vertices(G, which(igraph::degree(G) == 0))
        }
        if (plot) {
          ## layouts are computed without the (signed) edge weights
          if (is.function(layout))
            layout <- layout(igraph::delete_edge_attr(G, "weight"))
          plot(G, layout = layout)
        }
      }
      if (output == "corrplot") {
        G <- net
        if (plot) {
          if (ncol(net) > 100)
            colnames(net) <- rownames(net) <- rep(" ", ncol(net))
          diag(net) <- 0
          corrplot::corrplot(as.matrix(net), method = "color", is.corr = FALSE, tl.pos = "td", cl.pos = "n", tl.cex = 0.5, type = "upper")
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
    print = function(model = paste0("A ", self$who_am_I, ".\n")) {
      cat(model)
      cat("===========================================================================\n")
      print(as.data.frame(round(self$criteria, digits = 3), row.names = ""))
      cat("===========================================================================\n")
      cat("* Useful fields\n")
      cat("    $model_par, $posterior_par / $var_par, $clustering \n")
      cat("    $loglik, $BIC, $ICL, $objective, $nb_param, $criteria\n")
      cat("* Useful S3 methods\n")
      cat("    print(), summary(), plot(), coef(), sigma(), fitted(), predict() \n")
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    ## Groups of free parameters and the number of independent observations
    ## informing each one (see the penalty_terms field): the family's own
    ## groups, then the zero-inflation coefficients (d0 per variable, informed
    ## by its n observations) and the cluster proportions (informed by the p
    ## variables) when the model has them.
    penalty_groups = function() {
      rbind(private$family_penalty_groups(),
            if (private$zi) penalty_group("kappa", self$p * self$d0, self$n),
            if (!private$known_clusters) penalty_group("alpha", self$q - 1, self$p))
    },
    family_penalty_groups = function() penalty_group(character(0), integer(0), numeric(0)),
    ## Dimension of the network the EBIC is computed on. Overridden by the
    ## mean-block family, whose network relates the p variables.
    network_dim = function() self$q,
    family_name       = NA, # "normal-block-var" or "normal-block-mean", for who_am_I
    known_clusters    = FALSE, # TRUE for the classes with a known clustering
    zi                = FALSE, # whether the model is zero-inflated
    res_covariance    = NA, # shape of the residual covariance
    kappa             = NA, # matrix of zero-inflation probabilities
    B0                = NA, # zero-inflation regression matrix
    ZI_cond_mean      = NA, # fixed contribution of the ZI component to the log-likelihood
    B                 = NA, # regression matrix
    C                 = NA, # the matrix of posterior probabilities (tau) or group affectation
    Omega             = NA, # precision matrix for clusters or variables
    alpha             = NA, # vector of clusters probabilities
    optimizer         = NA, # a link to the function that perform the optimization
    ll_list           = NA, # list of log-likelihoods or ELBOs
    sparsity_         = NA, # scalar controlling the overall sparsity
    weights           = NA, # sparsity weights specific to each pairs of group
    approx            = NA, # use approximation/heuristic approach or not
    clustering_approx = NA, # name of the clustering heuristic, key into clustering_methods
    niter             = NA, # number of EM iterations required by the inference, if applicable
    niter_max         = NA, # niter cap passed to the last optimize() call (for plot_loglik()'s
    # and warn_if_not_converged()'s diagnostic)
    threshold         = NA, # convergence threshold passed to the last optimize() call (idem)
    warm_started      = FALSE, # set by warm_start_from() and by split()/merge(): tells
    # optim_initialize() to reuse the current B/dm1/Omega (and
    # C/M/S/alpha or gamma/mu) instead of (re-)deriving them
    # from the heuristic clustering.
    # Deliberately NOT inferred from "are these fields non-NA",
    # which would also be true for split()/merge() clones (those
    # keep their own, different, already-tested initialization path).

    ## Control list for the deliberately-truncated trial fits of
    ## best_of_inits()/candidates_split()/candidates_merge(). fixed_point_niter
    ## is only read by the mean-block family; the variance-block one ignores it.
    trial_control = function(trial_niter) {
      list(niter = trial_niter, threshold = 1e-4, fixed_point_niter = 5)
    },

    ## Short (V)EM runs of the split()/merge() candidates; a candidate whose
    ## run fails (e.g. a degenerate block) is discarded rather than aborting
    ## the search.
    trial_fits = function(candidates, trial_niter) {
      ok <- vapply(candidates, function(cand) {
        tryCatch({
          cand$optimize(private$trial_control(trial_niter), warn = FALSE)
          TRUE
        }, error = function(e) FALSE)
      }, logical(1))
      candidates[ok]
    },

    ## Shared tail of split() and merge(): applies the new clustering and its
    ## reseeded parameters to the model itself or to a clone, marked as warm
    ## started so that optim_initialize() reuses them.
    reseeded = function(in_place, C, ..., weights = NULL) {
      target <- if (in_place) self else self$clone()
      target$update(C = C, ..., alpha = colMeans(C), warm_started = TRUE)
      if (!is.null(weights)) target$sparsity_weights <- weights
      invisible(target)
    },

    ## Memberships after splitting cluster `index`: `split2` moves to a new
    ## cluster q + 1, `split1` stays, both kept slightly away from 0.
    split_memberships = function(index, split1, split2) {
      new_C <- cbind(private$C, .Machine$double.eps)
      new_C[split1, index] <- new_C[split1, index] - .Machine$double.eps
      new_C[split2, self$q + 1] <- new_C[split2, index]
      new_C[split2, index] <- .Machine$double.eps
      new_C / rowSums(new_C)
    },

    ## Memberships after merging the clusters `indices` (sorted) into the first
    merge_memberships = function(indices) {
      new_C <- private$C[, -indices[2], drop = FALSE]
      new_C[, indices[1]] <- private$C[, indices[1]] + private$C[, indices[2]]
      new_C
    },

    ## Ranks candidate cluster pairs for candidates_merge(), highest first.
    ## Overridden by each family (see NormalBlockVarBase/NormalBlockMeanBase).
    merge_score = function(pairs) rep(0, length(pairs)),
    ## family-specific default of split()'s rule
    default_split_rule = "noise",

    ## Which members of cluster `index` split() moves to the new cluster.
    ## "kmeans": a 2-means on the members' columns of the family's
    ## split_target() (the data its clustering is about); "noise": the
    ## members whose noise_variance() is at most the cluster's median. The
    ## 2-means falls back to the noise rule if it cannot produce two
    ## non-empty groups (e.g. fewer than 2 distinct columns).
    split_members = function(index, rule = NULL) {
      if (is.null(rule)) rule <- private$default_split_rule
      rule <- match.arg(rule, c("kmeans", "noise"))
      cl   <- self$clustering == index
      if (rule == "kmeans" && sum(cl) > 1) {
        target <- private$split_target()[, cl, drop = FALSE]
        km <- tryCatch(stats::kmeans(t(target), centers = 2, nstart = 10),
                       error = function(e) NULL)
        if (!is.null(km) && length(unique(km$cluster)) == 2) {
          moved <- logical(self$p)
          moved[cl] <- km$cluster == 2
          return(moved)
        }
      }
      var <- private$noise_variance()
      (var <= median(var[cl])) & cl
    },

    ## Warns when the (V)EM hit the niter cap without reaching threshold
    ## (heuristic fits have no ll_list/niter to check, so stay silent; a
    ## negative threshold asks for exactly niter iterations).
    warn_if_not_converged = function() {
      if (is.na(private$niter) || private$niter < private$niter_max || private$threshold < 0)
        return(invisible())
      last_increment <- abs(diff(private$ll_list))[private$niter]
      if (last_increment >= private$threshold)
        warning(sprintf(
          "%s: (V)EM stopped at the niter cap (%d) without reaching the convergence threshold (last increment = %.3g, threshold = %.1e). Consider raising `niter` in NB_control(), especially with many clusters -- see plot_loglik() to check.",
          self$who_am_I, private$niter_max, last_increment, private$threshold), call. = FALSE)
      invisible()
    },

    ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Zero-inflation probabilities, never revisited by the (V)EM: only their
    ## fixed log-likelihood contribution (ZI_cond_mean) reaches the recursion.
    ## They depend on nothing but the data, so the fit itself lives on NormalBlockData
    fit_zi_component = function() {
      zi <- self$data$zi_fit()
      private$B0           <- zi$B0
      private$kappa        <- zi$kappa
      private$ZI_cond_mean <- zi$ZI_cond_mean
    },

    ## Methods for integrated (V)EM inference --------------
    ## Each concrete subclass overrides EM_optimize() to call its
    ## Rcpp/Armadillo core (see src/exports.cpp and
    ## inst/normal_block_models.qmd); optim_initialize() (heuristic
    ## initialization) stays in R and supplies the starting values.
    ## MLE of MV Normal distribution
    multivariate_normal_inference = function() self$data$ols_fit(),

    optim_initialize = function() {},

    ## Converts a per-variable quantity from the internal (rescaled) fitting
    ## scale back to Y's original units: power = 1 for additive quantities
    ## (B, fitted values), power = -2 for inverse-variance quantities (dm1).
    rescale_to_original = function(M, power = 1) {
      factor <- self$data$Y_scale^power
      if (is.matrix(M)) M * matrix(factor, nrow(M), ncol(M), byrow = TRUE) else M * factor
    },

    ## Precision matrix from a covariance estimate: plain inversion, or
    ## graphical lasso when penalized
    get_Omega = function(Sigma, weights = self$sparsity_weights) {
      if (private$sparsity_ == 0) {
        chol2inv(chol(Sigma))
      } else {
        glasso_omega(Sigma, private$sparsity_ * weights)
      }
    },


    ## Registry of clustering heuristics turning the residuals R (n x p) into
    ## an initial clustering of the p variables into self$q groups, selected
    ## by NB_control(clustering_init = ...); kmeans is the default.
    ## spectral clusters the eigenvectors of cov(R) (top q, each row rescaled
    ## to unit L2 norm, the classic Ng-Jordan-Weiss normalization).
    ## Each one is the collection-level path (R/utils.R) for a single q.
    clustering_methods = list(
      kmeans   = function(R, q) kmeans_clustering_path(R, q)[[1]],
      ward2    = function(R, q) ward2_clustering_path(R, q)[[1]],
      sbm      = function(R, q) sbm_clustering_path(R, q)[[1]],
      spectral = function(R, q) spectral_clustering_path(R, q)[[1]]
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
    #' @field d0 number of zero-inflation covariates (dimensions in X0)
    d0 = function() self$data$d0,
    #' @field who_am_I a description of the model
    who_am_I = function() {
      paste(c(if (private$zi) "zero-inflated", private$res_covariance, private$family_name, "model with",
              if (private$known_clusters) "fixed clusters" else paste(self$q, "unknown clusters")),
            collapse = " ")
    },
    #' @field q number of blocks
    q = function() as.integer(ncol(private$C)),
    #' @field n_edges number of edges of the network (non null coefficient of the sparse precision matrix Omega)
    n_edges  = function() sum(private$Omega[upper.tri(private$Omega, diag = FALSE)] != 0),
    #' @field objective evolution of the objective function during (V)EM algorithm
    objective = function() private$ll_list[-1],
    #' @field loglik (or its variational lower bound)
    loglik = function() if (private$approx) NA else private$ll_list[[length(private$ll_list)]] + self$sparsity_term,
    #' @field deviance (or its variational lower bound)
    deviance = function() -2 * self$loglik,
    #' @field entropy entropy of the (variational) posterior distribution of
    #' the clustering, `-sum(tau * log(tau))`: 0 when the clustering is
    #' known. Only the discrete partition is completed, as in the ICL of
    #' Biernacki, Celeux and Govaert (2000); the latent vectors W of the
    #' variance-block family are integrated out, not completed.
    entropy    = function() {
      if (private$known_clusters) 0 else if (private$approx) NA else -sum(xlogx(private$C))
    },
    #' @field penalty_terms data frame with one row per group of parameters:
    #' its name (`group`), its number of free parameters (`nu`) and the
    #' number of independent observations informing it (`N`), so that the
    #' group costs `nu * log(N)` in the BIC. As in the ICL of Daudin, Picard
    #' and Robin (2008) for SBMs and Keribin et al. (2015) for latent block
    #' models, the cluster proportions are charged `log(p)` (the variables are
    #' the units being clustered), parameters shared by every cell of a block
    #' `log(n * p)`, and the remaining ones `log(n)`.
    penalty_terms = function() private$penalty_groups(),
    #' @field nb_param number of free parameters in the model, `sum(penalty_terms$nu)`
    nb_param = function() as.integer(sum(private$penalty_groups()$nu)),
    #' @field penalty BIC penalty, `sum(nu * log(N))` over `penalty_terms`
    penalty = function() {
      pen <- private$penalty_groups()
      sum(pen$nu * log(pen$N))
    },
    #' @field BIC Bayesian information criterion, `deviance + penalty` (lower
    #' is better); built on the ELBO when the clustering is unknown and the
    #' likelihood is not available in closed form.
    BIC = function() self$deviance + self$penalty,
    #' @field ICL integrated completed likelihood criterion, `BIC + 2 * entropy`
    #' (lower is better); equal to the BIC when the clustering is known.
    ICL        = function() self$BIC + 2 * self$entropy,
    #' @field EBIC extended BIC of Foygel and Drton (2010) with gamma = 1/2,
    #' `BIC + 2 * n_edges * log(dim)`, where `dim` is the dimension of the
    #' network: q for the variance-block family (network between clusters), p
    #' for the mean-block family (network between variables).
    EBIC   = function() self$BIC + 2 * ifelse(self$n_edges > 0, self$n_edges * log(private$network_dim()), 0),
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
    #' @field sparsity_term (sparsity_term term in log-likelihood due to sparsity)
    sparsity_term = function() self$sparsity * sum(abs(self$sparsity_weights * private$Omega)),
    #' @field memberships the p x q matrix of cluster memberships: indicators
    #' for a known clustering, posterior probabilities (tau) otherwise
    memberships = function() private$C,
    #' @field clustering the cluster label of each variable (MAP for an
    #' unknown clustering), named after the columns of Y
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
