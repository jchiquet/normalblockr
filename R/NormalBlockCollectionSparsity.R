## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionSparsity ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection of Models over a Sparsity Path
#'
#' Shared scaffolding for [NormalBlockVarCollectionSparsity] and
#' [NormalBlockMeanCollectionSparsity]: the penalty grid, the warm-started
#' path traversal, the penalty lookup and the criteria plot. Model selection,
#' StARS included (stability selection, Liu et al. 2010), is shared: a family
#' only says on which covariance estimate the default penalty grid is read
#' (`private$penalty_scale()`), how to refit its path on a subsample with the
#' clustering held fixed (`private$subsample_path()`) and how to name the
#' nodes of its network (`private$network_nodes()`).
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how collections are created and fitted.
#' @keywords internal
NormalBlockCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockCollectionSparsity",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data = NA,

    #' @description Create a collection over a path of sparsity penalties.
    #' By default, `control$n_sparsity_penalties` values on a log scale, from
    #' the largest absolute (weighted) entry of a covariance estimate, which
    #' removes every edge, down to `control$min_ratio` times it.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param blocks either a clustering matrix (known, fixed clustering) or a single integer (number of blocks to infer)
    #' @param zero_inflation whether the models of the path are zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new collection object
    initialize = function(mydata, blocks, zero_inflation = FALSE, control = NB_control()) {
      self$data    <- mydata
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      private$blocks_ <- blocks
      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)

      if (!is.null(control$sparsity_penalties)) {
        stopifnot("All penalties must be strictly positive" =
                    (min(control$sparsity_penalties) > 0))
        sparsity <- control$sparsity_penalties
      } else {
        scale <- private$penalty_scale(mydata, blocks, zero_inflation, control)
        ## a zero weight leaves its entry unpenalized
        scale <- scale[is.finite(scale)]
        if (length(scale) > 0) {
          max_pen  <- max(scale)
          sparsity <- 10^seq(log10(max_pen), log10(max_pen * control$min_ratio),
                             len = control$n_sparsity_penalties)
        } else {
          sparsity <- rep(0, len = control$n_sparsity_penalties)
        }
      }

      private$sparsity_ <- sort(sparsity, decreasing = TRUE)
      self$models <- map(private$sparsity_, function(lambda)
        get_model(mydata, blocks, lambda, zero_inflation, control, model = private$family))
      private$progress_field <- "sparsity"
      private$progress_label <- "penalty"
    },

    #' @description optimizes every model in the sparsity path, warm-starting
    #' each one (after the first) from the previous, adjacent penalty's
    #' converged parameters (see the family's base class's `warm_start_from()`)
    #' instead of re-deriving everything from the heuristic clustering, the
    #' way the generic [NormalBlockCollection] `optimize()` would. `blocks`
    #' (hence q) is fixed across the whole path, only the sparsity penalty
    #' changes, so the warm start is always between models of matching shape.
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      previous <- NULL
      self$models <- lapply(self$models, function(model) {
        ## a line rewritten in place only makes sense in a console
        if (control$verbose && interactive()) {
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
          flush.console()
        }
        if (!is.null(previous)) model$warm_start_from(previous)
        model$optimize(control)
        previous <<- model
        model
      })
      invisible(self)
    },

    #' @description returns the NormalBlockVarKnownClusters model corresponding to given penalty
    #' @param sparsity sparsity penalty asked by user
    #' @return A NormalBlockVarKnownClusters (sparse) object with given value penalty
    get_model = function(sparsity) {
      if (!(sparsity %in% private$sparsity_)) {
        sparsity <-  private$sparsity_[[which.min(abs(private$sparsity_ - sparsity))]]
        message("No model with this penalty in the collection. Returning model with closest penalty: ",
                sparsity,  " Collection penalty values can be found via $sparsity")
      }
      self$models[[which(private$sparsity_ == sparsity)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to perform the
    #' selection: "BIC" (default), "EBIC", "ICL" or "StARS".
    #' @param stability if `crit = "StARS"`, the level of stability required.
    #' @return a model of the collection's family
    get_best_model = function(crit = c("BIC", "EBIC", "ICL", "StARS"), stability = 0.9) {
      crit <- match.arg(crit)
      if (crit != "StARS") return(self$models[[private$best_id(crit)]]$clone())
      if (is.null(private$stab_path)) self$stability_selection()
      max_stab <- max(self$stability)
      if (max_stab < stability) {
        message("No model reaches the required stability ", stability,
                ", returning model with highest stability: ", max_stab)
        stability <- max_stab
      }
      ## the sparsest network whose edges are at least that stable
      id <- match(min(private$sparsity_[self$stability >= stability]), private$sparsity_)
      self$models[[id]]$clone()
    },

    #' @description Computes the stability path (StARS, Liu et al. 2010): the
    #' penalty path is refitted on subsamples of the observations, with the
    #' clustering held fixed, and the selection frequency of each edge is
    #' recorded at every penalty.
    #' @param subsamples a list of vectors describing the subsamples. The number of vectors (or list length) determines the number of subsamples used in the stability selection. Automatically set to 10 subsamples with size `10*sqrt(n)` if `n >= 144` and `0.8*n` otherwise following Liu et al. (2010) recommendations.
    #' @param n_subsamples number of subsamples to create if the subsamples are not given
    stability_selection = function(subsamples = NULL, n_subsamples = 10) {
      if (is.null(subsamples)) {
        size <- round(ifelse(self$data$n >= 144, 10 * sqrt(self$data$n), 0.8 * self$data$n))
        subsamples <- replicate(n_subsamples, sample.int(self$data$n, size), simplify = FALSE)
      }
      control_stabs <- self$control
      control_stabs$sparsity_penalties <- private$sparsity_
      control_stabs$verbose <- FALSE
      if (self$control$verbose) cat("\nStability selection over", length(subsamples), "subsamples\n")

      nodes <- private$network_nodes()
      upper <- upper.tri(diag(length(nodes)))
      nets <- lapply(subsamples, function(rows) {
        sub_data <- self$data$subset(rows, check_constant = FALSE)
        path <- private$subsample_path(sub_data, control_stabs)
        path$optimize(control_stabs)
        do.call(cbind, lapply(path$models, function(m) as.matrix(m$latent_network("support"))[upper]))
      })
      prob <- Reduce("+", nets) / length(subsamples)
      colnames(prob) <- private$sparsity_
      edges <- which(upper, arr.ind = TRUE)
      edges <- edges[order(edges[, "col"], edges[, "row"]), , drop = FALSE]  # column-major, as `[upper]`
      node1 <- rep(nodes[edges[, "row"]], times = ncol(prob))
      node2 <- rep(nodes[edges[, "col"]], times = ncol(prob))
      private$stab_path <- data.frame(
        Edge    = paste0(node1, "|", node2),
        Penalty = rep(private$sparsity_, each = nrow(prob)),
        Prob    = as.vector(prob),
        Node1   = node1,
        Node2   = node2
      )
      invisible(subsamples)
    },

    #' @description Plot the criteria of the models of the path against the
    #' penalty, with a dashed line at the penalty selected by each criterion
    #' (see also [plot.NormalBlockCollection()]).
    #' @param criteria the criteria to plot, among `"deviance"`, `"BIC"`,
    #' `"EBIC"` and `"ICL"`; all of them by default.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "EBIC", "ICL"), log.x = TRUE) {
      stopifnot(!anyNA(self$criteria[criteria]))
      p <- private$plot_criteria_path("sparsity", criteria)
      if (log.x) p <- p + ggplot2::scale_x_log10()
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    sparsity_ = NA,   # penalty values in the collection (lambda)
    blocks_   = NA,   # blocks (either a scalar or an indicator matrix)
    scope = function() paste("over sparsity penalties, with",
                             if (is.matrix(private$blocks_)) "fixed clusters" else "fixed q"),
    stab_path = NULL, # stability path, filled in by stability_selection()
    family    = NA_character_, # "var" or "mean", set by the subclasses
    ## family hook for the default penalty grid: the absolute entries of a
    ## covariance estimate divided by their sparsity weights
    penalty_scale  = function(data, blocks, zero_inflation, control) stop("not implemented for this family"),
    ## family hooks for stability_selection(): the same penalty path refitted
    ## on a subsample with the clustering held fixed, and the network's nodes
    subsample_path = function(data, control) stop("not implemented for this family"),
    network_nodes  = function() stop("not implemented for this family")
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field sparsity_details number, range and ratio of the penalties of the path
    sparsity_details = function()
      list("n_penalties" = length(self$sparsity), "min_ratio" = min(self$sparsity)/max(self$sparsity),
           "min_sparsity" = min(self$sparsity), "max_sparsity" = max(self$sparsity)),
    #' @field q number of blocks
    q = function() ifelse(is.matrix(private$blocks_), ncol(private$blocks_), private$blocks_),
    #' @field blocks group matrix or number of blocks
    blocks = function() private$blocks_,
    #' @field sparsity list of sparsity penalties
    sparsity = function() private$sparsity_,
    #' @field criteria a data frame with the criteria of every model of the
    #' path, plus its StARS stability once `stability_selection()` has run
    criteria = function() {
      res <- super$criteria
      res$stability <- self$stability
      res
    },
    #' @field stability_path selection probability of every edge at every
    #' penalty (StARS), `NULL` until `stability_selection()` has run
    stability_path = function() private$stab_path,
    #' @field stability mean edge stability along the penalty path,
    #' `1 - mean(4 p (1 - p))` over the edges; `NA` until
    #' `stability_selection()` has run
    stability = function() {
      if (is.null(private$stab_path)) return(rep(NA_real_, length(private$sparsity_)))
      sapply(private$sparsity_, function(lambda) {
        prob <- private$stab_path$Prob[private$stab_path$Penalty == lambda]
        1 - mean(4 * prob * (1 - prob))
      })
    }
  )
)
