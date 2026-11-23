## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection of Models over a Range of Cluster Counts
#'
#' Shared scaffolding for [NormalBlockVarCollectionClusters] and
#' [NormalBlockMeanCollectionClusters]: the models of the collection, one per
#' q in increasing order, model lookup and selection, the criteria plot and
#' the split/merge `refine()` search. Concrete subclasses only name the model
#' family.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how collections are created and fitted.
#' @keywords internal
NormalBlockCollectionClusters <- R6::R6Class(
  classname = "NormalBlockCollectionClusters",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a collection over a range of q, sorted increasingly.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list values of q (number of clusters) in the collection
    #' @param zero_inflation whether the models of the collection are zero-inflated
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new collection object
    initialize = function(mydata, q_list, zero_inflation = FALSE,
                          sparsity = 0, control = NB_control()) {
      check_q_list(q_list, mydata$p)
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      ord <- order(q_list)
      inits <- clustering_inits(mydata, q_list, private$family, zero_inflation, control)[ord]
      self$models <- map2(q_list[ord], inits, function(q, init) {
        control$clustering_init <- init
        get_model(mydata, q, sparsity = sparsity, zero_inflation = zero_inflation,
                  model = private$family, control = control)
      })
      private$progress_field <- "q"
      private$progress_label <- "number of clusters"
    },

    #' @description The model of the collection with q clusters.
    #' @param q a number of clusters among those of the collection (`$q_list`)
    #' @return a fitted model, of the class of the collection's models
    get_model = function(q) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      self$models[[which(self$q_list == q)]]
    },

    #' @description The best model of the collection for a criterion.
    #' @param crit the criterion minimized: "ICL" (the default), "BIC",
    #' "EBIC", or "deviance" (which always selects the largest q)
    #' @return a copy of the selected fitted model
    get_best_model = function(crit = c("ICL", "BIC", "EBIC", "deviance")) {
      crit <- match.arg(crit)
      id <- private$best_id(crit)
      self$models[[id]]$clone()
    },

    #' @description Plot the criteria of the models of the collection
    #' against q, with a dashed line at the q selected by each criterion (see
    #' also [plot.NormalBlockCollection()]).
    #' @param criteria the criteria to plot, among `"deviance"`, `"ICL"`,
    #' `"BIC"` and `"EBIC"`; all of them by default.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "ICL", "BIC", "EBIC")) {
      stopifnot(!anyNA(self$criteria[criteria]))
      private$plot_criteria_path("q", criteria)
    },

    #' @description optimizes every model in the collection, then -- if
    #' `control$refine` is `TRUE` (see [NB_control()], the default) --
    #' calls `$refine()` automatically, reporting in one line how many models
    #' it improved.
    #' @param control optimization parameters (niter, threshold, verbose)
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      super$optimize(control)
      ## heuristic fits have no likelihood for refine() to compare
      if (isTRUE(self$control$refine) && !isTRUE(self$control$heuristic)) {
        before <- self$criteria$deviance
        self$refine(verbose = FALSE)
        if (isTRUE(self$control$verbose))
          cat("\t refine(): improved", sum(self$criteria$deviance < before),
              "of", length(before), "models\n")
      }
      invisible(self)
    },

    #' @description Tries to improve every model in the collection with a
    #' short split-and-reoptimize trial seeded from its smaller-q neighbor
    #' (`"split"`), a short merge-and-reoptimize trial seeded from its
    #' larger-q neighbor (`"merge"`), or both (the default); a candidate
    #' replaces the original only if it strictly lowers the deviance, so this
    #' can only improve (or leave unchanged) each model it touches. Only
    #' contiguous q pairs (`q` and `q -/+ 1`, both present in the collection)
    #' are refined. Run by default after fitting a collection, see
    #' [NB_control()]'s `refine` for the rationale.
    #' @param trial_niter number of EM iterations used for the cheap trial
    #' candidates (passed to `candidates_split()`/`candidates_merge()`)
    #' before fully re-optimizing only the best one.
    #' @param max_candidates passed to `candidates_merge()` (ignored for
    #' `"split"`, which is never combinatorial in q) -- see its
    #' documentation.
    #' @param directions which neighbor(s) to seed refinement candidates
    #' from: `"split"` (smaller-q neighbor), `"merge"` (larger-q neighbor),
    #' or both (the default).
    #' @param split_rule how a cluster is split to seed the q + 1 model, see
    #' the `rule` argument of the models' `split()`: `"kmeans"`, `"noise"`,
    #' or `NULL` (default) for the family's own default (`"kmeans"` for the
    #' variance-block family, `"noise"` for the mean-block one).
    #' @param verbose whether to print, for each q attempted, whether the
    #' candidate from that neighbor improved on it. Defaults to
    #' `control$verbose` (the value set at construction, see
    #' [NB_control()]).
    #' @return invisibly returns `self`; improved models replace the
    #' originals in `$models` in place.
    refine = function(trial_niter = 2, max_candidates = 30,
                      directions = c("split", "merge"), verbose = self$control$verbose,
                      split_rule = NULL) {
      stopifnot(all(directions %in% c("split", "merge")))
      stopifnot("refine() compares deviances, which heuristic fits (NB_control(heuristic = TRUE)) do not have" =
                  !isTRUE(self$control$heuristic))
      private$refined <- TRUE
      if (!is.null(split_rule)) split_rule <- match.arg(split_rule, c("kmeans", "noise"))
      ord <- order(self$q_list)
      q_sorted <- self$q_list[ord]
      n <- length(ord)

      try_candidate <- function(idx_target, candidates, label) {
        if (length(candidates) == 0) return(invisible())
        best_candidate <- candidates[[which.min(map_dbl(candidates, "ICL"))]]
        if (inherits(try(best_candidate$optimize(), silent = TRUE), "try-error")) return(invisible())
        deviance_before <- self$models[[idx_target]]$deviance
        if (best_candidate$deviance < deviance_before) {
          if (verbose) cat("\t refine: q =", self$models[[idx_target]]$q, "-- deviance", round(deviance_before, 2),
                            "->", round(best_candidate$deviance, 2), "(improved via", label, ")\n")
          self$models[[idx_target]] <- best_candidate
        } else if (verbose) {
          cat("\t refine: q =", self$models[[idx_target]]$q, "-- no improvement found (", label, ")\n")
        }
      }

      if ("split" %in% directions) {
        for (i in 2:n) {
          if (q_sorted[i] != q_sorted[i - 1] + 1) next # non-contiguous, nothing to seed from
          candidates <- self$models[[ord[i - 1]]]$candidates_split(trial_niter = trial_niter, split_rule = split_rule)
          try_candidate(ord[i], candidates, "split")
        }
      }
      if ("merge" %in% directions) {
        for (i in (n - 1):1) {
          if (q_sorted[i] != q_sorted[i + 1] - 1) next # non-contiguous, nothing to seed from
          candidates <- self$models[[ord[i + 1]]]$candidates_merge(max_candidates = max_candidates, trial_niter = trial_niter)
          try_candidate(ord[i], candidates, "merge")
        }
      }
      invisible(self)
    },

    #' @description User-friendly print method: as for any collection, plus a
    #' warning line when the deviance increases with q somewhere (see
    #' `deviance_violations`).
    print = function() {
      super$print()
      viol <- self$deviance_violations
      if (length(viol) > 0)
        cat("* Warning: the deviance increases with q at q =", paste(viol, collapse = ", "),
            "\n    (nested models: a sign of poor local optima",
            if (private$refined) "that refine() did not repair).\n" else "). Try $refine().\n")
      invisible(self)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    family = NA_character_, # "var" or "mean", set by the subclasses
    refined = FALSE,        # whether refine() has run, for print()
    scope = function() "over q"
  ),
  active = list(
    #' @field q_list number of blocks
    q_list = function() map_dbl(self$models, "q"),
    #' @field deviance_violations the values of q whose deviance is higher
    #' than at the previous q of the collection. The models being nested, each
    #' one signals a poor local optimum, which `refine()` can repair.
    deviance_violations = function() {
      ord <- order(self$q_list)
      dev <- self$criteria$deviance[ord]
      self$q_list[ord][-1][which(diff(dev) > 1e-6 * pmax(1, abs(dev[-1])))]
    }
  )
)
