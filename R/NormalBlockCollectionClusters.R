## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection of Models over a Range of Cluster Counts
#'
#' Shared scaffolding for [NormalBlockVarCollectionClusters] and
#' [NormalBlockMeanCollectionClusters]: everything that does not depend on the
#' model family (model lookup, model selection, the criteria plot and the
#' split/merge `refine()` search). Concrete subclasses only build `self$models`
#' in their `initialize()` and name themselves through `who_am_I`.
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
    #' @description returns the unknown-clusters model corresponding to given q
    #' @param q number of blocks asked by user
    #' @return A unknown-clusters object with given value q
    get_model = function(q) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      self$models[[which(self$q_list == q)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "ICL" or "BIC". "ICL" is the default criterion
    #' @return a [`unknown-clusters`] object
    get_best_model = function(crit = c("ICL", "BIC", "EBIC", "deviance")) {
      crit <- match.arg(crit)
      id <- private$best_id(crit)
      self$models[[id]]$clone()
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a [`Networkfamily`])
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", "BIC", "ICL")`. Defaults to all of them.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "ICL", "BIC", "EBIC")) {
      stopifnot(!anyNA(self$criteria[criteria]))
      private$plot_criteria_path("q", criteria)
    },

    #' @description optimizes every model in the collection, then -- if
    #' `control$refine` is `TRUE` (see [NB_control()], default `FALSE`) --
    #' calls [refine()] automatically.
    #' @param control optimization parameters (niter, threshold, verbose)
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      super$optimize(control)
      if (isTRUE(self$control$refine)) self$refine()
      invisible(self)
    },

    #' @description Tries to improve every model in the collection with a
    #' short split-and-reoptimize trial seeded from its smaller-q neighbor
    #' (`"split"`), a short merge-and-reoptimize trial seeded from its
    #' larger-q neighbor (`"merge"`), or both (the default); a candidate
    #' replaces the original only if it strictly lowers the deviance, so this
    #' can only improve (or leave unchanged) each model it touches. Only
    #' contiguous q pairs (`q` and `q -/+ 1`, both present in the collection)
    #' are refined. See `inst/methods_initialization_and_refine.md` for the
    #' rationale and empirical evidence.
    #' @param trial_niter number of EM iterations used for the cheap trial
    #' candidates (passed to `candidates_split()`/`candidates_merge()`)
    #' before fully re-optimizing only the best one.
    #' @param max_candidates passed to `candidates_merge()` (ignored for
    #' `"split"`, which is never combinatorial in q) -- see its
    #' documentation.
    #' @param directions which neighbor(s) to seed refinement candidates
    #' from: `"split"` (smaller-q neighbor), `"merge"` (larger-q neighbor),
    #' or both (the default).
    #' @param verbose whether to print, for each q attempted, whether the
    #' candidate from that neighbor improved on it. Defaults to
    #' `control$verbose` (the value set at construction, see
    #' [NB_control()]).
    #' @return invisibly returns `self`; improved models replace the
    #' originals in `$models` in place.
    refine = function(trial_niter = 2, max_candidates = 30,
                      directions = c("split", "merge"), verbose = self$control$verbose) {
      stopifnot(all(directions %in% c("split", "merge")))
      ord <- order(self$q_list)
      q_sorted <- self$q_list[ord]
      n <- length(ord)

      try_candidate <- function(idx_target, candidates, label) {
        if (length(candidates) == 0) return(invisible())
        best_candidate <- candidates[[which.min(map_dbl(candidates, "ICL"))]]
        best_candidate$optimize()
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
          candidates <- self$models[[ord[i - 1]]]$candidates_split(trial_niter = trial_niter)
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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field q_list number of blocks
    q_list = function() map_dbl(self$models, "q")
  )
)
