## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionClusters #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' R6 class for a collection of normal-block models with different number of clusters (q) and a fixed sparsity level.
#' @export
NormalBlockCollectionClusters <- R6::R6Class(
  classname = "NormalBlockCollectionClusters",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockCollectionClusters`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockCollectionClusters`] object
    initialize = function(mydata, q_list, zero_inflation = FALSE,
                          sparsity = 0, control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(q_list) == length(unique(q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(q_list) <= ncol(mydata$Y))

      self$control <- control
      self$control$zero_inflation <- zero_inflation

      ## A single wide SBM exploration over the whole q_list range visits
      ## every intermediate block count on its way there, so it gives a
      ## clustering for every q in q_list at a fraction of the cost of letting
      ## each model run its own independent exploration (clustering_init =
      ## "sbm" otherwise repeats the explore step once per q -- see
      ## sbm_path_for_collection()/sbm_clustering_path() in R/utils.R).
      sbm_path <- sbm_path_for_collection(mydata, q_list, zero_inflation, control)

      # instantiates an NormalBlockUnknownClusters model for each q in nb_blocks
      this_control <- control
      self$models <- map(rank(q_list),
          function(r) {
            this_control$clustering_init <-
              if (!is.null(sbm_path)) sbm_path[[as.character(q_list[r])]]
              ## a list means one explicit clustering per q; anything else
              ## (a heuristic name, or NULL) applies identically to every q
              else if (is.list(control$clustering_init)) control$clustering_init[[r]]
              else control$clustering_init
            model <- get_model(mydata,
                               q_list[r],
                               sparsity = sparsity,
                               zero_inflation = zero_inflation,
                               control = this_control)
        })
      private$progress_field <- "q"
      private$progress_label <- "number of blocks"
    },

    #' @description returns the NormalBlockUnknownClusters model corresponding to given q
    #' @param q number of blocks asked by user
    #' @return A NormalBlockUnknownClusters object with given value q
    get_model = function(q) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      self$models[[which(self$q_list == q)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "ICL" or "BIC". "ICL" is the default criterion
    #' @return a [`NormalBlockUnknownClusters`] object
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
      private$plot_criteria_path("q", criteria, "ICL")
    },

    #' @description optimizes every model in the collection, then -- if
    #' `control$refine` is `TRUE` (see [NB_control()], default `FALSE`) --
    #' calls [refine()] automatically.
    #' @param control optimization parameters (niter, threshold, verbose)
    optimize = function(control = list(niter = 100, threshold = 1e-4, verbose = TRUE)) {
      super$optimize(control)
      if (isTRUE(self$control$refine)) self$refine()
      invisible(self)
    },

    #' @description Tries to improve every model in the collection with a
    #' short split-and-reoptimize trial seeded from its smaller-q neighbor
    #' (`"split"`), a short merge-and-reoptimize trial seeded from its
    #' larger-q neighbor (`"merge"`), or both (the default); a candidate
    #' replaces the original only if it strictly lowers the deviance. Useful
    #' because every model in the collection is cold-started independently
    #' (a heuristic clustering on the residuals, see
    #' `private$clustering_methods` in [NormalBlockBase]) and can settle
    #' into a milder local optimum than an adjacent q's solution would, once
    #' split/merged by one cluster -- a gap that "deviance must not increase
    #' in q" alone cannot reliably catch (see Details).
    #'
    #' @details Two cheaper attempts at targeting *which* q to refine --
    #' flagging outright deviance increases, then flagging drops unusually
    #' small relative to neighboring q's -- found nothing to fix on a real
    #' dataset (`brca_rppa`) where a full split/merge search was known to do
    #' systematically better at almost every q: the gap was spread evenly
    #' across the whole range, leaving no local anomaly to detect.
    #' `refine()` instead tries unconditionally everywhere and discards
    #' whatever doesn't help. On `brca_rppa`, `"split"` alone already
    #' recovers part of that gap (~1.5-2x the cost of fitting the collection
    #' alone); adding `"merge"` closes substantially more of it -- the two
    #' directions catch different cases, since split-from-below and
    #' merge-from-above are different, independently cold-started starting
    #' points -- for about twice the cost again, still well under
    #' [SelectionNClusters]'s full forward/backward search. Because nothing
    #' is kept unless strictly better, this never makes a model worse: on a
    #' second dataset (`university`) where independent cold starts were
    #' already better than chaining models together, both directions found
    #' close to nothing to fix rather than degrading anything.
    #'
    #' Only contiguous q pairs (`q` and `q -/+ 1`, both present in the
    #' collection) are refined; gaps in `q_list` are left untouched on that
    #' side.
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
    q_list = function() map_dbl(self$models, "q"),
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){paste0(self$control$noise_covariance, " normal-block model with unknown q")}
  )
)
