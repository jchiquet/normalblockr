## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionClusters #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Mean-Block Models over a Range of Cluster Counts
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase])
#' with different numbers of clusters (q). Inherits its scaffolding
#' (`print()`/`summary()`/`plot()`/`optimize()`, the `criteria` table) from
#' [NormalBlockCollection], which despite its name is generic across
#' both model families -- unlike [NormalBlockVarCollectionClusters], there
#' is no SBM-path shortcut here: the shared clustering-heuristic registry's
#' cov()/correlation-based methods are ill-suited to the mean-block family
#' (see [NormalBlockMeanBase]'s own default), so each q is fit
#' independently.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:5, model = "mean")
#' models$plot(c("BIC", "ICL"))
#' models$get_best_model()
#' @export
NormalBlockMeanCollectionClusters <- R6::R6Class(
  classname = "NormalBlockMeanCollectionClusters",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockMeanCollectionClusters`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockMeanCollectionClusters`] object
    initialize = function(mydata, q_list, sparsity = 0, control = NB_control()) {
      stopifnot("each nb_blocks value can only be present once in nb_blocks" =
                  length(q_list) == length(unique(q_list)))
      stopifnot("There cannot be more blocks than there are entities to cluster." =
                  max(q_list) <= ncol(mydata$Y))

      self$control <- control

      this_control <- control
      self$models <- map(rank(q_list),
          function(r) {
            ## a list means one explicit clustering per q; anything else
            ## (a heuristic name, or NULL) applies identically to every q
            this_control$clustering_init <-
              if (is.list(control$clustering_init)) control$clustering_init[[r]]
              else control$clustering_init
            get_model(mydata, q_list[r], sparsity = sparsity,
                     model = "mean", control = this_control)
        })
      private$progress_field <- "q"
      private$progress_label <- "number of blocks"
    },

    #' @description returns the [NormalBlockMeanUnknownClusters] model corresponding to given q
    #' @param q number of blocks asked by user
    #' @return A [NormalBlockMeanUnknownClusters] object with given value q
    get_model = function(q) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      self$models[[which(self$q_list == q)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "ICL" or "BIC". "ICL" is the default criterion
    #' @return a [`NormalBlockMeanUnknownClusters`] object
    get_best_model = function(crit = c("ICL", "BIC", "EBIC", "deviance")) {
      crit <- match.arg(crit)
      id <- private$best_id(crit)
      self$models[[id]]$clone()
    },

    #' @description Display various goodness-of-fit criteria against q
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
    #' are refined. Same mechanism as
    #' [NormalBlockVarCollectionClusters]'s `refine()`, backed by
    #' [NormalBlockMeanBase]'s `candidates_split()`/`candidates_merge()`.
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
    who_am_I = function() "normal-block-mean model with unknown q"
  )
)
