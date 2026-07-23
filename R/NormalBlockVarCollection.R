## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollection #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a collection of normal-block models
#'
#' Shared scaffolding for the collections explored by [get_model()]/
#' [normal_block()]: a sweep over sparsity penalties ([`NormalBlockVarCollectionSparsity`]),
#' over the number of clusters ([`NormalBlockVarCollectionClusters`]), or over both
#' ([`NormalBlockVarCollectionClustersSparsity`]). Concrete subclasses set
#' `private$progress_field`/`private$progress_label` in their `initialize()`
#' and provide their own `get_best_model()`, delegating the (row of
#' `self$criteria` minimizing a criterion) lookup to `private$best_id()`.
NormalBlockVarCollection <- R6::R6Class(
  classname = "NormalBlockVarCollection",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field models list of models (or sub-collections) explored by the collection
    models  = NULL,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,

    #' @description optimizes every model (or sub-collection) in the collection
    #' @param control optimization parameters (niter and threshold). When
    #' `control$clustering_init` is `"best_of_inits"`, each leaf model is fit
    #' via its own `best_of_inits()` instead of a plain `optimize()`.
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      boi <- uses_best_of_inits(control)
      self$models <- lapply(self$models, function(model) {
        if (control$verbose)
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
        flush.console()
        if (boi && !is.null(model$best_of_inits)) {
          model <- model$best_of_inits(control = control)
        } else {
          model$optimize(control)
        }
        model
      })
      invisible(self)
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    progress_field = NA_character_, # name of the field reported by optimize()'s progress message
    progress_label = NA_character_, # human-readable label for that field

    ## Row of crit_df minimizing `crit`, after checking that the criterion is
    ## well-defined for the whole collection. Centralizes the `length(...) >
    ## 1` guard so that the id is always defined, even for a collection
    ## reduced to a single model. `crit_df` defaults to (a fresh) self$criteria
    ## but can be passed in explicitly to reuse one already computed by the
    ## caller (see plot_criteria_path() below) instead of rebuilding it --
    ## self$criteria rebuilds the whole collection's criteria data frame from
    ## every model's own `criteria` field on every access, not a free lookup.
    best_id = function(crit, check_inference = TRUE, crit_df = self$criteria) {
      if (check_inference)
        stopifnot("Log-likelihood based criteria do not apply to the heuristic method" =
                    self$models[[1]]$inference_method == "integrated")
      stopifnot(!anyNA(crit_df[[crit]]))
      id <- 1
      if (length(crit_df[[crit]]) > 1) id <- which.min(crit_df[[crit]])
      id
    },

    ## Shared chart for plot() in NormalBlockVarCollectionClusters (x_var = "q") and
    ## NormalBlockVarCollectionSparsity (x_var = "sparsity"): one line per
    ## criterion against x_var, with a dashed vline at each criterion's own
    ## best model, colour-matched to that criterion's line (so which vline
    ## belongs to which criterion is unambiguous without a second legend --
    ## a single grey vline used to mark only one hardcoded criterion's best
    ## model regardless of how many were actually plotted, which made it
    ## impossible to tell, e.g., the BIC-best from the ICL-best at a glance).
    ## "deviance" is excluded by default: it is monotonic along x_var (more
    ## blocks/less penalty always fits at least as well), so it has no
    ## interior best model the way a penalized criterion does.
    ##
    ## self$criteria is computed once, locally (crit_df), and reused for both
    ## the line plot and every vline's position: vline x-positions are read
    ## directly off crit_df via best_id(), not via get_best_model(), which
    ## would (a) rebuild self$criteria all over again and (b) clone() a whole
    ## fitted model just to read a single scalar (q or sparsity) off it.
    plot_criteria_path = function(x_var, criteria, vline_crit = setdiff(criteria, "deviance")) {
      vline_crit <- intersect(criteria, vline_crit)
      crit_df <- self$criteria
      dplot <- crit_df %>%
        dplyr::select(dplyr::all_of(c(x_var, criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -dplyr::all_of(x_var)) %>%
        dplyr::group_by(criterion)
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = .data[[x_var]], y = value, group = criterion, colour = criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label = "Model selection criteria", subtitle = "Lower is better") +
        ggplot2::theme_bw()
      if (length(vline_crit) > 0) {
        dvlines <- tibble::tibble(
          criterion = vline_crit,
          x         = sapply(vline_crit, function(crit) crit_df[[x_var]][[private$best_id(crit, crit_df = crit_df)]])
        )
        p <- p + ggplot2::geom_vline(data = dvlines, ggplot2::aes(xintercept = x, colour = criterion),
                                      linetype = "dashed", alpha = 0.6, show.legend = FALSE)
      }
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field criteria a data frame with the values of some criteria for the collection of models
    criteria = function() purrr::map_df(self$models, "criteria")
  )
)
