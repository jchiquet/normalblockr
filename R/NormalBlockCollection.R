## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollection #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 abstract class for a collection of normal-block models
#'
#' Shared scaffolding for the collections explored by [get_model()]/
#' [normal_block()]: a sweep over sparsity penalties ([`NormalBlockCollectionSparsity`]),
#' over the number of clusters ([`NormalBlockCollectionClusters`]), or over both
#' ([`NormalBlockCollectionClustersSparsity`]). Concrete subclasses set
#' `private$progress_field`/`private$progress_label` in their `initialize()`
#' and provide their own `get_best_model()`, delegating the (row of
#' `self$criteria` minimizing a criterion) lookup to `private$best_id()`.
NormalBlockCollection <- R6::R6Class(
  classname = "NormalBlockCollection",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field models list of models (or sub-collections) explored by the collection
    models  = NULL,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,

    #' @description optimizes every model (or sub-collection) in the collection
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 100, threshold = 1e-4, verbose = TRUE)) {
      self$models <- lapply(self$models, function(model) {
        if (control$verbose)
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
        flush.console()
        model$optimize(control)
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

    ## Row of self$criteria minimizing `crit`, after checking that the
    ## criterion is well-defined for the whole collection. Centralizes the
    ## `length(...) > 1` guard so that the id is always defined, even for a
    ## collection reduced to a single model.
    best_id = function(crit, check_inference = TRUE) {
      if (check_inference)
        stopifnot("Log-likelihood based criteria do not apply to the heuristic method" =
                    self$models[[1]]$inference_method == "integrated")
      stopifnot(!anyNA(self$criteria[[crit]]))
      id <- 1
      if (length(self$criteria[[crit]]) > 1) id <- which.min(self$criteria[[crit]])
      id
    },

    ## Shared chart for plot() in NormalBlockCollectionClusters (x_var = "q") and
    ## NormalBlockCollectionSparsity (x_var = "sparsity"): one line per
    ## criterion against x_var, with a dashed vline marking the best model
    ## under vline_crit. Was previously duplicated near-verbatim in both.
    plot_criteria_path = function(x_var, criteria, vline_crit) {
      vlines <- sapply(intersect(criteria, vline_crit), function(crit) self$get_best_model(crit)[[x_var]])
      dplot <- self$criteria %>%
        dplyr::select(dplyr::all_of(c(x_var, criteria))) %>%
        tidyr::gather(key = "criterion", value = "value", -dplyr::all_of(x_var)) %>%
        dplyr::group_by(criterion)
      ggplot2::ggplot(dplot, ggplot2::aes(x = .data[[x_var]], y = value, group = criterion, colour = criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label = "Model selection criteria", subtitle = "Lower is better") +
        ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
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
