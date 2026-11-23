## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollection #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection of Normal-Block Models
#'
#' Shared scaffolding for the collections explored by [get_model()]/
#' [normal_block()]: a sweep over sparsity penalties ([`NormalBlockVarCollectionSparsity`]),
#' over the number of clusters ([`NormalBlockVarCollectionClusters`]), or over both
#' ([`NormalBlockVarCollectionClustersSparsity`]). Concrete subclasses set
#' `private$progress_field`/`private$progress_label` in their `initialize()`
#' and provide their own `get_best_model()`, delegating the (row of
#' `self$criteria` minimizing a criterion) lookup to `private$best_id()`.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how collections (NormalBlockVarCollectionClusters,
#' # NormalBlockVarCollectionSparsity, NormalBlockVarCollectionClustersSparsity)
#' # are actually created and fitted.
#' @keywords internal
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
    #' @param control optimization parameters (niter and threshold). When
    #' `control$clustering_init` is `"best_of_inits"`, each leaf model is fit
    #' via its own `best_of_inits()` instead of a plain `optimize()`.
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      boi <- uses_best_of_inits(control)
      self$models <- lapply(self$models, function(model) {
        ## a line rewritten in place only makes sense in a console
        if (control$verbose && interactive()) {
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
          flush.console()
        }
        if (boi && !is.null(model$best_of_inits)) {
          model <- model$best_of_inits(control = control)
        } else {
          model$optimize(control)
        }
        model
      })
      invisible(self)
    },

    #' @description User-friendly print method: model type and the range of
    #' q/sparsity explored. See `summary()` for the full criteria table.
    print = function() {
      crit <- self$criteria
      cat("A", self$who_am_I, "\n")
      cat("===========================================================================\n")
      cat(" ", nrow(crit), "model(s) explored\n")
      if (length(unique(crit$q)) > 1)
        cat("    q ranging from", min(crit$q), "to", max(crit$q), "\n")
      if (length(unique(crit$sparsity)) > 1)
        cat("    sparsity ranging from", signif(min(crit$sparsity), 3),
            "to", signif(max(crit$sparsity), 3), "\n")
      cat("===========================================================================\n")
      cat("* Useful fields\n")
      cat("    $models, $criteria\n")
      cat("* Useful methods\n")
      cat("    print(), summary(), plot(), $get_best_model()\n")
    },

    #' @description Summarize the collection: model type, full criteria
    #' table, and the range of q/sparsity explored.
    #' @return An object of class `summary.NormalBlockCollection`,
    #' printed with a dedicated [print.summary.NormalBlockCollection()]
    #' method.
    summary = function() {
      crit <- self$criteria
      res <- list(
        who_am_I       = self$who_am_I,
        criteria       = crit,
        q_range        = if (length(unique(crit$q)) > 1) range(crit$q) else unique(crit$q),
        sparsity_range = if (length(unique(crit$sparsity)) > 1) range(crit$sparsity) else unique(crit$sparsity)
      )
      class(res) <- "summary.NormalBlockCollection"
      res
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    progress_field = NA_character_, # name of the field reported by optimize()'s progress message
    progress_label = NA_character_, # human-readable label for that field
    scope = function() "", # what varies across the collection (see who_am_I)

    ## Row of crit_df minimizing `crit`. `crit_df` can be passed in to reuse
    ## one already computed by the caller instead of rebuilding self$criteria.
    best_id = function(crit, check_inference = TRUE, crit_df = self$criteria) {
      if (check_inference)
        stopifnot("Log-likelihood based criteria do not apply to the heuristic method" =
                    self$models[[1]]$inference_method == "integrated")
      stopifnot(!anyNA(crit_df[[crit]]))
      id <- 1
      if (length(crit_df[[crit]]) > 1) id <- which.min(crit_df[[crit]])
      id
    },

    ## Shared chart for plot(): one line per criterion against x_var, with a
    ## dashed, colour-matched vline at each criterion's own best model.
    ## "deviance" is excluded by default: monotonic along x_var, so it has no
    ## interior best model the way a penalized criterion does.
    plot_criteria_path = function(x_var, criteria, vline_crit = setdiff(criteria, "deviance")) {
      vline_crit <- intersect(criteria, vline_crit)
      crit_df <- self$criteria
      dplot <- data.frame(x         = rep(crit_df[[x_var]], length(criteria)),
                          criterion = rep(criteria, each = nrow(crit_df)),
                          value     = unlist(crit_df[criteria], use.names = FALSE))
      p <- ggplot2::ggplot(dplot, ggplot2::aes(x = .data$x, y = .data$value, group = .data$criterion,
                                               colour = .data$criterion)) +
        ggplot2::geom_line() + ggplot2::geom_point() +
        ggplot2::ggtitle(label = "Model selection criteria", subtitle = "Lower is better") +
        ggplot2::labs(x = x_var) +
        ggplot2::theme_bw()
      if (length(vline_crit) > 0) {
        ## the model selected by each criterion: a dashed line, a large point
        ## on the criterion's own curve, and a label naming the criterion and
        ## the selected value, stacked from the top so that criteria
        ## selecting the same model stay readable
        best <- vapply(vline_crit, function(crit) private$best_id(crit, crit_df = crit_df), integer(1))
        xs   <- crit_df[[x_var]][best]
        what <- if (x_var == "q") paste("q =", xs) else paste("penalty =", signif(xs, 2))
        rng  <- range(dplot$value, na.rm = TRUE)
        ## labels on the side of the line with more room
        pos  <- (match(xs, sort(unique(crit_df[[x_var]]))) - 1) / max(1, length(unique(crit_df[[x_var]])) - 1)
        dbest <- data.frame(
          criterion = vline_crit,
          x         = xs,
          value     = mapply(function(crit, i) crit_df[[crit]][[i]], vline_crit, best),
          label     = paste0(vline_crit, ": ", what),
          y_label   = rng[2] - (seq_along(vline_crit) - 1) * 0.08 * diff(rng),
          hjust     = ifelse(pos > 0.5, 1, 0)
        )
        p <- p +
          ggplot2::geom_vline(data = dbest, ggplot2::aes(xintercept = .data$x, colour = .data$criterion),
                              linetype = "dashed", alpha = 0.6, show.legend = FALSE) +
          ggplot2::geom_point(data = dbest, ggplot2::aes(x = .data$x, y = .data$value, colour = .data$criterion),
                              size = 3.5, shape = 21, fill = "white", stroke = 1.2, inherit.aes = FALSE,
                              show.legend = FALSE) +
          ggplot2::geom_label(data = dbest, ggplot2::aes(x = .data$x, y = .data$y_label, label = .data$label,
                                                         colour = .data$criterion, hjust = .data$hjust),
                              size = 3, fontface = "bold", inherit.aes = FALSE, show.legend = FALSE) +
          ggplot2::labs(caption = "Dashed lines and circled points: model selected by each criterion")
      }
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a description of the collection: the kind of models it
    #' holds and what varies across them
    who_am_I = function() {
      m <- self$models[[1]]
      while (inherits(m, "NormalBlockCollection")) m <- m$models[[1]]
      paste("collection of", sub(" model with .*$", "", m$who_am_I), "models", private$scope())
    },
    #' @field criteria a data frame with the values of some criteria for the collection of models
    criteria = function() do.call(rbind, map(self$models, "criteria")),
    #' @field loglik not defined for a collection: accessing it raises an informative error. Use `logLik()` for every model's log-likelihood, or `$get_best_model()$loglik` for a single one.
    loglik = function() {
      stop(
        "`$loglik` is not defined for a collection of models (which one?) -- ",
        "use `logLik(<collection>)` for every model's log-likelihood, or ",
        "`<collection>$get_best_model()$loglik` for a single model.",
        call. = FALSE
      )
    }
  )
)
