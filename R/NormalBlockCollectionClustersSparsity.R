## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionClustersSparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection over Cluster Counts and Sparsity Levels
#'
#' Shared scaffolding for [NormalBlockVarCollectionClustersSparsity] and
#' [NormalBlockMeanCollectionClustersSparsity]: a collection of sparsity
#' sub-collections, one per q in increasing order. Everything lives here
#' (two-key model lookup, model selection over both axes, the criteria
#' heatmap); subclasses only name the model family.
#' @examples
#' # An internal abstract base class, never instantiated directly -- see
#' # normal_block() for how collections are created and fitted.
#' @keywords internal
NormalBlockCollectionClustersSparsity <- R6::R6Class(
  classname = "NormalBlockCollectionClustersSparsity",
  inherit   = NormalBlockCollection,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a collection over a range of q, sorted
    #' increasingly, with a sparsity path for each q.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list values of q (number of clusters) in the collection
    #' @param zero_inflation whether the models of the collection are zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new collection object
    initialize = function(mydata, q_list, zero_inflation = FALSE, control = NB_control()) {
      check_q_list(q_list, mydata$p)
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      ord <- order(q_list)
      inits <- clustering_inits(mydata, q_list, private$family, zero_inflation, control)[ord]
      path_class <- get(paste0("NormalBlock", if (private$family == "var") "Var" else "Mean",
                               "CollectionSparsity"))
      self$models <- map2(q_list[ord], inits, function(q, init) {
        control$clustering_init <- init
        path_class$new(mydata, q, zero_inflation, control)
      })
      private$progress_field <- "q"
      private$progress_label <- "number of clusters"
    },

    #' @description returns a collection of models corresponding to given q
    #' or one single model if penalty is also given
    #' @param q number of blocks asked by user.
    #' @param sparsity sparsity penalty asked by user
    #' @return either a sparsity sub-collection or a single model object
    get_model = function(q, sparsity = NA) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      model <- self$models[[which(self$q_list == q)]]
      if (!is.na(sparsity)) {
        if (!(sparsity %in% model$sparsity)) {
          sparsity <-  model$sparsity[[which.min(abs(model$sparsity - sparsity))]]
          message("No model with this penalty in the collection. Returning model with closest penalty: ",
                  sparsity,  " Collection penalty values can be found via $sparsity")
        }
        model <- model$models[[which(model$sparsity == sparsity)]]
      }
      model
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "EBIC" or "ICL". "ICL" is the default criterion
    #' @return a single fitted model object
    get_best_model = function(crit = c("ICL", "BIC", "EBIC")) {
      crit <- match.arg(crit)
      id <- private$best_id(crit, check_inference = FALSE)
      best_pen <- self$criteria$sparsity[[id]]
      best_q   <- self$criteria$q[[id]]
      self$get_model(best_q, best_pen)$clone()
    },

    #' @description Plot one criterion as a heatmap over q and the penalty
    #' (see also [plot.NormalBlockCollection()]).
    #' @param criterion the criterion to plot, among `"deviance"` (the
    #' default), `"ICL"`, `"BIC"` and `"EBIC"`.
    #' @param n_intervals number of intervals into which the range of the
    #' penalties is split
    #' @return a [`ggplot2::ggplot`] heatmap
    plot = function(criterion = c("deviance", "ICL", "BIC", "EBIC"),
                    n_intervals = NULL) {
      criterion   <- match.arg(criterion)
      ## cut() needs at least two intervals; a short penalty path would
      ## otherwise round down to zero and make it fail
      if (is.null(n_intervals))
        n_intervals <- max(2, round(0.1 * length(unique(self$criteria$sparsity))))
      crit <- self$criteria
      crit$pen_binned <- cut(crit$sparsity, breaks = n_intervals)
      df <- stats::aggregate(crit[criterion], by = crit[c("pen_binned", "q")], FUN = mean)
      names(df)[3] <- "avg_crit"
      p  <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pen_binned, y = .data$q, fill = .data$avg_crit)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(label    = criterion,
                         subtitle = "Lower is better" ) +
        ggplot2::labs(x = "Penalties (Binned)", y = "q", fill = paste0("Average ", criterion),
                      title = criterion) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    family = NA_character_, # "var" or "mean", set by the subclasses
    scope = function() "over q and sparsity penalties"
  ),
  active = list(
    #' @field q_list number of blocks
    q_list = function() map_dbl(self$models, "q"),
    #' @field sparsity list of penalties used for each q
    sparsity = function(){
      crit <- self$criteria
      stats::aggregate(list(sparsity = crit$sparsity), by = list(q = crit$q),
                       FUN = function(s) paste(round(s, 2), collapse = ", "))
    }
  )
)
