## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockCollectionClustersSparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Base Class for a Collection over Cluster Counts and Sparsity Levels
#'
#' Shared scaffolding for [NormalBlockVarCollectionClustersSparsity] and
#' [NormalBlockMeanCollectionClustersSparsity]: a collection of sparsity
#' sub-collections, one per q. Everything family-agnostic (two-key model
#' lookup, model selection over both axes, the criteria heatmap) lives here;
#' subclasses only build `self$models` and name themselves.
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

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits 
    #' @param criterion The criteria to plot in `c("deviance", BIC", "EBIC", "ICL")`. Defaults deviance.
    #' @param n_intervals number of intervals into which the penalties range should be split
    #' @return a [`ggplot2::ggplot`] heatmap
    plot = function(criterion = c("deviance", "ICL", "BIC", "EBIC"),
                    n_intervals = NULL) {
      criterion   <- match.arg(criterion)
      ## cut() needs at least two intervals; a short penalty path would
      ## otherwise round down to zero and make it fail
      if (is.null(n_intervals))
        n_intervals <- max(2, round(0.1 * length(unique(self$criteria$sparsity))))
      df <- self$criteria %>% dplyr::mutate(pen_binned = cut(sparsity, breaks = n_intervals)) %>%
        dplyr::group_by(pen_binned, q) %>% dplyr::summarize(avg_crit = mean(.data[[criterion]]), .groups = "drop")
      p  <- ggplot2::ggplot(df, ggplot2::aes(x = pen_binned, y = q, fill = avg_crit)) +
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
  active = list(
    #' @field q_list number of blocks
    q_list = function() map_dbl(self$models, "q"),
    #' @field sparsity list of penalties used for each q
    sparsity = function(){
      self$criteria %>%
        dplyr::group_by(q) %>%
        dplyr::summarize(sparsity = paste(round(sparsity, 2), collapse = ", "))
    }
  )
)
