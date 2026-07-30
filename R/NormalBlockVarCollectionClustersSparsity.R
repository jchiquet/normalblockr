## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionClustersSparsity ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block Models over Cluster Counts and Sparsity Levels
#'
#' R6 class for a collection of normal-block models with different number of
#' clusters (q) and different sparsity levels.
#' @export
NormalBlockVarCollectionClustersSparsity <- R6::R6Class(
  classname = "NormalBlockVarCollectionClustersSparsity",
  inherit   = NormalBlockVarCollection,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarCollectionClustersSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param zero_inflation boolean to specify whether data is zero-inflated
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockVarCollectionClustersSparsity`] object
    initialize = function(mydata, q_list, zero_inflation = FALSE,
                          control = NB_control()) {

      ## Store user-defined fields
      self$control <- control
      self$control$zero_inflation <- zero_inflation
      control_ <- control

      ## One wide SBM exploration over q_list replaces one independent
      ## exploration per q (see sbm_clustering_path()).
      sbm_path <- sbm_path_for_collection(mydata, q_list, zero_inflation, control)

      self$models <- purrr::map(seq_along(q_list), function(rank) {
        if (!is.null(sbm_path))
          control_$clustering_init <- sbm_path[[as.character(q_list[rank])]]
        else if (is.list(control$clustering_init)) # one explicit clustering per q
          control_$clustering_init <- control$clustering_init[[rank]]
        ## else: a heuristic name (or NULL) applies identically to every q,
        ## already carried over since control_ <- control above
        model <- NormalBlockVarCollectionSparsity$new(mydata, q_list[rank], zero_inflation, control_)
      })
      private$progress_field <- "q"
      private$progress_label <- "number of blocks"
    },

    #' @description returns a collection of models corresponding to given q
    #' or one single model if penalty is also given
    #' @param q number of blocks asked by user.
    #' @param sparsity sparsity penalty asked by user
    #' @return either a NormalBlockVarCollectionSparsity or a NormalBlockVarUnknownClusters object
    get_model = function(q, sparsity = NA) {
      stopifnot("No such model in the collection. Acceptable values can be found via $q" = q %in% self$q_list)
      model <- self$models[[which(self$q_list == q)]]
      if (!is.na(sparsity)) {
        if (!(sparsity %in% model$sparsity)) {
          sparsity <-  model$sparsity[[which.min(abs(model$sparsity - sparsity))]]
          cat(paste0("No model with this penalty in the collection. Returning model with closest penalty: ",
                     sparsity,  " Collection penalty values can be found via $sparsity \n"))
        }
        model <- model$models[[which(model$sparsity == sparsity)]]
      }
      model
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to performed the selection.
    #' Either "BIC", "EBIC" or "ICL". "ICL" is the default criterion
    #' @return a [`NormalBlockVarUnknownClusters`] object
    get_best_model = function(crit = c("ICL", "BIC", "EBIC")) {
      crit <- match.arg(crit)
      id <- private$best_id(crit, check_inference = FALSE)
      best_pen <- self$criteria$sparsity[[id]]
      best_q   <- self$criteria$q[[id]]
      self$get_model(best_q, best_pen)$clone()
    },

    #' @description Display various outputs (goodness-of-fit criteria, robustness, diagnostic) associated with a collection of network fits (a [`Networkfamily`])
    #' @param criterion The criteria to plot in `c("deviance", BIC", "EBIC", "ICL")`. Defaults deviance.
    #' @param n_intervals number of intervals into which the penalties range should be split
    #' @return a [`ggplot`] heatmap
    plot = function(criterion = c("deviance", "ICL", "BIC", "EBIC"),
                    n_intervals = NULL) {
      criterion   <- match.arg(criterion)
      if(is.null(n_intervals)) n_intervals <- round(0.1 * length(unique(self$criteria$sparsity )))
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
    #' @field q_list  number of blocks
    q_list = function() map_dbl(self$models, "q"),
    #' @field sparsity list of penalties used for each q
    sparsity = function(){
      self$criteria %>%
        dplyr::group_by(q) %>%
        dplyr::summarize(sparsity = paste(round(sparsity, 2), collapse = ", "))
    },
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){
      paste("Collection of ",
            ifelse(self$control$zero_inflation, " zero-inflated", ""),
            self$control$noise_covariance,
            "normal-block models with different values of q and different penalties.")
    }

  )
)
