## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS SelectionNClusters ############################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for selecting the number of clusters (q) by forward/backward
#' split-and-merge exploration, comparing models with the ICL
#'
#' Kept internal/unexported for reference: empirically (see
#' [NormalBlockCollectionClusters]'s `refine()` documentation) chaining every
#' q from a single starting point this way is no longer recommended over
#' [NormalBlockCollectionClusters]'s independent per-q cold starts followed
#' by `refine()` -- the latter matches or beats this class's quality at a
#' fraction of the cost on every real dataset tested, and has no equivalent
#' to this class's failure mode (an early bad split/merge propagating to
#' every larger q it explores from there).
#' @import tibble
#' @keywords internal
SelectionNClusters <- R6::R6Class(
  classname = "SelectionNClusters",

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field best_models list of models explored so far indexed by the cluster sizes
    best_models = NULL,
    #' @field ICL_explored list of ICL reached so far, indexed by the cluster sizes
    ICL_explored = NULL,
    #' @field n_clusters_range the range of the cluster sizes to browse
    n_clusters_range = NULL,
    #' @field control store the list of user-defined model settings and optimization parameters
    control = NA,

    #' @description Explores and selects the optimal number of classes
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param n_clusters_range a size-2 vector with the range of cluster size to browse
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return The best model in the ICL sense
    initialize = function(mydata, n_clusters_range, zero_inflation = FALSE,
                          sparsity = 0, control = NB_control()) {

      self$control <- control
      self$control$zero_inflation <- zero_inflation

      n_clusters_range <- sort(n_clusters_range)
      stopifnot("bad formatting of n_clusters_range" =
                  all(n_clusters_range > 0, length(n_clusters_range) == 2,
                      n_clusters_range[2] <= ncol(mydata$Y)))
      self$n_clusters_range <- n_clusters_range

      # instantiates the first model
      model0 <- get_model(mydata,
                          n_clusters_range[1],
                          sparsity = sparsity,
                          zero_inflation = zero_inflation,
                          control = control)
      model0$optimize()

      self$best_models <- tibble(
        n_clusters     = model0$q,
        split_explored = FALSE,
        merge_explored = TRUE ,
        ICL            = model0$ICL,
        model          = list(model0)
      )
      self$ICL_explored <- data.frame(n_clusters = model0$q, ICL = model0$ICL)
    },

    #' @description perform model selection with forward/backward exploration
    #' with split and merge strategy
    fit = function() {
      self$explore_forward()
      self$explore_backward()
      invisible(self)
    },

    #' @description perform model selection with forward/backward exploration
    #' with split and merge strategy
    #' @param model a normal-block model
    #' @param strategy a character, either "split" or "merge"
    #' @param max_training maximal of model fully trained at each step
    #' of the exploration
    train_best_candidates = function(model, strategy, max_training = 3) {
      if (strategy == "merge")
        candidates <- model$candidates_merge()
      else
        candidates <- model$candidates_split()

      ibest <- candidates %>% map_dbl("ICL") %>% order()
      best_candidates <- candidates[na.omit(ibest[1:max_training])]
      map(best_candidates, function(model) model$optimize())
      invisible(best_candidates)
    },

    #' @description perform forward exploration with a split strategy
    explore_forward = function() {

      while (1) {
        ## Select the model from which to split
        icurrent <- with(self$best_models, which(split_explored == FALSE))

        ## Exit if all model explored
        if (length(icurrent) == 0) break

        ## Exit if maximum nb of clusters is reached
        current <- self$best_models$model[[icurrent[1]]]
        self$best_models$split_explored[icurrent[1]] <- TRUE
        if (current$q >= self$n_clusters_range[2]) break

        ## Exit if no candidate is found
        candidates <- self$train_best_candidates(current, "split")
        if (length(candidates) == 0) break

        cat("Explore by spliting a model with", current$q, "clusters \r")
        candidates_ICL <- map_dbl(candidates, "ICL")
        best_model <- candidates[[which.min(candidates_ICL)]]

        icurrent_best <- with(self$best_models, which(n_clusters == best_model$q))
        if (length(icurrent_best) == 0) {
          self$best_models <- self$best_models %>% add_row(
            n_clusters     = best_model$q,
            split_explored = FALSE,
            merge_explored = FALSE,
            ICL            = best_model$ICL,
            model = list(best_model)
          )
          self$ICL_explored <- rbind(
            self$ICL_explored,
            data.frame(n_clusters = best_model$q, ICL = candidates_ICL)
          )
        } else if (best_model$ICL < self$best_models[icurrent_best, "ICL"]) {
          self$best_models$ICL[icurrent_best] <- best_model$ICL
          self$best_models$model[[icurrent_best]] <- best_model
        }
      }
      invisible(self)
    },

    #' @description perform backward exploration with a merge strategy
    explore_backward = function() {

      while (1) {
        ## Select the model from which to split
        ocurrent <- order(self$best_models$n_clusters, decreasing = TRUE)
        icurrent <- na.omit(with(self$best_models, which(merge_explored == FALSE))[ocurrent])

        ## Exit if all model explored
        if (length(icurrent) == 0) break
        current <- self$best_models$model[[icurrent[1]]]
        self$best_models$merge_explored[icurrent[1]] <- TRUE
        cat("Explore by merging a model with", current$q, "clusters \r")

        ## Exit if the minimum nb of clusters is reached
        if (current$q <= 2) break

        ## Exit if no candidate is found
        candidates <- self$train_best_candidates(current, "merge")
        no_candidate <- length(candidates) == 0
        if (no_candidate) break

        candidates_ICL <- map_dbl(candidates, "ICL")
        best_model <- candidates[[which.min(candidates_ICL)]]
        icurrent_best <- with(self$best_models, which(n_clusters == best_model$q))
        if (best_model$ICL < self$best_models[icurrent_best, "ICL"]) {
          self$best_models$ICL[icurrent_best] <- best_model$ICL
          self$best_models$model[[icurrent_best]] <- best_model
        }

        self$ICL_explored <- rbind(
          self$ICL_explored,
          data.frame(n_clusters = best_model$q, ICL = candidates_ICL)
        )
      }
      invisible(self)
    },


    #' @description Display the ICL for all the best models explored per number of cluster, and the winner
    #' @return a [`ggplot2::ggplot`] graph
    plot = function() {
      data_min <- self$ICL_explored %>% dplyr::group_by(n_clusters) %>% dplyr::summarise(min_ICL = min(ICL))
      vlines <- data_min$n_clusters[which.min(data_min$min_ICL)]
      p <- ggplot2::ggplot(self$ICL_explored) +
        ggplot2::geom_point(ggplot2::aes(x = n_clusters, y = ICL), alpha = .75) +
        ggplot2::geom_line(data = data_min, ggplot2::aes(x = n_clusters, y = min_ICL), color = "red") +
        ggplot2::ggtitle(label    = "Model selection criteria",
                         subtitle = "Lower is better" ) + ggplot2::xlab("number of clusters") +
        ggplot2::scale_x_continuous(breaks = scales::pretty_breaks()) +
        ggplot2::theme_bw() + ggplot2::geom_vline(xintercept = vlines, linetype = "dashed", alpha = 0.25)
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field best_model best model explored so far in term of ICL
    best_model = function() self$best_models$model[[which.min(self$best_models$ICL)]]
  )
)
