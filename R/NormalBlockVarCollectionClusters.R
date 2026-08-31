## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionClusters #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block Models over a Range of Cluster Counts
#'
#' R6 class for a collection of normal-block models with different number of
#' clusters (q) and a fixed sparsity level.
#' @export
NormalBlockVarCollectionClusters <- R6::R6Class(
  classname = "NormalBlockVarCollectionClusters",
  inherit   = NormalBlockCollectionClusters,

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @description Create a new [`NormalBlockVarCollectionClusters`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param q_list list of q values (number of groups) in the collection
    #' @param zero_inflation whether the models in the collection should be zero-inflated or not
    #' @param sparsity sparsity penalty on the network density
    #' @param control structured list of more specific parameters, to generate with NB_control
    #' @return A new [`NormalBlockVarCollectionClusters`] object
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

      # instantiates an NormalBlockVarUnknownClusters model for each q in nb_blocks
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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){paste0(self$control$noise_covariance, " normal-block-var model with unknown q")}
  )
)
