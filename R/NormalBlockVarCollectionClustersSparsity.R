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
  inherit   = NormalBlockCollectionClustersSparsity,
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
      sbm_path <- clustering_path_for_family(mydata, q_list, "var", zero_inflation, control)

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
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field who_am_I a method to print what model is being fitted
    who_am_I  = function(){
      paste("Collection of ",
            ifelse(self$control$zero_inflation, " zero-inflated", ""),
            self$control$noise_covariance,
            "normal-block-var models with different values of q and different penalties.")
    }

  )
)
