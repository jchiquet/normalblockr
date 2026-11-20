## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionClusters #################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block-Var Models over a Range of Cluster Counts
#'
#' R6 class for a collection of normal-block models with different number of
#' clusters (q) and a fixed sparsity level. The q-independent part of the
#' initial clustering heuristic is computed once for the whole collection.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:4, control = NB_control(verbose = FALSE))
#' models$get_best_model("ICL")$q
#' @export
NormalBlockVarCollectionClusters <- R6::R6Class(
  classname = "NormalBlockVarCollectionClusters",
  inherit   = NormalBlockCollectionClusters,
  private   = list(family = "var")
)
