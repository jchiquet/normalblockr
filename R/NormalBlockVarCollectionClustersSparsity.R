## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockVarCollectionClustersSparsity ###############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Collection of Normal-Block-Var Models over Cluster Counts and Sparsity Levels
#'
#' R6 class for a collection of normal-block models with different number of
#' clusters (q) and different sparsity levels, i.e. one
#' [NormalBlockVarCollectionSparsity] per q.
#' @examples
#' ex <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:3, sparsity = TRUE,
#'                        control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
#' models$get_best_model("BIC")$q
#' @export
NormalBlockVarCollectionClustersSparsity <- R6::R6Class(
  classname = "NormalBlockVarCollectionClustersSparsity",
  inherit   = NormalBlockCollectionClustersSparsity,
  private   = list(family = "var")
)
