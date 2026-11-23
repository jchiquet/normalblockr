## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionClustersSparsity ##############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Normal-Block-Mean Models over Cluster Counts and Sparsity Levels
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase]) over
#' both a range of cluster counts (q) and a sparsity path, i.e. one
#' [NormalBlockMeanCollectionSparsity] per q. Mirrors
#' [NormalBlockVarCollectionClustersSparsity].
#' @examples
#' ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:4, sparsity = TRUE, model = "mean",
#'                        control = NB_control(n_sparsity_penalties = 4))
#' plot(models, criterion = "BIC")
#' models$get_best_model("BIC")
#' @export
NormalBlockMeanCollectionClustersSparsity <- R6::R6Class(
  classname = "NormalBlockMeanCollectionClustersSparsity",
  inherit   = NormalBlockCollectionClustersSparsity,
  private   = list(family = "mean")
)
