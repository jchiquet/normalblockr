## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionClusters ############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Normal-Block-Mean Models over a Range of Cluster Counts
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase])
#' with different numbers of clusters (q). As for
#' [NormalBlockVarCollectionClusters], the q-independent part of the initial
#' clustering heuristic is computed once for the whole collection.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 2:5, model = "mean")
#' plot(models, criteria = c("BIC", "ICL"))
#' models$get_best_model()
#' @export
NormalBlockMeanCollectionClusters <- R6::R6Class(
  classname = "NormalBlockMeanCollectionClusters",
  inherit   = NormalBlockCollectionClusters,
  private   = list(family = "mean")
)
