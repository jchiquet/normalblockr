## =========================================================================================
##
## PUBLIC S3 METHODS FOR NormalBlockVarCollection
##
## =========================================================================================

isNBcollection <- function(object) inherits(object, "NormalBlockVarCollection")

#' @title Extract Log-Likelihood of a Collection of Normal-Block Models
#' @description Returns the log-likelihood of every model in a collection of
#' normal-block models (see [NormalBlockVarCollectionClusters],
#' [NormalBlockVarCollectionSparsity], [NormalBlockVarCollectionClustersSparsity]).
#' @param object An object inheriting from NormalBlockVarCollection.
#' @param ... not used, only here for S3 compatibility
#' @return A numeric vector of log-likelihood values, one per model in the collection.
#' @importFrom stats logLik
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
#' logLik(models)
logLik.NormalBlockVarCollection <- function(object, ...) {
  stopifnot(isNBcollection(object))
  object$criteria$loglik
}

#' @title Bayesian Information Criterion for a Collection of Normal-Block Models
#' @description Returns the (variational) BIC of every model in a collection
#' of normal-block models.
#' @param object An object inheriting from NormalBlockVarCollection.
#' @param ... not used, only here for S3 compatibility
#' @return A numeric vector of BIC values, one per model in the collection.
#' @importFrom stats BIC
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
#' BIC(models)
BIC.NormalBlockVarCollection <- function(object, ...) {
  stopifnot(isNBcollection(object))
  object$criteria$BIC
}
