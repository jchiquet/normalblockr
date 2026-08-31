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

#' @title Print a Collection of Normal-Block Models
#' @description Print a short summary of a collection of normal-block
#' models: model type and the range of q/sparsity explored. See
#' [summary.NormalBlockVarCollection()] for the full criteria table.
#' @param x An object inheriting from NormalBlockVarCollection.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
#' print(models)
print.NormalBlockVarCollection <- function(x, ...) {
  stopifnot(isNBcollection(x))
  x$print()
  invisible(x)
}

#' @title Summarize a Collection of Normal-Block Models
#' @description Summarizes a collection of normal-block models: model type,
#' the full criteria table, and the range of q/sparsity explored.
#' @param object An object inheriting from NormalBlockVarCollection.
#' @param ... not used, only here for S3 compatibility
#' @return An object of class `summary.NormalBlockVarCollection` (a list
#' with the collection's `who_am_I`, `criteria`, `q_range` and
#' `sparsity_range`), printed with a dedicated
#' [print.summary.NormalBlockVarCollection()] method.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
#' summary(models)
summary.NormalBlockVarCollection <- function(object, ...) {
  stopifnot(isNBcollection(object))
  object$summary()
}

#' @title Print a Collection Summary
#' @description Print method for objects returned by [summary.NormalBlockVarCollection()].
#' @param x An object of class `summary.NormalBlockVarCollection`.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
print.summary.NormalBlockVarCollection <- function(x, ...) {
  cat("A", x$who_am_I, "\n")
  cat("===========================================================================\n")
  if (length(x$q_range) > 1)
    cat("    q ranging from", x$q_range[1], "to", x$q_range[2], "\n")
  if (length(x$sparsity_range) > 1)
    cat("    sparsity ranging from", signif(x$sparsity_range[1], 3),
        "to", signif(x$sparsity_range[2], 3), "\n")
  cat("===========================================================================\n")
  print(as.data.frame(round(x$criteria, digits = 3), row.names = rep("", nrow(x$criteria))))
  invisible(x)
}
