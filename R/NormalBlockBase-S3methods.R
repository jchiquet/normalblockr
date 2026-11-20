## ===================================================================================
##
## PUBLIC S3 METHODS FOR NormalBlockBase (all families)
##
## ===================================================================================

#' @title Check if an Object is a Normal-Block Model
#' @description Checks if a model is of class [NormalBlockBase()], i.e. any
#' normal-block model (variance-block or mean-block family).
#' @param object An R object.
#' @return A boolean telling whether object inherits from the NormalBlockBase class.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' isNB(model)
#' @export
isNB <- function(object) {inherits(object, "NormalBlockBase")}

#' @title Extract Model Coefficients
#' @description Extract the regression coefficients B of a normal-block model.
#' For the variance-block models, one column per variable, on the original
#' scale of Y even when the variables were rescaled for the fit. For the
#' mean-block models, one column per cluster, on the scale of the fit (the
#' original one by default), since a profile is shared by variables of
#' different scales.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A d x p (variance-block) or d x q (mean-block) matrix.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' coef(model)[, 1:5]
#' @export
coef.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  B <- object$model_par$B
  if (inherits(object, "NormalBlockVarBase")) {
    B <- B * matrix(object$data$Y_scale, nrow(B), ncol(B), byrow = TRUE)
    cols <- colnames(object$data$Y)
  } else {
    cols <- cluster_names(ncol(B))
  }
  with_dimnames(B, colnames(object$data$X), cols)
}

#' @title Extract the Covariance Matrix
#' @description Extract the covariance matrix `Omega^-1`. For the
#' variance-block models, the q x q covariance of the latent block signals,
#' on the scale of the fit (it has no original-scale counterpart when the
#' variables of a block were rescaled by different factors). For the
#' mean-block models, the p x p residual covariance between variables, on
#' the original scale of Y.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return The covariance matrix, of size q x q or p x p depending on the model.
#' @importFrom stats sigma
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' sigma(model)
#' @export
sigma.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  S <- chol2inv(chol(object$model_par$Omega))
  if (inherits(object, "NormalBlockMeanBase")) S <- S * tcrossprod(object$data$Y_scale)
  with_dimnames(S, network_nodes_names(object), network_nodes_names(object))
}

#' @title Extract Fitted Values
#' @description Extract the fitted values of a normal-block model, on the
#' original scale of Y. For the variance-block models, they include the
#' posterior mean of the latent block signals, `X B + E[W] C'`; for the
#' mean-block models, they are `X B C'` (with the posterior cluster
#' probabilities when the clustering is unknown). Zero-inflated models set
#' the observed zeros to 0.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A matrix of fitted values extracted from the object.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' head(fitted(model)[, 1:5])
#' @export
fitted.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  with_dimnames(object$fitted, rownames(object$data$Y), colnames(object$data$Y))
}


#' @title Predict Method for Normal-Block Models
#' @description Predicts the mean of the observations Y for new covariates X,
#' on the original scale of Y: `B' X` for the variance-block models (the
#' latent block signals, centered, average out; unlike [fitted()], which adds
#' their posterior mean for the training observations), `C B' X` for the
#' mean-block models. For zero-inflated models, the mean of the Gaussian
#' component.
#' @param object An object of class NormalBlockBase.
#' @param new_X New set of covariates.
#' @param ... not used, only here for S3 compatibility
#' @return A n*p prediction matrix for new observations
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' predict(model, ex_data$X[1:3, , drop = FALSE])[, 1:5]
#' @export
predict.NormalBlockBase <- function(object, new_X, ...){
  stopifnot(isNB(object))
  with_dimnames(object$predict(new_X), rownames(new_X), colnames(object$data$Y))
}

#' @title Extract Log-Likelihood of a Normal-Block Model
#' @description Returns the (variational) log-likelihood of a fitted
#' normal-block model as a `"logLik"` object, compatible with [stats::AIC()]
#' and [stats::BIC()].
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return An object of class `"logLik"`. The numeric value is the
#' log-likelihood or its variational lower bound (ELBO), on the scale of the
#' data the model was fitted on (rescaled or not, see [NormalBlockData]); use
#' the `heldout_loglik()` method for a log-likelihood on Y's original scale.
#' Attributes `df` and `nobs` hold the number of parameters and observations.
#' @importFrom stats logLik
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' logLik(model)
logLik.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  structure(object$loglik, class = "logLik", df = object$nb_param, nobs = object$n)
}

#' @title Bayesian Information Criterion for a Normal-Block Model
#' @description Extracts the (variational) BIC of a fitted normal-block
#' model, computed as `deviance + penalty` (lower is better). The penalty
#' charges each group of parameters the log of the number of observations
#' informing it (see the `penalty_terms` field of [NormalBlockBase]), so it
#' differs from `stats::BIC()`'s default `log(nobs) * df` computed from
#' [logLik()].
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A scalar: the (variational) BIC.
#' @importFrom stats BIC
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' BIC(model)
BIC.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  object$BIC
}

#' @title Print a Normal-Block Model
#' @description Print a short summary of a fitted normal-block model: model
#' type, goodness-of-fit criteria, and the useful fields/methods to explore
#' it further.
#' @param x An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' print(model)
print.NormalBlockBase <- function(x, ...) {
  stopifnot(isNB(x))
  x$print()
  invisible(x)
}

#' @title Summarize a Normal-Block Model
#' @description Summarizes a fitted normal-block model: model type,
#' goodness-of-fit criteria, cluster sizes, and the density of the inferred
#' network between blocks.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return An object of class `summary.NormalBlockBase` (a list with the
#' model's `who_am_I`, `criteria`, `cluster_sizes` and network `density`),
#' printed with a dedicated [print.summary.NormalBlockBase()] method.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' summary(model)
summary.NormalBlockBase <- function(object, ...) {
  stopifnot(isNB(object))
  ## the network lives on Omega: q x q (variance-block) or p x p (mean-block)
  size <- ncol(object$model_par$Omega)
  res <- list(
    who_am_I      = object$who_am_I,
    criteria      = object$criteria,
    cluster_sizes = object$cluster_sizes,
    n_edges       = object$n_edges,
    density       = if (size > 1) object$n_edges / choose(size, 2) else NA
  )
  class(res) <- "summary.NormalBlockBase"
  res
}

#' @title Print a Normal-Block Model Summary
#' @description Print method for objects returned by [summary.NormalBlockBase()].
#' @param x An object of class `summary.NormalBlockBase`.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' print(summary(model))
#' @export
print.summary.NormalBlockBase <- function(x, ...) {
  cat(paste0("A ", x$who_am_I, ".\n"))
  cat("===========================================================================\n")
  print(as.data.frame(round(x$criteria, digits = 3), row.names = ""))
  cat("===========================================================================\n")
  cat_wrapped(paste(x$cluster_sizes, collapse = ", "), label = "* Cluster sizes: ")
  cat("* Network: ", x$n_edges, " edge(s), density = ",
      if (is.na(x$density)) NA else round(x$density, 3), "\n", sep = "")
  invisible(x)
}

#' @title Plot a Normal-Block Model
#' @description Plots the evolution of the objective (log-likelihood or
#' ELBO) across the (V)EM iterations of the last call to `optimize()`, see
#' `$plot_loglik()`.
#' @param x An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns the [ggplot2::ggplot] object; called for its
#' side effect of plotting.
#' @export
#' @examples
#' ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#' model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
#' plot(model)
plot.NormalBlockBase <- function(x, ...) {
  stopifnot(isNB(x))
  p <- x$plot()
  print(p)
  invisible(p)
}

## names of the clusters, as in the network plots
cluster_names <- function(q) paste0("Cluster_", seq_len(q))

## names of the nodes of a model's network: the clusters (variance-block) or
## the variables (mean-block)
network_nodes_names <- function(model) {
  if (inherits(model, "NormalBlockVarBase")) cluster_names(model$q) else colnames(model$data$Y)
}

## names the dimensions of a matrix when they match its size
with_dimnames <- function(M, rows = NULL, cols = NULL) {
  M <- as.matrix(M)
  if (length(rows) == nrow(M)) rownames(M) <- rows
  if (length(cols) == ncol(M)) colnames(M) <- cols
  M
}
