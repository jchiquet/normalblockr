## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CROSS-VALIDATED LOG-LIKELIHOOD #####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Cross-Validated Log-Likelihood of Normal-Block Models
#'
#' Splits the observations into folds, fits [normal_block()] on all but one
#' fold and evaluates the held-out log-likelihood of the left-out fold (see
#' `heldout_loglik()` in [NormalBlockBase]), for every model of the
#' collection. It selects the number of clusters without a penalty, and
#' compares models of different families, shapes of Sigma or scalings.
#'
#' @param data a [NormalBlockData] object
#' @param blocks as in [normal_block()]: a range of numbers of clusters, a
#' single one, or a known clustering (p x q indicator matrix)
#' @param model,zero_inflation,sparsity as in [normal_block()]. A sparsity
#' path (`sparsity = TRUE`) is not supported: its penalty grid depends on
#' the data, hence on the fold; pick the penalty with StARS instead (see
#' [NormalBlockVarCollectionSparsity]).
#' @param folds the number of folds, or a vector of `data$n` fold labels.
#' With a number, the folds are drawn at random, stratified on the distinct
#' rows of the design X when it is categorical (e.g. a factor coded by
#' `model.matrix()`), so that every training set sees every category.
#' @param control a list of controls, see [NB_control()]
#' @return an object of class `normal_block_cv`: a data frame with one row
#' per model (`q`), the total held-out log-likelihood over the `n`
#' observations (`loglik`, larger is better) and its standard error
#' (`se`, from the spread of the per-observation values); the
#' per-observation values are kept in the `per_observation` attribute
#' (n x number of models) and the folds in `folds`. Its `best` attribute
#' is the q chosen by the one-standard-error rule, the smallest q within one
#' standard error of the maximum; `best_max` is the maximum itself, which
#' tends to over-select.
#' @examples
#' ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' set.seed(1)
#' cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
#' cv
#' @export
normal_block_cv <- function(data, blocks, model = c("var", "mean"), zero_inflation = FALSE,
                            sparsity = 0, folds = 5, control = NB_control(verbose = FALSE)) {
  model <- match.arg(model)
  stopifnot("data must be a NormalBlockData object" = inherits(data, "NormalBlockData"),
            "a sparsity path is not supported (data-dependent grid): use StARS" = !isTRUE(sparsity),
            "heuristic fits have no likelihood" = !isTRUE(control$heuristic))
  n <- data$n
  if (length(folds) == 1) {
    ## stratified on the distinct rows of X when the design is categorical,
    ## so that every training set sees every category
    design <- match(do.call(paste, as.data.frame(data$X)), unique(do.call(paste, as.data.frame(data$X))))
    K <- folds
    folds <- integer(n)
    if (max(design) <= n / 2) {
      for (g in unique(design)) folds[design == g] <- sample(rep_len(sample(K), sum(design == g)))
    } else {
      folds <- sample(rep_len(seq_len(K), n))
    }
  }
  stopifnot("folds must give one label per observation" = length(folds) == n,
            "at least two folds are needed" = length(unique(folds)) >= 2)

  ## the data of a fold, with the same scale request
  rows_data <- function(rows, k, what) {
    tryCatch(data$subset(rows, check_constant = what == "training"), error = function(e)
      stop("fold ", k, " (", what, " set): ", conditionMessage(e), call. = FALSE))
  }

  per_obs <- NULL
  for (k in unique(folds)) {
    test <- folds == k
    fit <- normal_block(rows_data(!test, k, "training"), blocks, sparsity = sparsity,
                        zero_inflation = zero_inflation, control = control, model = model)
    models <- if (isNBcollection(fit)) fit$models else list(fit)
    if (is.null(per_obs)) {
      per_obs <- matrix(NA_real_, n, length(models),
                        dimnames = list(NULL, paste0("q", vapply(models, function(m) m$q, numeric(1)))))
    }
    test_data <- rows_data(test, k, "test")
    per_obs[test, ] <- vapply(models, function(m) m$heldout_loglik(test_data, per_observation = TRUE),
                              numeric(sum(test)))
  }
  res <- data.frame(q      = as.numeric(sub("^q", "", colnames(per_obs))),
                    loglik = colSums(per_obs),
                    se     = apply(per_obs, 2, stats::sd) * sqrt(n),
                    row.names = NULL)
  i_max <- which.max(res$loglik)
  best  <- min(res$q[res$loglik >= res$loglik[i_max] - res$se[i_max]])
  structure(res, per_observation = per_obs, folds = folds, best = best, best_max = res$q[i_max],
            class = c("normal_block_cv", "data.frame"))
}

#' @title Print a Cross-Validated Log-Likelihood
#' @param x an object of class `normal_block_cv`, from [normal_block_cv()]
#' @param ... not used, only here for S3 compatibility
#' @return Invisibly returns `x`.
#' @examples
#' ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' set.seed(1)
#' cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
#' print(cv)
#' @export
print.normal_block_cv <- function(x, ...) {
  cat("Cross-validated log-likelihood,", length(unique(attr(x, "folds"))), "folds",
      "(larger is better)\n")
  print(as.data.frame(unclass(x))[c("q", "loglik", "se")], row.names = FALSE)
  cat("Best: q =", attr(x, "best"), "(one-standard-error rule; maximum at q =",
      paste0(attr(x, "best_max"), ")\n"))
  invisible(x)
}

#' @title Plot a Cross-Validated Log-Likelihood
#' @description The held-out log-likelihood against the number of clusters,
#' with one standard error on either side.
#' @param x an object of class `normal_block_cv`, from [normal_block_cv()]
#' @param ... not used, only here for S3 compatibility
#' @return a [`ggplot2::ggplot`] graph
#' @examples
#' ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' set.seed(1)
#' cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
#' plot(cv)
#' @export
plot.normal_block_cv <- function(x, ...) {
  d <- as.data.frame(unclass(x))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$q, y = .data$loglik)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$loglik - .data$se, ymax = .data$loglik + .data$se),
                           width = 0.2) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::labs(x = "number of clusters q", y = "held-out log-likelihood") +
    ggplot2::theme_bw()
}
