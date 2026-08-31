## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockMeanCollectionSparsity ####################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Collection of Mean-Block Models over a Sparsity Path
#'
#' R6 class for a collection of mean-block models ([NormalBlockMeanBase]) with
#' a fixed clustering (or a fixed number of blocks) and different sparsity
#' levels applied to the p x p precision matrix of the variables. Mirrors
#' [NormalBlockVarCollectionSparsity], minus the StARS/stability selection
#' path, which relies on `fixed_tau` -- not supported by the mean-block VEM.
#' @examples
#' ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
#' data <- NormalBlockData$new(ex$Y, ex$X)
#' models <- normal_block(data, blocks = 3, sparsity = TRUE, model = "mean",
#'                        control = NB_control(n_sparsity_penalties = 5))
#' models$plot(c("BIC", "EBIC"))
#' @export
NormalBlockMeanCollectionSparsity <- R6::R6Class(
  classname = "NormalBlockMeanCollectionSparsity",
  inherit   = NormalBlockCollection,
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field data object of NormalBlockData class, with responses and design matrix
    data = NA,

    #' @description Create a new [`NormalBlockMeanCollectionSparsity`] object.
    #' @param mydata object of NormalBlockData class, with responses and design matrix
    #' @param blocks either a clustering matrix (known, fixed clustering) or a single integer (number of blocks to infer)
    #' @param control structured list of parameters to handle sparsity control
    #' @return A new [`NormalBlockMeanCollectionSparsity`] object
    initialize = function(mydata, blocks, control = NB_control()) {
      self$data    <- mydata
      self$control <- control
      private$blocks_ <- blocks

      stopifnot("blocks must be either a clustering matrix or a fixed number of blocks" =
                  is.matrix(blocks) | length(blocks) == 1)

      if (!is.null(control$sparsity_penalties)) {
        stopifnot("All penalties must be strictly positive" =
                    (min(control$sparsity_penalties) > 0))
        sparsity <- control$sparsity_penalties
      } else {
        ## The variance-block path reads its scale off a short unpenalized fit;
        ## here that fit is exactly what sparsity is meant to rescue when
        ## n <= p (Sigma singular), so the scale is taken directly from the OLS
        ## residual covariance -- the quantity the fitted Omega inverts anyway.
        Sigma    <- stats::cov(ols_residuals(mydata))
        weights  <- matrix(1, mydata$p, mydata$p)
        diag(weights) <- 0
        if (!is.null(control$sparsity_weights)) weights <- control$sparsity_weights
        diag_pen <- max(diag(weights)) > 0
        scale    <- abs((Sigma / weights)[upper.tri(Sigma, diag = diag_pen)])
        scale    <- scale[is.finite(scale)]
        if (length(scale) > 0) {
          max_pen  <- max(scale)
          sparsity <- 10^seq(log10(max_pen), log10(max_pen * control$min_ratio),
                             len = control$n_sparsity_penalties)
        } else {
          sparsity <- rep(0, len = control$n_sparsity_penalties)
        }
      }

      private$sparsity_ <- sort(sparsity, decreasing = TRUE)
      self$models <- map(private$sparsity_, function(lambda)
        get_model(mydata, blocks, sparsity = lambda, model = "mean", control = control)
      )
      private$progress_field <- "sparsity"
      private$progress_label <- "penalty"
    },

    #' @description optimizes every model in the sparsity path, warm-starting
    #' each one (after the first) from the previous, adjacent penalty's
    #' converged parameters (see [NormalBlockMeanBase]'s `warm_start_from()`).
    #' @param control optimization parameters (niter and threshold)
    optimize = function(control = list(niter = 500, threshold = 1e-4, verbose = TRUE)) {
      previous <- NULL
      self$models <- lapply(self$models, function(model) {
        if (control$verbose)
          cat("\t", private$progress_label, "=", model[[private$progress_field]], "          \r")
        flush.console()
        if (!is.null(previous)) model$warm_start_from(previous)
        model$optimize(control)
        previous <<- model
        model
      })
      invisible(self)
    },

    #' @description returns the model corresponding to a given penalty
    #' @param sparsity sparsity penalty asked by user
    #' @return a [`NormalBlockMeanBase`] object with the given penalty
    get_model = function(sparsity) {
      if (!(sparsity %in% private$sparsity_)) {
        sparsity <- private$sparsity_[[which.min(abs(private$sparsity_ - sparsity))]]
        message("No model with this penalty in the collection. Returning model with closest penalty: ",
                sparsity, " Collection penalty values can be found via $sparsity")
      }
      self$models[[which(private$sparsity_ == sparsity)]]
    },

    #' @description Extract best model in the collection
    #' @param crit a character for the criterion used to perform the selection,
    #' either "BIC", "EBIC" or "ICL". Default is BIC
    #' @return a [`NormalBlockMeanBase`] object
    get_best_model = function(crit = c("BIC", "EBIC", "ICL")) {
      crit <- match.arg(crit)
      self$models[[private$best_id(crit)]]$clone()
    },

    #' @description Display goodness-of-fit criteria along the sparsity path
    #' @param criteria vector of characters. The criteria to plot in `c("deviance", "BIC", "EBIC", "ICL")`. Defaults to all of them.
    #' @param log.x logical: should the x-axis be represented in log-scale? Default is `TRUE`.
    #' @return a [`ggplot2::ggplot`] graph
    plot = function(criteria = c("deviance", "BIC", "EBIC", "ICL"), log.x = TRUE) {
      stopifnot(!is.null(self$criteria[criteria]))
      p <- private$plot_criteria_path("sparsity", criteria)
      if (log.x) p <- p + ggplot2::coord_trans(x = "log10")
      p
    }
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PRIVATE MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  private = list(
    sparsity_ = NA, # penalty values in the collection (lambda)
    blocks_   = NA  # blocks (either a scalar or an indicator matrix)
  ),

  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ##  ACTIVE BINDINGS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  active = list(
    #' @field q number of blocks
    q = function() ifelse(is.matrix(private$blocks_), ncol(private$blocks_), private$blocks_),
    #' @field blocks group matrix or number of blocks
    blocks = function() private$blocks_,
    #' @field sparsity list of sparsity penalties
    sparsity = function() private$sparsity_,
    #' @field who_am_I a method to print what model is being fitted
    who_am_I = function() "normal-block-mean model with sparsity path"
  )
)
