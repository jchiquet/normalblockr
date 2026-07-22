## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  CLASS NormalBlockData ##################################
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' R6 class for a generic normal model
#' @param Y the matrix of responses.
#' @param X design matrix.
#' @param X0 zero-inflation design matrix, if applicable.
#' @param formula describes the relationship between Y and X, and X0 if applicable, useful if not all of X's or X0's covariates should be used, should be formatted ~ X1 + X2... | Z1 + Z2... with the Normal formula before the | and the ZI formula after the |
#' @param scale whether to rescale each column of Y by its own standard
#' deviation before fitting (default TRUE). Columns are *not* centered: the
#' model's own intercept (the constant or group-indicator columns a user is
#' expected to include in X) already absorbs each variable's mean, so
#' centering here would be redundant -- unless X has no such column, in
#' which case the residual model is misspecified regardless of this setting.
#' Rescaling matters because the normal-block model assumes a single shared
#' covariance value within each block (`Cov(Y_j, Y_j') = Var(W_k)` for j, j'
#' in block k); on the raw scale, that assumption is swamped whenever
#' variables in the same true block have very different baseline variances
#' (e.g. species with very different total abundance). The resulting `B`,
#' `Omega` and clustering are therefore properties of the *scaled* data, not
#' directly convertible back to the original units (the per-column scaling
#' factors are not the same within a block, so there is no single way to
#' "unscale" a shared block covariance back to a p x p matrix).
#' @export
NormalBlockData <- R6::R6Class(
  classname = "NormalBlockData",
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## PUBLIC MEMBERS ----
  ## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  public = list(
    #' @field Y the matrix of responses (rescaled column-wise if `scale = TRUE`)
    Y  = NULL,
    #' @field Y_scale the per-column standard deviation Y was divided by (all
    #' 1's if `scale = FALSE`)
    Y_scale = NULL,
    #' @field X the matrix of covariates
    X = NULL,
    #' @field X0 the matrix of zero-inflation covariates, if applicable
    X0 = NULL,
    #' @field formula describes the relationship between Y and X, and X0 if applicable, useful if not all of X's or X0's covariates should be used, should be formatted ~ X1 + X2... | Z1 + Z2... with the Normal formula before the | and the ZI formula after the |
    formula = NULL,
    #' @field n sample size
    n  = NULL,
    #' @field d number of covariates
    d = NULL,
    #' @field d0 number of zero-inflation covariates, if applicable
    d0 = NULL,
    #' @field p number of variables
    p = NULL,
    #' @field XtXm1 inverse of XtX, useful for inference
    XtXm1 = NULL,
    #' @field XtY useful for inference
    XtY = NULL,
    #' @field npY total number of non zeros in Y
    npY = NULL,
    #' @field nY total number of non zeros for each column/variable in Y
    nY = NULL,
    #' @field zeros where are the zero in Y
    zeros = NULL,
    #' @field zeros_bar where are the non-zeros in Y
    zeros_bar = NULL,

    #' @description Create a new [`NormalBlockData`] object.
    #' @param Y the matrix of responses (called Y in the model).
    #' @param X design matrix (called X in the model).
    #' @param X0 zero-inflation design matrix, if applicable.
    #' @param formula describes the relationship between Y and X, useful if not all of X's covariates should be used.
    #' @param scale whether to rescale each column of Y by its own standard
    #' deviation (no centering). Default TRUE -- see the class-level
    #' documentation for the rationale and its limits.
    initialize = function(Y, X, X0 = NULL, formula = NULL, scale = TRUE) {
      stopifnot("Y and X must be matrices" = all(is.matrix(Y), is.matrix(X)))
      stopifnot("Y and X must have the same number of rows" = (nrow(Y) == nrow(X)))
      self$Y_scale <- if (scale) pmax(apply(Y, 2, sd), .Machine$double.eps) else rep(1, ncol(Y))
      self$Y <- Y / matrix(self$Y_scale, nrow(Y), ncol(Y), byrow = TRUE)
      self$n <- nrow(Y)
      self$p <- ncol(Y)
      self$formula <- formula
      fm_zi <- NULL
      if(is.null(formula)){self$X <- X
      }else{
        stopifnot("The formula should start with ~" = (formula[[1]] == "~"))
        if(formula[[2]][[1]] == "|"){
          fm <- as.formula(paste0("~", formula[[2]][2]))
          fm_zi <- as.formula(paste0("~", formula[[2]][3]))
        }else{
          fm <- formula}
        stopifnot("Covariates given in formula must be present in X" = (length(setdiff(all.vars(terms(fm)), colnames(X))) ==0))
        self$X <- model.matrix(fm, as.data.frame(X))
      }
      if(is.null(X0)){X0 <- matrix(rep(1, self$n))}
      if(is.null(fm_zi)){
        self$X0 <- X0
      }else{
        self$X0 <- model.matrix(fm_zi, as.data.frame(X0))
        stopifnot("Zero-inflation covariates given in the formula must be present in X0" = (length(setdiff(all.vars(terms(fm_zi)), colnames(X0))) ==0))
      }
      self$d0 <- ncol(self$X0)
      self$d <- ncol(self$X)
      self$XtXm1 <- solve(crossprod(self$X))
      self$XtY   <- crossprod(self$X, self$Y)
      self$zeros     <- 1 * (self$Y == 0)
      self$zeros_bar <- 1 * (self$Y != 0)
      self$npY <- sum(self$zeros_bar)
      self$nY  <- colSums(self$zeros_bar)
    }
  )
)
