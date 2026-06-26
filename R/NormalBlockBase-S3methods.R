## =========================================================================================
##
## PUBLIC S3 METHODS FOR NormalBlockBase
##
## =========================================================================================

#' @title Checks if a model is of class [NormalBlockBase()]
#' @description Checks if a model is of class [NormalBlockBase()]
#' @param object An R object.
#' @return A boolean telling whether object inherits from the NormalBlockBase class.
#' @export
isNB <- function(object) {inherits(object, "NormalBlockBase")}

#' @title Extracts model coefficients from objects returned by [NormalBlockBase()] and its variants
#' @description Extract coefficients from an NormalBlockBase object.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A matrix of coefficients extracted from the NormalBlockBase model.
#' @export
coef.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  object$model_par$B
}

#' @title Extracts the latent-block covariance matrix from objects returned by [NormalBlockBase()] and its variants
#' @description Extract the covariance matrix between latent blocks (the inverse of `Omegaq`) from a NormalBlockBase object.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return The q x q covariance matrix between latent blocks.
#' @importFrom stats sigma
#' @export
sigma.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  solve(object$model_par$Omegaq)
}

#' @title Extracts fitted values from objects returned by [NormalBlockBase()] and its variants
#' @description Extract fitted values from an NormalBlockBase object.
#' @param object An object of class NormalBlockBase.
#' @param ... not used, only here for S3 compatibility
#' @return A matrix of fitted values extracted from the object.
#' @export
fitted.NormalBlockBase <- function(object, ...){
  stopifnot(isNB(object))
  object$fitted
}


#' @title Predicts observations Y for new covariates X.
#' @description Predicts observations Y for new covariates X.
#' @param object An object of class NormalBlockBase.
#' @param new_X New set of covariates.
#' @param ... not used, only here for S3 compatibility
#' @return A n*p prediction matrix for new observations
#' @export
predict.NormalBlockBase <- function(object, new_X, ...){
  stopifnot(isNB(object))
  object$predict(new_X)
}
