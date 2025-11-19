## =========================================================================================
##
## PUBLIC S3 METHODS FOR NB
##
## =========================================================================================

#' @title Checks if a model is of class [NB()]
#' @description Checks if a model is of class [NB()]
#' @param object An R object.
#' @return A boolean telling whether object inherits from the NB class.
#' @export
isNB <- function(object) {inherits(object, "NB")}

#' @title Extracts model coefficients from objects returned by [NB()] and its variants
#' @description Extract coefficients from an NB object.
#' @param object An object of class NB.
#' @return A matrix of coefficients extracted from the NB model.
#' @export
coef.NB <- function(object, ...){
  stopifnot(isNB(object))
  object$model_par$B
}

#' @title Extracts model coefficients from objects returned by [NB()] and its variants
#' @description Extract coefficients from an NB object.
#' @param object An object of class NB.
#' @return A matrix of coefficients extracted from the NB model.
#' @export
sigma.NB <- function(object, ...){
  stopifnot(isNB(object))
  solve(object$model_par$Omegaq)
}

#' @title Extracts fitted values from objects returned by [NB()] and its variants
#' @description Extract fitted values from an NB object.
#' @param object An object of class NB.
#' @return A matrix of Fitted values extracted from the object object.
#' @export
fitted.NB <- function(object, ...){
  stopifnot(isNB(object))
  object$fitted
}


#' @title Predicts observations Y for new covariates X.
#' @description Predicts observations Y for new covariates X.
#' @param object An object of class NB.
#' @param new_X New set of covariates.
#' @return A n*p prediction matrix for new observations
#' @export
predict.NB <- function(object, new_X){
  stopifnot(isNB(object))
  object$fitted
}

