###############################################################################
## Regression test for zi_weighted_fit() (R/utils.R): a variable observed at
## very few non-zero points relative to X's degrees of freedom can be fit
## exactly on one WLS iterate, driving its residual sum of squares -- and so
## dm1 = nY/ssq -- to Inf; on the next iterate that Inf weight propagates
## into XtWX and MASS::ginv()'s svd() errors on the non-finite input (found
## via inst/onema/study_ONEMA.R: a rare species seen at a single station).
set.seed(1)
n <- 50; p <- 5; d <- 1
Y <- matrix(rexp(n * p), n, p)
Y[, 1] <- 0
Y[1, 1] <- 5 # a single non-zero observation for variable 1
X <- matrix(1, n, d)
data <- NormalBlockData$new(Y, X)

test_that("zi_weighted_fit() never produces non-finite dm1/B, even for a variable with a single non-zero observation", {
  fit <- zi_weighted_fit(data)
  expect_true(all(is.finite(fit$dm1)))
  expect_true(all(is.finite(fit$B)))
})

test_that("a ZI model with such a variable optimizes without erroring", {
  model <- ZINormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  expect_no_error(model$optimize(control = NB_control(verbose = FALSE, niter = 5)))
})
