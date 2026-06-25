###############################################################################
## Regression tests for NormalBlockData's scale = TRUE default (R/NormalBlockData.R):
## each column of Y is divided by its own standard deviation (no centering --
## the model's own intercept in X is expected to absorb the mean) so that
## the normal-block model's shared within-block covariance assumption isn't
## swamped by variables with very different baseline variances.
set.seed(1)
n <- 40; p <- 6
Y <- matrix(rnorm(n * p), n, p) * matrix(c(1, 10, 100, 0.1, 1000, 1), n, p, byrow = TRUE)
X <- matrix(1, n, 1)

test_that("scale = TRUE (default) rescales each column to unit variance, without centering", {
  data <- NormalBlockData$new(Y, X)
  expect_equal(apply(data$Y, 2, sd), rep(1, p), tolerance = 1e-10)
  expect_equal(data$Y_scale, apply(Y, 2, sd))
  ## not centered: dividing by sd alone does not zero out the column means
  expect_false(isTRUE(all.equal(colMeans(data$Y), rep(0, p))))
  expect_equal(data$Y, Y / matrix(data$Y_scale, n, p, byrow = TRUE))
})

test_that("scale = FALSE leaves Y untouched", {
  data <- NormalBlockData$new(Y, X, scale = FALSE)
  expect_equal(data$Y, Y)
  expect_equal(data$Y_scale, rep(1, p))
})

test_that("scaling does not change which entries are zero (zeros/zeros_bar unaffected)", {
  Yz <- Y; Yz[1:5, 1] <- 0
  data_scaled   <- NormalBlockData$new(Yz, X, scale = TRUE)
  data_unscaled <- NormalBlockData$new(Yz, X, scale = FALSE)
  expect_equal(data_scaled$zeros_bar, data_unscaled$zeros_bar)
  expect_equal(data_scaled$nY, data_unscaled$nY)
})

test_that("a near-constant column does not produce Inf/NaN in Y or Y_scale", {
  Yconst <- Y; Yconst[, 1] <- 1 # zero variance
  data <- NormalBlockData$new(Yconst, X)
  expect_true(all(is.finite(data$Y)))
  expect_true(all(is.finite(data$Y_scale)))
})
