###############################################################################
## Exercises get_model()/normal_block()'s dispatch logic: every combination of
## `blocks`/`sparsity`/`zero_inflation` should reach the class the dispatcher
## is meant to build. Numerical accuracy of the fit is not the point here
## (covered by the other test files) -- niter is kept tiny on purpose.
set.seed(7)
ex   <- generate_normal_block_data(n = 40, p = 8, d = 1, q = 2)
data <- NormalBlockData$new(ex$Y, ex$X)
C    <- ex$parameters$C
fast <- NB_control(niter = 2, threshold = -1, verbose = FALSE)

test_that("normal_block() dispatches known-clusters models (matrix blocks)", {
  m1 <- normal_block(data, blocks = C, zero_inflation = FALSE, control = fast)
  expect_s3_class(m1, "NormalBlockKnownClusters")
  m2 <- normal_block(data, blocks = C, zero_inflation = TRUE, control = fast)
  expect_s3_class(m2, "ZINormalBlockKnownClusters")
})

test_that("normal_block() dispatches unknown-clusters models (single integer blocks)", {
  m1 <- normal_block(data, blocks = 2, zero_inflation = FALSE, control = fast)
  expect_s3_class(m1, "NormalBlockUnknownClusters")
  m2 <- normal_block(data, blocks = 2, zero_inflation = TRUE, control = fast)
  expect_s3_class(m2, "ZINormalBlockUnknownClusters")
})

test_that("normal_block() dispatches NormalBlockChangingSparsity for sparsity = TRUE", {
  control <- fast
  control$sparsity_penalties <- c(0.1, 0.05)
  m1 <- normal_block(data, blocks = 2, sparsity = TRUE, zero_inflation = FALSE, control = control)
  expect_s3_class(m1, "NormalBlockChangingSparsity")
  ## zero_inflation propagates to the models in the collection, not to the
  ## collection's own class name
  m2 <- normal_block(data, blocks = C, sparsity = TRUE, zero_inflation = TRUE, control = control)
  expect_s3_class(m2, "NormalBlockChangingSparsity")
  expect_s3_class(m2$models[[1]], "ZINormalBlockKnownClusters")
})

test_that("normal_block() dispatches NormalBlockUnknownQ for a vector of candidate block counts", {
  m1 <- normal_block(data, blocks = c(2, 3), zero_inflation = FALSE, control = fast)
  expect_s3_class(m1, "NormalBlockUnknownQ")
  m2 <- normal_block(data, blocks = c(2, 3), zero_inflation = TRUE, control = fast)
  expect_s3_class(m2, "NormalBlockUnknownQ")
  expect_s3_class(m2$models[[1]], "ZINormalBlockUnknownClusters")
})

test_that("normal_block() dispatches NormalBlockUnknownQChangingSparsity for sparsity = TRUE and a vector of block counts", {
  control <- fast
  control$sparsity_penalties <- c(0.1, 0.05)
  m1 <- normal_block(data, blocks = c(2, 3), sparsity = TRUE, zero_inflation = FALSE, control = control)
  expect_s3_class(m1, "NormalBlockUnknownQChangingSparsity")
})
