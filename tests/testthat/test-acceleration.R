###############################################################################
## NB_control(acceleration = ): SQUAREM on top of (V)EM (src/normal_block_em_base.h)
## can be turned off to run plain (V)EM, in both families. Plain (V)EM must be
## monotone, reach the same optimum when the clustering is known, and need at
## least as many iterations.
testdata_var  <- readRDS("testdata/testdata_normal.RDS")
testdata_mean <- readRDS("testdata/testdata_normal_mean_block.RDS")
data_var  <- NormalBlockData$new(testdata_var$Y, testdata_var$X)
data_mean <- NormalBlockData$new(testdata_mean$Y, testdata_mean$X)

fit_both <- function(data, blocks, model) {
  lapply(c(acc = TRUE, plain = FALSE), function(acc) {
    set.seed(1)
    normal_block(data, blocks, model = model,
                 control = NB_control(verbose = FALSE, acceleration = acc))
  })
}

test_that("NB_control() carries the acceleration switch, on by default", {
  expect_true(NB_control()$acceleration)
  expect_false(NB_control(acceleration = FALSE)$acceleration)
})

test_that("variance-block, known clusters: plain EM is monotone and reaches the accelerated optimum in more iterations", {
  fits <- fit_both(data_var, testdata_var$parameters$C, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_equal(fits$plain$loglik, fits$acc$loglik, tolerance = 1e-6)
  expect_gt(fits$plain$criteria$niter, fits$acc$criteria$niter)
})

test_that("variance-block, unknown clusters: the switch reaches the VEM", {
  fits <- fit_both(data_var, 3, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_gt(fits$plain$criteria$niter, fits$acc$criteria$niter)
})

test_that("mean-block: plain (V)EM is monotone and reaches the same optimum", {
  known <- fit_both(data_mean, testdata_mean$parameters$C, "mean")
  expect_true(all(diff(known$plain$objective) > -1e-8))
  expect_equal(known$plain$loglik, known$acc$loglik, tolerance = 1e-6)

  unknown <- fit_both(data_mean, 3, "mean")
  expect_true(all(diff(unknown$plain$objective) > -1e-8))
  expect_equal(unknown$plain$loglik, unknown$acc$loglik, tolerance = 1e-6)
})

test_that("zero-inflated models accept the switch", {
  zi <- readRDS("testdata/testdata_normal_zi.RDS")
  data_zi <- NormalBlockData$new(zi$Y, zi$X)
  fits <- fit_both(data_zi, zi$parameters$C, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_equal(fits$plain$loglik, fits$acc$loglik, tolerance = 1e-5)
})
