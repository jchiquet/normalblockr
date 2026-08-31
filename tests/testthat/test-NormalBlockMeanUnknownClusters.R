###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
q <- ncol(C)
B <- testdata$parameters$B


test_that("normal block mean with unknown clusters - basic tests", {
  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockMeanUnknownClusters$new(data, q, 0.1)
  model$optimize()
  expect_lt(model$BIC, 11224)
})

test_that("NB_control(blas_threads = ) is applied during optimize() and restored afterwards", {
  skip_if_not_installed("RhpcBLASctl")
  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockMeanUnknownClusters$new(data, q, 0.1)

  RhpcBLASctl::blas_set_num_threads(3)
  on.exit(RhpcBLASctl::blas_set_num_threads(RhpcBLASctl::blas_get_num_procs()))
  model$optimize(control = NB_control(niter = 5, threshold = -1, verbose = FALSE, blas_threads = 1))
  expect_equal(RhpcBLASctl::blas_get_num_procs(), 3)
})
