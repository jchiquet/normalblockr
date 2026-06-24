###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_covardep_zi.RDS")
Y  <- testdata$Y
X  <- testdata$X
X0 <- testdata$X0
C  <- testdata$parameters$C ; q <- ncol(C)
data  <- NormalBlockData$new(Y, X, X0 = X0)

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  ## Diagonal model
  ## clustering_approx pinned to "ward2": the thresholds below were tuned for
  ## it specifically, and other heuristics (e.g. the "kmeans" default) land in
  ## a different, worse local optimum on this particular small dataset (see
  ## test-clustering-heuristics.R for the broader point that no heuristic
  ## dominates across datasets/models).
  model <- ZINormalBlockUnknownClusters$new(data, q, control = NB_control(clustering_approx = "ward2"))
  model$optimize()
  expect_lt(model$BIC, 4300)
  expect_gt(model$loglik, -2000)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.55)

  model <- ZINormalBlockKnownClusters$new(data, C, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -2000)
  expect_lt(Metrics::rmse(model$fitted, Y), 0.55)
})
