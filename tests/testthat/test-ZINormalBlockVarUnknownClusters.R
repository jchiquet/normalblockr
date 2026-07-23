###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)

###############################################################################

test_that("zero inflated normal block with diagonal residual covariance and known clusters", {
  ## Diagonal model
  model <- ZINormalBlockVarUnknownClusters$new(data, q)
  model$optimize()
  expect_lt(model$BIC, 5700)
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)

  model <- ZINormalBlockVarUnknownClusters$new(data, q, sparsity = 2)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)

})

test_that("zero inflated normal block with spherical residual covariance and known clusters", {
  ## Spherical model
  ctrl <- NB_control(noise_covariance = "spherical")
  model <- ZINormalBlockVarUnknownClusters$new(data, q, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)

  model <- ZINormalBlockVarUnknownClusters$new(data, q, sparsity = 0.1, control = ctrl)
  model$optimize()
  expect_gt(model$loglik, -2700)
  expect_lt(Metrics::rmse(model$fitted, Y), 3)
})

test_that("zero inflated normal block with known clusters, heuristic", {
  model <- ZINormalBlockVarUnknownClusters$new(data, q, sparsity = 2,
                                  control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_lt(Metrics::rmse(model$fitted, Y), 3)
})
