###############################################################################
###############################################################################
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C
B <- testdata$parameters$B


test_that("normal block mean with known clusters - basic tests", {
  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockMeanKnownClusters$new(data, C, 0.1)
  model$optimize()
  expect_lt(model$BIC, 15400)
})
