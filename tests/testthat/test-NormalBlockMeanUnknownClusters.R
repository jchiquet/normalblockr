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
