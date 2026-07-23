###############################################################################
###############################################################################

testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)

test_that("Robustness of starting clustering with NormalBlockVarBase", {
  ## Diagonal model
  model_kmeans <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "kmeans"))
  model_kmeans$optimize()

  model_ward2 <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "ward2"))
  model_ward2$optimize()

  model_sbm <- NormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2650)
  expect_gt(model_ward2$loglik, -2650)
  expect_gt(model_sbm$loglik, -2650)

  # expect_gte(model_kmeans$loglik, model_sbm$loglik)
  # expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gt(aricode::ARI(model_kmeans$clustering, model_sbm$clustering), 0.75)
  expect_gt(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.75)
})


testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)

test_that("Robustness of starting clustering with ZINB", {
  ## Diagonal model
  model_ward2 <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "ward2"))
  model_ward2$optimize()

  model_kmeans <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "kmeans"))
  model_kmeans$optimize()

  model_sbm <- ZINormalBlockVarUnknownClusters$new(data, q, control = NB_control(clustering_init = "sbm"))
  model_sbm$optimize()

  expect_gt(model_kmeans$loglik, -2750)
  expect_gt(model_ward2$loglik, -2750)
  expect_gt(model_sbm$loglik, -2750)

  expect_gte(model_kmeans$loglik, model_sbm$loglik)
  expect_gte(model_ward2$loglik, model_sbm$loglik)

  expect_gte(aricode::ARI(model_kmeans$clustering, model_sbm$clustering)  , 0.49)
  expect_gte(aricode::ARI(model_kmeans$clustering, model_ward2$clustering), 0.99)
})
