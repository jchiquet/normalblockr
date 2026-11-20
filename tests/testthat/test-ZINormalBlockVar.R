###############################################################################
## Zero-inflated variance-block models (ZINormalBlockVar{Known,Unknown}Clusters),
## with an intercept-only zero-inflation (given through a formula) or one
## depending on covariates X0: diagonal and spherical residual covariance,
## with and without sparsity, known clustering or recovered one.
## Frozen fixtures (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X ; colnames(X) <- "X1"
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X, formula = ~ 0 + X1 | 1)

test_that("intercept-only zero-inflation: every combination fits, an unknown clustering is recovered", {
  for (noise in c("diagonal", "spherical")) for (known in c(TRUE, FALSE)) for (sparsity in c(0, 2)) {
    set.seed(1)
    model <- normal_block(data, if (known) C else q, sparsity = sparsity, zero_inflation = TRUE,
                          control = NB_control(verbose = FALSE, noise_covariance = noise))
    expect_s3_class(model, if (known) "ZINormalBlockVarKnownClusters" else "ZINormalBlockVarUnknownClusters")
    expect_gt(model$loglik, -2050)
    expect_lt(model$BIC, 4250)
    expect_lt(sqrt(mean((model$fitted - Y)^2)), 0.9)
    expect_true(all(model$fitted[data$zeros == 1] == 0))
    if (!known) expect_equal(aricode::ARI(model$clustering, get_clusters(C)), 1)
  }
})

test_that("zero-inflation depending on covariates X0", {
  td <- readRDS("testdata/testdata_normal_covardep_zi.RDS")
  d  <- NormalBlockData$new(td$Y, td$X, X0 = td$X0)
  Ctrue <- td$parameters$C
  for (sparsity in c(0, 2)) {
    model <- normal_block(d, Ctrue, sparsity = sparsity, zero_inflation = TRUE,
                          control = NB_control(verbose = FALSE))
    expect_gt(model$loglik, -1450)
    expect_lt(model$BIC, 3150)
    expect_lt(sqrt(mean((model$fitted - td$Y)^2)), 0.75)
  }
  ## kmeans lands in a worse local optimum on these small data
  set.seed(1)
  model <- normal_block(d, ncol(Ctrue), zero_inflation = TRUE,
                        control = NB_control(verbose = FALSE, clustering_init = "ward2"))
  expect_gt(model$loglik, -1450)
  expect_lt(model$BIC, 3200)
  expect_lt(sqrt(mean((model$fitted - td$Y)^2)), 0.8)
  expect_equal(dim(model$model_par$B0), c(ncol(td$X0), ncol(td$Y)))
})
