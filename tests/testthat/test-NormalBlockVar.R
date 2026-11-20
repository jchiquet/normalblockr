###############################################################################
## Variance-block models (NormalBlockVar{Known,Unknown}Clusters): diagonal and
## spherical residual covariance, with and without sparsity, known clustering
## or recovered one, and setting $sparsity then re-optimizing.
## Frozen fixture (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)
rmse <- function(model) sqrt(mean((model$fitted - Y)^2))

test_that("every combination fits and an unknown clustering is recovered", {
  min_loglik <- c(diagonal = -950, spherical = -1020)
  for (noise in c("diagonal", "spherical")) for (known in c(TRUE, FALSE)) for (sparsity in c(0, 0.05)) {
    set.seed(1)
    model <- normal_block(data, if (known) C else q, sparsity = sparsity,
                          control = NB_control(verbose = FALSE, noise_covariance = noise))
    expect_s3_class(model, if (known) "NormalBlockVarKnownClusters" else "NormalBlockVarUnknownClusters")
    expect_equal(model$sparsity, sparsity)
    expect_gt(model$loglik, min_loglik[[noise]])
    expect_lt(rmse(model), 1.1)
    if (sparsity == 0) expect_true(all(diff(model$objective) > -1e-5))
    if (!known) expect_equal(aricode::ARI(model$clustering, get_clusters(C)), 1)
  }
})

test_that("setting $sparsity and re-optimizing reaches the direct sparse fit", {
  direct <- NormalBlockVarKnownClusters$new(data, C, sparsity = 0.05)
  direct$optimize()
  model <- NormalBlockVarKnownClusters$new(data, C)
  model$optimize()
  model$sparsity <- 0.05
  expect_equal(model$sparsity, 0.05)
  model$optimize()
  expect_equal(model$loglik, direct$loglik, tolerance = 1e-2)
  expect_equal(model$model_par$Omega, direct$model_par$Omega, tolerance = 1e-2)
  expect_error(model$sparsity <- -1, "positive")
})
