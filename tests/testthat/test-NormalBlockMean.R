###############################################################################
## Mean-block models (NormalBlockMean{Known,Unknown}Clusters): every shape of
## Sigma, a penalized full Sigma, the recovery of an unknown clustering, and
## fixed_tau. See test-cpp-normal-block-mean.R for the numerical accuracy of
## the recursion.
## Frozen fixture (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)
data <- NormalBlockData$new(Y, X)

test_that("every shape of Sigma fits, monotonically, and an unknown clustering is recovered", {
  max_BIC <- c(diagonal = 17200, spherical = 17100, full = 18000)
  min_ARI <- c(diagonal = 0.85, spherical = 0.7, full = 0.85)
  for (noise in names(max_BIC)) for (known in c(TRUE, FALSE)) {
    set.seed(1)
    model <- normal_block(data, if (known) C else q, model = "mean",
                          control = NB_control(verbose = FALSE, noise_covariance = noise))
    expect_s3_class(model, if (known) "NormalBlockMeanKnownClusters" else "NormalBlockMeanUnknownClusters")
    expect_lt(model$BIC, max_BIC[[noise]])
    expect_true(all(diff(model$objective) > -1e-6))
    expect_equal(dim(model$model_par$B), c(ncol(X), q))
    if (noise != "full") expect_equal(model$model_par$Omega, diag(diag(model$model_par$Omega)))
    if (!known) expect_gt(aricode::ARI(model$clustering, get_clusters(C)), min_ARI[[noise]])
  }
})

test_that("a penalized full Sigma fits, with or without a known clustering", {
  known <- NormalBlockMeanKnownClusters$new(data, C, sparsity = 0.1)
  known$optimize()
  expect_lt(known$BIC, 17500)
  set.seed(1)
  unknown <- NormalBlockMeanUnknownClusters$new(data, q, sparsity = 0.1)
  unknown$optimize()
  expect_lt(unknown$BIC, 17500)
  ## the same fit on column-rescaled data recovers the clusters much worse
  expect_gt(aricode::ARI(unknown$clustering, get_clusters(C)), 0.8)
})

test_that("NB_control(fixed_tau = TRUE) leaves tau at its initial value", {
  init <- NormalBlockMeanUnknownClusters$new(data, q, control = NB_control(verbose = FALSE)
            )$.__enclos_env__$private$optim_initialize()
  model <- NormalBlockMeanUnknownClusters$new(data, q,
             control = NB_control(verbose = FALSE, clustering_init = get_clusters(init$tau)))
  model$optimize(NB_control(verbose = FALSE, fixed_tau = TRUE))
  expect_equal(model$var_par$tau, init$tau)
})
