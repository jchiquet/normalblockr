###############################################################################
## normal_block_sequential() fits a mean-block model, then a variance-block one
## on its residuals. The positive control below carries two genuinely distinct
## structures (unrelated partitions in the mean and in the covariance); the
## estimator must recover both.

test_that("normal_block_sequential() recovers two distinct structures", {
  set.seed(42)
  n <- 150; p <- 40; d <- 3; q_mean <- 5; q_var <- 4
  X  <- cbind(1, matrix(rnorm(n * (d - 1)), n, d - 1))
  Cm <- normalblockr:::as_indicator(sample(rep(1:q_mean, length.out = p)))
  Cv <- normalblockr:::as_indicator(sample(rep(1:q_var,  length.out = p)))
  Y  <- X %*% matrix(rnorm(d * q_mean, sd = 3), d, q_mean) %*% t(Cm) +
        matrix(rnorm(n * q_var, sd = 1.5), n, q_var) %*% t(Cv) +
        matrix(rnorm(n * p, sd = 0.5), n, p)

  fit <- normal_block_sequential(NormalBlockData$new(Y, X), blocks_mean = 3:7, blocks_var = 2:6)

  expect_s3_class(fit, "normal_block_sequential")
  expect_equal(fit$mean$q, q_mean)
  expect_equal(fit$var$q,  q_var)
  expect_gt(aricode::ARI(fit$mean$clustering, normalblockr:::get_clusters(Cm)), 0.95)
  expect_gt(aricode::ARI(fit$var$clustering,  normalblockr:::get_clusters(Cv)), 0.95)
  ## the two structures are unrelated by construction, and must stay so
  expect_lt(aricode::ARI(fit$mean$clustering, fit$var$clustering), 0.2)
  expect_equal(fit$residuals, Y - fitted(fit$mean), ignore_attr = TRUE)
})

test_that("normal_block_sequential() accepts fixed clusterings and prints", {
  testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
  data <- NormalBlockData$new(testdata$Y, testdata$X)
  C    <- testdata$parameters$C

  set.seed(1)
  fit <- normal_block_sequential(data, blocks_mean = C, blocks_var = 2)
  expect_s3_class(fit$mean, "NormalBlockMeanKnownClusters")
  expect_equal(fit$var$q, 2)
  expect_output(print(fit), "sequential mean-then-variance")
  capture.output(expect_invisible(print(fit)))
})

test_that("the zero-inflated chain runs and shares one mask across both stages", {
  set.seed(7)
  ex <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
  Y  <- ex$Y; Y[runif(length(Y)) < 0.2] <- 0
  data <- NormalBlockData$new(Y, ex$X)

  set.seed(1)
  fit <- normal_block_sequential(data, blocks_mean = 3, blocks_var = 2,
                                 zero_inflation = TRUE)

  expect_s3_class(fit, "normal_block_sequential")
  expect_true(inherits(fit$mean, "ZINormalBlockMeanUnknownClusters"))
  expect_true(inherits(fit$var, "ZINormalBlockVarUnknownClusters"))
  ## the residuals carry no zero of their own, so only the inherited mask can
  ## have told the second stage where the structural zeros are
  expect_equal(fit$var$data$zeros, data$zeros)
  expect_true(all(fit$residuals[data$zeros == 1] == 0))
})

test_that("normal_block_sequential() can start from a mean-block fit, without refitting it", {
  set.seed(3)
  ex   <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
  data <- NormalBlockData$new(ex$Y, ex$X)
  ctrl <- NB_control(verbose = FALSE)
  set.seed(1)
  coll <- normal_block(data, blocks = 2:4, model = "mean", control = ctrl)

  ## from a collection: its best model by crit is the first stage, unchanged
  set.seed(4)
  fit <- normal_block_sequential(data, blocks_mean = coll, blocks_var = 2, crit = "BIC")
  best <- coll$get_best_model("BIC")
  expect_equal(fit$mean$clustering, best$clustering)
  expect_equal(fit$mean$loglik, best$loglik)
  expect_equal(fit$residuals, ex$Y - fitted(best), ignore_attr = TRUE)

  ## from a single model: the same second stage as when the first is refitted
  single <- normal_block(data, blocks = 3, model = "mean", control = ctrl)
  set.seed(5)
  from_fit <- normal_block_sequential(data, blocks_mean = single, blocks_var = 2)
  expect_equal(from_fit$mean$loglik, single$loglik)
  expect_equal(from_fit$var$q, 2)

  ## a fit on other data, or with another zero-inflation setting, is refused
  other <- NormalBlockData$new(ex$Y[1:40, ], ex$X[1:40, , drop = FALSE])
  expect_error(normal_block_sequential(other, blocks_mean = single, blocks_var = 2), "obtained on `data`")
  expect_error(normal_block_sequential(data, blocks_mean = single, blocks_var = 2, zero_inflation = TRUE),
               "zero_inflation")
})

test_that("the second stage gets the scale request of the data", {
  set.seed(6)
  ex <- generate_normal_block_mean_data(n = 60, p = 12, d = 1, q = 2)
  Y  <- ex$Y * matrix(c(1, 5), 60, 12, byrow = TRUE)
  for (scale in list(NULL, TRUE, FALSE)) {
    data <- NormalBlockData$new(Y, ex$X, scale = scale)
    set.seed(1)
    fit <- normal_block_sequential(data, blocks_mean = ex$parameters$C, blocks_var = 2)
    ## the residuals are on Y's original scale whatever the request
    expect_equal(fit$residuals, Y - fitted(fit$mean), ignore_attr = TRUE)
    ## the variance-block family rescales unless told not to, the mean-block one only when told to
    expect_equal(all(fit$var$data$Y_scale == 1), isFALSE(scale))
    expect_equal(all(fit$mean$data$Y_scale == 1), !isTRUE(scale))
  }
})
