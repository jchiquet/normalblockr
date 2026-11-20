###############################################################################
## NormalBlockMeanCollectionSparsity (R/NormalBlockMeanCollectionSparsity.R):
## the sparsity path over the p x p Omega under noise_covariance = "full",
## get_best_model()/plot(), the n > p requirement a full Sigma imposes unless
## sparsity regularizes it, and the noise_covariance/sparsity consistency
## checks. StARS is in test-stability-selection.R.
## Frozen fixture (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

fast <- NB_control(verbose = FALSE, n_sparsity_penalties = 5)
data <- NormalBlockData$new(Y, X)
set.seed(1)
coll <- normal_block(data, blocks = q, sparsity = TRUE, model = "mean", control = fast)

test_that("mean-block sparsity path sparsifies Omega monotonically", {
  expect_s3_class(coll, "NormalBlockMeanCollectionSparsity")
  expect_length(coll$models, 5)
  expect_true(all(diff(coll$sparsity) < 0))       # decreasing penalties
  expect_true(all(diff(coll$criteria$n_edges) >= 0)) # denser as the penalty drops
  expect_equal(coll$criteria$n_edges[1], 0)       # the largest penalty empties the network
})

test_that("a known clustering is accepted along the sparsity path", {
  coll <- normal_block(data, blocks = C, sparsity = TRUE, model = "mean", control = fast)
  expect_equal(coll$q, q)
  expect_true(all(sapply(coll$models, function(m) inherits(m, "NormalBlockMeanKnownClusters"))))
})

test_that("get_best_model()/plot() work on the mean-block sparsity path", {
  for (crit in c("BIC", "EBIC", "ICL")) {
    best <- coll$get_best_model(crit)
    expect_s3_class(best, "NormalBlockMeanUnknownClusters")
    expect_equal(best[[crit]], min(coll$criteria[[crit]]))
  }
  expect_s3_class(coll$plot(c("BIC", "EBIC")), "ggplot")
  expect_message(closest <- coll$get_model(-1), "closest penalty")
  expect_equal(closest$sparsity, min(coll$sparsity))
})

test_that("only a full, unpenalized Sigma requires n > p", {
  narrow <- NormalBlockData$new(Y[1:(ncol(Y) - 2), , drop = FALSE],
                                X[1:(ncol(Y) - 2), , drop = FALSE])
  expect_error(normal_block(narrow, blocks = q, model = "mean",
                            control = NB_control(verbose = FALSE, noise_covariance = "full")),
               "need n > p", fixed = TRUE)
  for (nc in c("diagonal", "spherical"))
    expect_no_error(normal_block(narrow, blocks = q, model = "mean",
                                 control = NB_control(verbose = FALSE, noise_covariance = nc)))
  expect_no_error(normal_block(narrow, blocks = q, sparsity = TRUE, model = "mean", control = fast))
})

test_that("the mean-block objective is penalized, so $loglik is the plain ELBO", {
  model <- normal_block(data, blocks = q, sparsity = 0.01, model = "mean", control = fast)
  expect_gt(model$sparsity_term, 0)
  ## $loglik adds the penalty back to the penalized trace, as for the
  ## variance-block family: it must land above the penalized objective
  expect_equal(model$loglik, tail(model$objective, 1) + model$sparsity_term)
  expect_gt(model$loglik, tail(model$objective, 1))
})

test_that("sparsity implies a full Sigma unless one is named explicitly", {
  shape <- function(m) m$.__enclos_env__$private$res_covariance

  expect_equal(shape(normal_block(data, blocks = q, model = "mean", control = fast)),
               "diagonal")
  expect_equal(shape(normal_block(data, blocks = q, sparsity = 0.3, model = "mean",
                                  control = fast)), "full")
  expect_true(all(sapply(coll$models, shape) == "full"))

  ## asking for both explicitly is contradictory, and says so
  expect_error(
    normal_block(data, blocks = q, sparsity = 0.3, model = "mean",
                 control = NB_control(verbose = FALSE, noise_covariance = "diagonal")),
    "needs noise_covariance = 'full'", fixed = TRUE)
})
