###############################################################################
## NormalBlockMeanCollectionSparsity/ClustersSparsity (R/NormalBlockMeanCollection*.R):
## the sparsity path over the p x p Omega under noise_covariance = "full",
## its interaction with q (crossed sparsity/clusters), get_best_model()/plot(),
## the n > p requirement a full Sigma imposes unless sparsity regularizes it,
## and the noise_covariance/sparsity consistency checks.
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal_mean_block.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

fast <- NB_control(verbose = FALSE, n_sparsity_penalties = 5, niter = 20)
data <- NormalBlockData$new(Y, X)

test_that("mean-block sparsity path sparsifies Omega monotonically", {
  coll <- normal_block(data, blocks = q, sparsity = TRUE, model = "mean", control = fast)

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
  coll <- normal_block(data, blocks = q, sparsity = TRUE, model = "mean", control = fast)
  for (crit in c("BIC", "EBIC", "ICL")) {
    best <- coll$get_best_model(crit)
    expect_s3_class(best, "NormalBlockMeanUnknownClusters")
    expect_equal(best[[crit]], min(coll$criteria[[crit]]))
  }
  expect_s3_class(coll$plot(c("BIC", "EBIC")), "ggplot")
  expect_message(coll$get_model(-1), "closest penalty")
})

test_that("mean-block models require n > p unless the penalty regularizes Omega", {
  narrow <- NormalBlockData$new(Y[1:(ncol(Y) - 2), , drop = FALSE],
                                X[1:(ncol(Y) - 2), , drop = FALSE])
  full <- NB_control(verbose = FALSE, n_sparsity_penalties = 5, niter = 20,
                     noise_covariance = "full")
  expect_error(normal_block(narrow, blocks = q, model = "mean", control = full),
               "need n > p", fixed = TRUE)
  expect_no_error(
    normal_block(narrow, blocks = q, sparsity = TRUE, model = "mean", control = fast))
})

test_that("the mean-block objective is penalized, so $loglik is the plain ELBO", {
  model <- normal_block(data, blocks = q, sparsity = 0.01, model = "mean", control = fast)
  expect_gt(model$sparsity_term, 0)
  ## $loglik adds the penalty back to the penalized trace, as for the
  ## variance-block family: it must land above the penalized objective
  expect_equal(model$loglik, tail(model$objective, 1) + model$sparsity_term)
  expect_gt(model$loglik, tail(model$objective, 1))
})

test_that("q and sparsity can be crossed for mean-block models", {
  coll <- normal_block(data, blocks = 2:3, sparsity = TRUE, model = "mean",
                       control = NB_control(verbose = FALSE, n_sparsity_penalties = 3, niter = 15))
  expect_s3_class(coll, "NormalBlockMeanCollectionClustersSparsity")
  expect_equal(nrow(coll$criteria), 6)
  expect_s3_class(coll$get_model(2), "NormalBlockMeanCollectionSparsity")
  expect_s3_class(coll$get_best_model("BIC"), "NormalBlockMeanBase")
  expect_s3_class(coll$plot("BIC"), "ggplot")
})

test_that("sparsity implies a full Sigma unless one is named explicitly", {
  shape <- function(m) m$.__enclos_env__$private$res_covariance

  expect_equal(shape(normal_block(data, blocks = q, model = "mean", control = fast)),
               "diagonal")
  expect_equal(shape(normal_block(data, blocks = q, sparsity = 0.3, model = "mean",
                                  control = fast)), "full")
  coll <- normal_block(data, blocks = q, sparsity = TRUE, model = "mean", control = fast)
  expect_true(all(sapply(coll$models, shape) == "full"))

  ## asking for both explicitly is contradictory, and says so
  expect_error(
    normal_block(data, blocks = q, sparsity = 0.3, model = "mean",
                 control = NB_control(verbose = FALSE, noise_covariance = "diagonal")),
    "needs noise_covariance = 'full'", fixed = TRUE)
})

## StARS, shared with the variance-block family (NormalBlockCollectionSparsity):
## the network relates the p variables, and the subsample paths hold the
## clustering fixed (known, or the MAP partition of the least penalized model).
test_that("StARS works on mean-block sparsity paths, known or unknown clusters", {
  set.seed(2)
  ex <- generate_normal_block_mean_data(n = 80, p = 12, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  ctrl <- NB_control(verbose = FALSE, n_sparsity_penalties = 4)
  for (blocks in list(ex$parameters$C, 3)) {
    path <- normal_block(d, blocks = blocks, sparsity = TRUE, model = "mean", control = ctrl)
    expect_true(all(is.na(path$stability)))
    set.seed(3)
    path$stability_selection(n_subsamples = 3)
    sp <- path$stability_path
    expect_equal(nrow(sp), choose(12, 2) * 4)
    expect_true(all(sp$Prob >= 0 & sp$Prob <= 1))
    expect_equal(sort(unique(sp$Penalty)), sort(path$sparsity))
    expect_match(sp$Edge[1], "^Y1\\|Y2$")
    expect_true(all(path$stability >= 0 & path$stability <= 1))
    expect_equal(path$criteria$stability, path$stability)
    best <- path$get_best_model("StARS", stability = 0)
    expect_true(best$sparsity %in% path$sparsity)
    expect_s3_class(best, "NormalBlockMeanBase")
  }
})

test_that("the subsample paths reuse the collection's penalties, so the stability path is labelled right", {
  set.seed(4)
  ex <- generate_normal_block_mean_data(n = 60, p = 10, d = 1, q = 2)
  path <- normal_block(NormalBlockData$new(ex$Y, ex$X), blocks = ex$parameters$C, sparsity = TRUE,
                       model = "mean", control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
  ## a subsample would otherwise derive its own penalty grid from its own
  ## covariance, and the probabilities would be filed under the wrong penalty
  ctrl <- path$control; ctrl$sparsity_penalties <- path$sparsity
  rows <- 1:40
  sub_data <- NormalBlockData$new(ex$Y[rows, ], ex$X[rows, , drop = FALSE])
  sub_path <- path$.__enclos_env__$private$subsample_path(sub_data, ctrl)
  expect_equal(sub_path$sparsity, path$sparsity)
})
