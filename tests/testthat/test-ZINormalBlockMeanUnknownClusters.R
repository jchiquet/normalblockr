set.seed(1234)
ex     <- generate_normal_block_mean_data(n = 100, p = 20, d = 2, q = 3)
truth  <- apply(ex$parameters$C, 1, which.max)
data0  <- NormalBlockData$new(ex$Y, ex$X)
Yzi    <- ex$Y; Yzi[runif(length(Yzi)) < 0.25] <- 0
datazi <- NormalBlockData$new(Yzi, ex$X)
ctrl   <- NB_control(verbose = FALSE)

test_that("with no zeros at all, the ZI model reduces to the plain diagonal one", {
  zi  <- normal_block(data0, 3, model = "mean", zero_inflation = TRUE, control = ctrl)
  ref <- normal_block(data0, 3, model = "mean", control = ctrl)
  expect_equal(zi$loglik, ref$loglik, tolerance = 1e-5)
  expect_equal(aricode::ARI(zi$clustering, ref$clustering), 1)
})

test_that("the VEM increases the ELBO at every iteration and recovers the clustering", {
  fit <- normal_block(datazi, 3, model = "mean", zero_inflation = TRUE, control = ctrl)
  expect_true(all(diff(fit$objective) > -1e-8))
  expect_gt(aricode::ARI(fit$clustering, truth), 0.8)
  expect_equal(dim(fit$fitted), dim(ex$Y))
  expect_true(all(fit$fitted[datazi$zeros == 1] == 0))
})

test_that("nb_param counts alpha and the zero-inflation coefficients", {
  fit <- normal_block(datazi, 3, model = "mean", zero_inflation = TRUE, control = ctrl)
  ## q*d coefficients + p variances + (q-1) alpha + p*d0 zero-inflation
  expect_equal(fit$nb_param, 3L * 2L + 20L + 2L + 20L * 1L)
})

test_that("fixed_tau leaves the clustering untouched", {
  fit <- normal_block(datazi, 3, model = "mean", zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE, fixed_tau = TRUE,
                                           clustering_init = truth))
  expect_equal(fit$clustering, truth)
})

test_that("a sparsity path is refused: it would need a full Sigma", {
  expect_error(normal_block(datazi, 3, sparsity = TRUE, model = "mean",
                            zero_inflation = TRUE, control = ctrl),
               "full Sigma")
})

test_that("a collection over a range of q runs, selects and refines", {
  coll <- normal_block(datazi, 2:5, model = "mean", zero_inflation = TRUE, control = ctrl)
  expect_s3_class(coll, "NormalBlockMeanCollectionClusters")
  expect_true(all(map_lgl(coll$models, inherits, "ZINormalBlockMeanUnknownClusters")))
  expect_equal(nrow(coll$criteria), 4L)
  expect_true(all(is.finite(coll$criteria$BIC)))
  expect_equal(coll$get_best_model("ICL")$q, 3L)

  ## deviance is non-increasing in q for nested models: a violation would
  ## measure an optimization failure, not a modelling one
  expect_true(all(diff(coll$criteria$deviance) <= 1e-6))

  coll$refine(verbose = FALSE)
  expect_true(all(is.finite(coll$criteria$BIC)))
})
