###############################################################################
## Zero-inflated mean-block models with an unknown clustering
## (ZINormalBlockMeanUnknownClusters).
## Frozen fixture (see testdata/test_data_generation.R): data simulated on the
## fly differ across BLAS/LAPACK builds, which made the ARI check below fail
## on some platforms only.
testdata <- readRDS("testdata/testdata_normal_mean_block_zi.RDS")
truth  <- apply(testdata$C, 1, which.max)
data0  <- NormalBlockData$new(testdata$Y, testdata$X)
datazi <- NormalBlockData$new(testdata$Yzi, testdata$X)
ctrl   <- NB_control(verbose = FALSE)
set.seed(1) # the kmeans initialization draws random starts
fit    <- normal_block(datazi, 3, model = "mean", zero_inflation = TRUE, control = ctrl)

test_that("with no zeros at all, the ZI model reduces to the plain diagonal one", {
  set.seed(1)
  zi  <- normal_block(data0, 3, model = "mean", zero_inflation = TRUE, control = ctrl)
  set.seed(1)
  ref <- normal_block(data0, 3, model = "mean", control = ctrl)
  expect_equal(zi$loglik, ref$loglik, tolerance = 1e-5)
  expect_equal(aricode::ARI(zi$clustering, ref$clustering), 1)
})

test_that("the VEM increases the ELBO at every iteration and recovers the clustering", {
  expect_true(all(diff(fit$objective) > -1e-8))
  expect_gt(aricode::ARI(fit$clustering, truth), 0.8)
  expect_equal(dim(fit$fitted), dim(testdata$Y))
  expect_true(all(fit$fitted[datazi$zeros == 1] == 0))
})

test_that("nb_param counts alpha and the zero-inflation coefficients", {
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
  ## a tight threshold, so that the deviances of these nested models are
  ## compared at convergence rather than where the optimizer stopped
  set.seed(1)
  coll <- normal_block(datazi, 2:5, model = "mean", zero_inflation = TRUE,
                       control = NB_control(verbose = FALSE, threshold = 1e-8, niter = 1000))
  expect_s3_class(coll, "NormalBlockMeanCollectionClusters")
  expect_true(all(vapply(coll$models, inherits, logical(1), "ZINormalBlockMeanUnknownClusters")))
  expect_equal(nrow(coll$criteria), 4L)
  expect_true(all(is.finite(coll$criteria$BIC)))
  expect_equal(coll$get_best_model("ICL")$q, 3L)

  expect_true(all(diff(coll$criteria$deviance) <= 1e-6))

  coll$refine(verbose = FALSE)
  expect_true(all(is.finite(coll$criteria$BIC)))
})
