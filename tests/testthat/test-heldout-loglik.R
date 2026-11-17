###############################################################################
## heldout_loglik() (R/NormalBlockBase.R, family hooks in NormalBlockVarBase /
## NormalBlockMeanBase) and normal_block_cv() (R/normal_block_cv.R).
## For a known clustering, the held-out log-likelihood of the training data
## itself is the model's exact log-likelihood -- for the variance-block family
## the marginal over W, here computed through Woodbury -- up to the Jacobian
## of the column rescaling (-n sum log Y_scale).

set.seed(1)
ex <- generate_normal_block_var_data(n = 60, p = 15, d = 2, q = 3)
C  <- ex$parameters$C
ctrl <- NB_control(verbose = FALSE)

test_that("variance-block: the held-out log-likelihood of the training data is the exact one", {
  for (sc in c(TRUE, FALSE)) for (nc in c("diagonal", "spherical")) {
    d <- NormalBlockData$new(ex$Y, ex$X, scale = sc)
    m <- normal_block(d, C, control = NB_control(verbose = FALSE, noise_covariance = nc))
    expect_equal(m$heldout_loglik(d), m$loglik - d$n * sum(log(m$data$Y_scale)), tolerance = 1e-8)
  }
})

test_that("zero-inflated variance-block: same, with the zero-inflation probabilities", {
  exz <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3, kappa = rep(0.2, 15))
  d <- NormalBlockData$new(exz$Y, exz$X)
  m <- normal_block(d, exz$parameters$C, zero_inflation = TRUE, control = ctrl)
  expect_equal(m$heldout_loglik(d), m$loglik - d$n * sum(log(m$data$Y_scale)), tolerance = 1e-8)
})

test_that("mean-block (diagonal, full, zero-inflated): same", {
  em <- generate_normal_block_mean_data(n = 60, p = 15, d = 2, q = 3)
  d  <- NormalBlockData$new(em$Y, em$X)
  for (nc in c("diagonal", "full")) {
    m <- normal_block(d, em$parameters$C, model = "mean",
                      control = NB_control(verbose = FALSE, noise_covariance = nc))
    expect_equal(m$heldout_loglik(d), m$loglik, tolerance = 1e-8)
  }
  Yz <- em$Y; Yz[sample(length(Yz), 100)] <- 0
  dz <- NormalBlockData$new(Yz, em$X)
  m  <- normal_block(dz, em$parameters$C, model = "mean", zero_inflation = TRUE, control = ctrl)
  expect_equal(m$heldout_loglik(dz), m$loglik, tolerance = 1e-8)
})

test_that("per-observation values add up, and invalid inputs are refused", {
  d <- NormalBlockData$new(ex$Y, ex$X)
  m <- normal_block(d, 3, control = ctrl)
  per_obs <- m$heldout_loglik(d, per_observation = TRUE)
  expect_length(per_obs, d$n)
  expect_equal(sum(per_obs), m$heldout_loglik(d))
  expect_error(m$heldout_loglik(NormalBlockData$new(ex$Y[, 1:5], ex$X)), "same variables")
  mh <- normal_block(d, C, control = NB_control(verbose = FALSE, heuristic = TRUE))
  expect_error(mh$heldout_loglik(d), "heuristic")
})

test_that("normal_block_cv() finds the true q on well-separated variance-block data", {
  d <- NormalBlockData$new(ex$Y, ex$X)
  set.seed(2)
  cv <- normal_block_cv(d, blocks = 1:5, folds = 3)
  expect_s3_class(cv, "normal_block_cv")
  expect_equal(cv$q, 1:5)
  ## the log-likelihood plateaus beyond q = 3; the one-SE rule stops there
  expect_equal(attr(cv, "best"), 3)
  expect_gte(attr(cv, "best_max"), 3)
  expect_equal(dim(attr(cv, "per_observation")), c(60, 5))
  expect_output(print(cv), "Best: q = 3")
  expect_s3_class(plot(cv), "ggplot")
})

test_that("normal_block_cv() takes explicit folds, the mean-block family and zero-inflation", {
  em <- generate_normal_block_mean_data(n = 60, p = 12, d = 1, q = 2)
  d  <- NormalBlockData$new(em$Y, em$X)
  folds <- rep(1:2, 30)
  cv <- normal_block_cv(d, blocks = 1:3, model = "mean", folds = folds)
  expect_identical(attr(cv, "folds"), folds)
  expect_true(all(is.finite(cv$loglik)))

  exz <- generate_normal_block_var_data(n = 60, p = 12, d = 1, q = 2, kappa = rep(0.2, 12))
  set.seed(3)
  cvz <- normal_block_cv(NormalBlockData$new(exz$Y, exz$X), blocks = 2, zero_inflation = TRUE, folds = 2)
  expect_equal(nrow(cvz), 1)
  expect_true(is.finite(cvz$loglik))

  expect_error(normal_block_cv(d, blocks = 2, sparsity = TRUE), "StARS")
})
