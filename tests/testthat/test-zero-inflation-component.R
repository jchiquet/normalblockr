###############################################################################
## The zero-inflation component itself: the p logistic regressions of each
## variable's zero pattern on X0 (NormalBlockData$zi_fit(), surfaced on a model
## as kappa/B0/ZI_cond_mean), and the masked least-squares fit the
## zero-inflated models start from (NormalBlockData$zi_ols_fit()). Fitting
## the former is gated on the internal zero_inflation flag, and B0 is d0 x p
## whatever d0 is.
###############################################################################
set.seed(11)
ex   <- generate_normal_block_var_data(n = 30, p = 8, d = 1, q = 2)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("non zero-inflated models skip the kappa/B0 fit", {
  model <- NormalBlockVarKnownClusters$new(data, ex$parameters$C, control = NB_control(verbose = FALSE))
  priv  <- model$.__enclos_env__$private
  expect_true(is.na(priv$kappa))
  expect_true(is.na(priv$B0))
  expect_true(is.na(priv$ZI_cond_mean))

  model2 <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  priv2  <- model2$.__enclos_env__$private
  expect_true(is.na(priv2$kappa))
  expect_true(is.na(priv2$B0))
})

test_that("zero-inflated models still fit kappa/B0 as before", {
  model <- ZINormalBlockVarKnownClusters$new(data, ex$parameters$C, control = NB_control(verbose = FALSE))
  priv  <- model$.__enclos_env__$private
  expect_false(anyNA(priv$kappa))
  expect_equal(dim(priv$kappa), c(data$n, data$p))
  expect_false(is.na(priv$ZI_cond_mean))

  model2 <- ZINormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  priv2  <- model2$.__enclos_env__$private
  expect_false(anyNA(priv2$kappa))
})

test_that("zero-inflation works with more than one zero-inflation covariate", {
  ## B0 used to be built as t(sapply(...)), which is d0 x p only because
  ## sapply collapses to a vector when d0 == 1; with d0 > 1 it came out
  ## transposed and X0 %*% B0 failed with "non-conformable arguments".
  set.seed(31)
  ex <- generate_normal_block_var_data(n = 90, p = 8, d = 1, q = 2, kappa = rep(0.3, 8))
  d  <- NormalBlockData$new(ex$Y, ex$X, X0 = cbind(1, rnorm(90)))
  expect_equal(d$d0, 2)
  expect_equal(dim(d$zi_fit()$B0), c(2L, 8L))
  expect_equal(dim(d$zi_fit()$kappa), c(90L, 8L))

  fit <- normal_block(d, blocks = 2, zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE))
  expect_true(is.finite(fit$loglik))
  ## kappa costs p * d0 parameters, so a second ZI covariate is counted
  d1  <- NormalBlockData$new(ex$Y, ex$X)
  fit1 <- normal_block(d1, blocks = 2, zero_inflation = TRUE,
                       control = NB_control(verbose = FALSE))
  expect_equal(fit$nb_param - fit1$nb_param, 8L)
})

test_that("B0 keeps its d0 x p orientation in the degenerate all-nonzero case", {
  set.seed(32)
  ex <- generate_normal_block_var_data(n = 50, p = 6, d = 1, q = 2) # no zeros at all
  d  <- NormalBlockData$new(ex$Y, ex$X, X0 = cbind(1, rnorm(50)))
  expect_equal(dim(d$zi_fit()$B0), c(2L, 6L))
  expect_true(all(is.infinite(d$zi_fit()$B0)))
  expect_true(all(d$zi_fit()$kappa == 0))
})

test_that("zi_fit() is memoized and matches the per-variable logistic fits", {
  set.seed(2)
  ex <- generate_normal_block_var_data(n = 80, p = 12, d = 1, q = 2, kappa = rep(0.3, 12))
  d  <- NormalBlockData$new(ex$Y, ex$X)

  zi <- d$zi_fit()
  expect_named(zi, c("B0", "kappa", "ZI_cond_mean"))
  expect_identical(d$zi_fit(), zi)

  B0 <- t(sapply(lapply(1:d$p, function(j) {
    glm(zeros ~ 0 + ., family = binomial(link = "logit"),
        data = data.frame(zeros = d$zeros[, j], d$X0))$coefficients
  }), unlist))
  expect_equal(zi$B0, B0, ignore_attr = TRUE)
  expect_equal(zi$kappa, stats::plogis(d$X0 %*% zi$B0))
  expect_equal(zi$ZI_cond_mean,
               sum(normalblockr:::xlogy(d$zeros, zi$kappa)) +
               sum(normalblockr:::xlogy(d$zeros_bar, 1 - zi$kappa)))
})

test_that("zi_ols_fit() is memoized and matches zi_weighted_fit()", {
  set.seed(21)
  ex <- generate_normal_block_var_data(n = 70, p = 12, d = 1, q = 2, kappa = rep(0.3, 12))
  d  <- NormalBlockData$new(ex$Y, ex$X)
  fit <- d$zi_ols_fit()
  expect_named(fit, c("B", "dm1", "R"))
  expect_equal(fit, normalblockr:::zi_weighted_fit(d))
  expect_identical(d$zi_ols_fit(), fit)
  expect_true(all(fit$R[d$zeros == 1] == 0))
})

test_that("the masked fit stays finite for a variable with a single non-zero observation", {
  ## fitted exactly on one iterate, its residual sum of squares (hence dm1)
  ## would otherwise turn infinite and break the next weighted solve
  set.seed(1)
  n <- 50; p <- 5
  Y <- matrix(rexp(n * p), n, p)
  Y[, 1] <- 0
  Y[1, 1] <- 5
  data <- NormalBlockData$new(Y, matrix(1, n, 1))
  fit <- normalblockr:::zi_weighted_fit(data)
  expect_true(all(is.finite(fit$dm1)))
  expect_true(all(is.finite(fit$B)))
  model <- ZINormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  expect_no_error(model$optimize(control = NB_control(verbose = FALSE, niter = 5, threshold = -1)))
})

test_that("the masked fit tolerates a design group without any non-zero observation of a variable", {
  ## singular weighted normal equations for that variable
  set.seed(8)
  n <- 20; p <- 3
  X <- cbind(A = rep(c(1, 0), each = n / 2), B = rep(c(0, 1), each = n / 2))
  Y <- matrix(rnorm(n * p, mean = 5, sd = 1), n, p)
  Y[X[, "B"] == 1, 1] <- 0
  model <- ZINormalBlockVarKnownClusters$new(NormalBlockData$new(Y, X), matrix(1, p, 1),
                                             control = NB_control(verbose = FALSE))
  expect_no_error(model$optimize(control = list(niter = 2, threshold = -1)))
  expect_false(anyNA(model$model_par$B))
})

test_that("spherical zero-inflated models start from the pooled masked precision", {
  set.seed(5)
  exz <- generate_normal_block_var_data(n = 40, p = 10, d = 1, q = 2, kappa = rep(0.3, 10))
  dz  <- NormalBlockData$new(exz$Y, exz$X)
  ctrl <- function(nc) NB_control(verbose = FALSE, noise_covariance = nc)
  ## on the data each model is fitted on (rescaled for the variance-block family)
  pooled <- function(m) {
    fit <- m$data$zi_ols_fit()
    sum(m$data$nY) / sum(m$data$nY / fit$dm1)
  }
  init <- function(m) m$.__enclos_env__$private$optim_initialize()

  for (m in list(ZINormalBlockVarKnownClusters$new(dz, exz$parameters$C, control = ctrl("spherical")),
                 ZINormalBlockVarUnknownClusters$new(dz, 2, control = ctrl("spherical"))))
    expect_equal(init(m)$dm1, rep(pooled(m), 10))
  m <- ZINormalBlockVarKnownClusters$new(dz, exz$parameters$C, control = ctrl("diagonal"))
  expect_equal(init(m)$dm1, m$data$zi_ols_fit()$dm1)

  for (m in list(ZINormalBlockMeanKnownClusters$new(dz, exz$parameters$C, control = ctrl("spherical")),
                 ZINormalBlockMeanUnknownClusters$new(dz, 2, control = ctrl("spherical"))))
    expect_equal(diag(init(m)$Omega), rep(pooled(m), 10))
})
