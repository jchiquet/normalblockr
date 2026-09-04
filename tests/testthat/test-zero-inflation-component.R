###############################################################################
## The zero-inflation component itself: the p logistic regressions of each
## variable's zero pattern on X0 (NormalBlockData$zi_fit(), surfaced on a model
## as kappa/B0/ZI_cond_mean). Fitting it is gated on the internal
## zero_inflation flag, and B0 is d0 x p whatever d0 is.
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
