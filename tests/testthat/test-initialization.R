###############################################################################
## Regression test: NormalBlockVarBase$initialize() used to fit a zero-inflation
## logistic regression (kappa/B0) unconditionally for every model, including
## plain (non zero-inflated) ones that never read kappa/B0/ZI_cond_mean
## afterwards. This is now gated on the (internal) zero_inflation flag set by
## the ZI subclasses, so non-ZI models should skip that fit entirely.
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
