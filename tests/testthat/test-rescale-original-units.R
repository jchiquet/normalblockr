###############################################################################
## Regression tests for converting B/dm1/fitted back to Y's original units
## (R/NormalBlockVarBase.R: private$rescale_to_original(), B_original,
## dm1_original, and the fitted bindings in the 4 leaf classes). Y's columns
## are rescaled by NormalBlockData(scale = TRUE) (default) but not centered;
## model_par$B/model_par$dm1 stay on that internal scale (warm_start_from()
## relies on this), while B_original/dm1_original/fitted are converted back.
set.seed(1)
ex <- generate_normal_block_var_data(n = 60, p = 12, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("B_original/dm1_original undo NormalBlockData's column-wise rescaling exactly", {
  model <- NormalBlockVarKnownClusters$new(data, ex$parameters$C)
  model$optimize(NB_control(verbose = FALSE))

  s <- data$Y_scale
  expect_equal(model$B_original, model$model_par$B * matrix(s, nrow(model$model_par$B), ncol(model$model_par$B), byrow = TRUE))
  expect_equal(model$dm1_original, model$model_par$dm1 / s^2)
})

test_that("fitted is on Y's original units: X %*% B_original reproduces it exactly", {
  model <- NormalBlockVarKnownClusters$new(data, ex$parameters$C, control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_equal(model$fitted, data$X %*% model$B_original)
})

test_that("model_par$B/model_par$dm1 stay on the internal (scaled) axis, unaffected", {
  model <- NormalBlockVarKnownClusters$new(data, ex$parameters$C, control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_false(isTRUE(all.equal(model$model_par$B, model$B_original)))
})

test_that("with scale = FALSE, B_original/dm1_original/fitted equal the internal-scale versions", {
  data_noscale <- NormalBlockData$new(ex$Y, ex$X, scale = FALSE)
  model <- NormalBlockVarKnownClusters$new(data_noscale, ex$parameters$C, control = NB_control(heuristic = TRUE))
  model$optimize()
  expect_equal(model$B_original, model$model_par$B)
  expect_equal(model$dm1_original, model$model_par$dm1)
})

test_that("warm_start_from() still copies model_par$B (internal scale), not B_original", {
  model1 <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model1$optimize(NB_control(verbose = FALSE))
  model2 <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model2$warm_start_from(model1)
  expect_equal(model2$.__enclos_env__$private$B, model1$model_par$B)
})
