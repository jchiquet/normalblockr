set.seed(1234)
ex     <- generate_normal_block_mean_data(n = 100, p = 20, d = 2, q = 3)
data0  <- NormalBlockData$new(ex$Y, ex$X)
Yzi    <- ex$Y; Yzi[runif(length(Yzi)) < 0.25] <- 0
datazi <- NormalBlockData$new(Yzi, ex$X)
ctrl   <- NB_control(verbose = FALSE)

test_that("with no zeros at all, the ZI model reduces to the plain diagonal one", {
  zi  <- normal_block(data0, ex$parameters$C, model = "mean", zero_inflation = TRUE, control = ctrl)
  ref <- normal_block(data0, ex$parameters$C, model = "mean", control = ctrl)
  expect_equal(zi$loglik, ref$loglik, tolerance = 1e-5)
  expect_equal(zi$model_par$B, ref$model_par$B, tolerance = 1e-4)
})

test_that("the EM increases the log-likelihood at every iteration", {
  fit <- normal_block(datazi, ex$parameters$C, model = "mean", zero_inflation = TRUE, control = ctrl)
  expect_true(all(diff(fit$objective) > -1e-8))
  expect_equal(dim(fit$model_par$B), c(2L, 3L))
  expect_equal(dim(fit$fitted), dim(ex$Y))
})

test_that("fitted values vanish on the masked entries", {
  fit <- normal_block(datazi, ex$parameters$C, model = "mean", zero_inflation = TRUE, control = ctrl)
  expect_true(all(fit$fitted[datazi$zeros == 1] == 0))
})

test_that("the spherical variant works and a full Sigma is refused", {
  fit <- normal_block(datazi, ex$parameters$C, model = "mean", zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE, noise_covariance = "spherical"))
  expect_true(all(diff(fit$objective) > -1e-8))
  expect_error(
    normal_block(datazi, ex$parameters$C, model = "mean", zero_inflation = TRUE,
                 control = NB_control(verbose = FALSE, noise_covariance = "full")),
    "diagonal")
  expect_error(
    normal_block(datazi, ex$parameters$C, sparsity = 0.1, model = "mean",
                 zero_inflation = TRUE, control = ctrl),
    "sparsity")
})
