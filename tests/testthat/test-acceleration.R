###############################################################################
## NB_control(acceleration = ): SQUAREM on top of (V)EM (src/normal_block_em_base.h)
## can be turned off to run plain (V)EM, in both families. Plain (V)EM must be
## monotone, reach the same optimum when the clustering is known, and need at
## least as many iterations.
testdata_var  <- readRDS("testdata/testdata_normal.RDS")
testdata_mean <- readRDS("testdata/testdata_normal_mean_block.RDS")
data_var  <- NormalBlockData$new(testdata_var$Y, testdata_var$X)
data_mean <- NormalBlockData$new(testdata_mean$Y, testdata_mean$X)

fit_both <- function(data, blocks, model) {
  lapply(c(acc = TRUE, plain = FALSE), function(acc) {
    set.seed(1)
    normal_block(data, blocks, model = model,
                 control = NB_control(verbose = FALSE, acceleration = acc))
  })
}

test_that("by default, SQUAREM runs for the variance-block family and a full Sigma, not for a diagonal one", {
  expect_null(NB_control()$acceleration)
  expect_false(NB_control(acceleration = FALSE)$acceleration)
  same_trace <- function(a, b) isTRUE(all.equal(a$objective, b$objective))
  fit <- function(data, blocks, model, acc, nc = NULL) {
    set.seed(1)
    normal_block(data, blocks, model = model,
                 control = NB_control(verbose = FALSE, acceleration = acc, noise_covariance = nc))
  }
  ## variance-block: the default is the accelerated fit
  C_var <- testdata_var$parameters$C
  expect_true(same_trace(fit(data_var, C_var, "var", NULL), fit(data_var, C_var, "var", TRUE)))
  ## mean-block: the rule itself (on these data the EM converges in two or
  ## three iterations, with or without SQUAREM, so traces cannot tell)
  C_mean <- testdata_mean$parameters$C
  uses <- function(nc, acc) {
    m <- NormalBlockMeanKnownClusters$new(data_mean, C_mean,
           control = NB_control(verbose = FALSE, noise_covariance = nc))
    m$.__enclos_env__$private$use_acceleration(NB_control(acceleration = acc))
  }
  expect_false(uses("diagonal", NULL))
  expect_false(uses("spherical", NULL))
  expect_true(uses("full", NULL))
  expect_true(uses("diagonal", TRUE))
  expect_false(uses("full", FALSE))
})

test_that("variance-block, known clusters: plain EM is monotone and reaches the accelerated optimum in more iterations", {
  fits <- fit_both(data_var, testdata_var$parameters$C, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_equal(fits$plain$loglik, fits$acc$loglik, tolerance = 1e-6)
  expect_gt(fits$plain$criteria$niter, fits$acc$criteria$niter)
})

test_that("variance-block, unknown clusters: the switch reaches the VEM", {
  fits <- fit_both(data_var, 3, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_gt(fits$plain$criteria$niter, fits$acc$criteria$niter)
})

test_that("mean-block: plain (V)EM is monotone and reaches the same optimum", {
  known <- fit_both(data_mean, testdata_mean$parameters$C, "mean")
  expect_true(all(diff(known$plain$objective) > -1e-8))
  expect_equal(known$plain$loglik, known$acc$loglik, tolerance = 1e-6)

  unknown <- fit_both(data_mean, 3, "mean")
  expect_true(all(diff(unknown$plain$objective) > -1e-8))
  expect_equal(unknown$plain$loglik, unknown$acc$loglik, tolerance = 1e-6)
})

test_that("zero-inflated models accept the switch", {
  zi <- readRDS("testdata/testdata_normal_zi.RDS")
  data_zi <- NormalBlockData$new(zi$Y, zi$X)
  fits <- fit_both(data_zi, zi$parameters$C, "var")
  expect_true(all(diff(fits$plain$objective) > -1e-8))
  expect_equal(fits$plain$loglik, fits$acc$loglik, tolerance = 1e-5)
})

test_that("mean-block, full Sigma: from the same start, SQUAREM reaches at least the plain recursion's ELBO", {
  Y <- testdata_mean$Y; C <- testdata_mean$parameters$C; q <- ncol(C)
  data <- data_mean
  ctrl <- NB_control(verbose = FALSE, clustering_init = get_clusters(C),
                     noise_covariance = "full")
  init <- NormalBlockMeanUnknownClusters$new(data, q, control = ctrl
            )$.__enclos_env__$private$optim_initialize()

  fit <- function(accelerate) {
    NormalBlockMeanUnknownClusters_fit(
      Y = data$Y, X = data$X, B0 = init$B, Omega0 = init$Omega, tau0 = init$tau,
      sparsity = 0, sparsity_weights = matrix(1, ncol(Y), ncol(Y)) - diag(ncol(Y)),
      noise_covariance = "full", fixed_point_niter = 5, niter = 50, threshold = 1e-6,
      accelerate = accelerate)
  }
  plain <- fit(FALSE)
  accel <- fit(TRUE)

  expect_gte(tail(accel$objective, 1), tail(plain$objective, 1) - 1e-6)
  expect_true(all(diff(accel$objective) >= -1e-6)) # still monotone
})
