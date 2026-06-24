###############################################################################
###############################################################################
## Use pre-saved zero-inflated testdata (seed are hard to handle in testthat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

niter     <- 5
threshold <- -1 # never trigger early stopping: forces exactly `niter` iterations

###############################################################################
###############################################################################
## These tests check that the Rcpp/RcppArmadillo zero-inflated (V)EM core
## (ZINormalBlockKnownClusters_fit / ZINormalBlockUnknownClusters_fit, src/exports.cpp)
## reproduces *closely* (up to numerical precision) the EM/VEM recursion
## already implemented in R (ZINormalBlockKnownClusters / ZINormalBlockUnknownClusters), starting from
## the same initial parameters, both unpenalized (sparsity = 0) and penalized
## (sparsity > 0, graphical lasso via glassoFast). The B-update (and, for
## unknown clusters, the M-update) is exactly quadratic despite the
## zero-inflation mask, so the C++ side solves it directly (see
## src/zi_closed_form_solvers.h) instead of going through an iterative
## optimizer; the R reference still uses nloptr/L-BFGS for these subproblems
## (R/ZINormalBlockKnownClusters.R, R/ZINormalBlockUnknownClusters.R), which only
## converges to its own optimizer tolerance -- hence the looser tolerance
## here than in the non zero-inflated tests (test-cpp-normal-block.R), which
## only involve exact linear algebra on both sides. The R6 heuristic
## initialization (private method EM_initialize) is reached into via R6's
## `.__enclos_env__` so that both implementations start from the very same
## point.

test_that("ZINormalBlockKnownClusters_fit matches ZINormalBlockKnownClusters (diagonal/spherical, unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)

  for (nc in c("diagonal", "spherical")) {
    for (sparsity in c(0, 0.05)) {
      model <- ZINormalBlockKnownClusters$new(data, C, sparsity = sparsity,
                                      control = NB_control(noise_covariance = nc, verbose = FALSE))
      init <- model$.__enclos_env__$private$EM_initialize()
      model$optimize(control = list(niter = niter, threshold = threshold))

      zi_cond_mean <- model$.__enclos_env__$private$ZI_cond_mean

      res <- ZINormalBlockKnownClusters_fit(Y = data$Y, X = data$X,
                                    zeros_bar = data$zeros_bar, zi_cond_mean = zi_cond_mean, C = C,
                                    B0 = init$B, dm1_0 = init$dm1, Omegaq0 = init$Omegaq,
                                    sparsity = sparsity, sparsity_weights = model$sparsity_weights,
                                    noise_covariance = nc, niter = niter, threshold = threshold)

      gamma_r <- model$posterior_par$gamma
      gamma_cpp <- lapply(seq_len(data$n), function(i) res$gamma[, , i])

      expect_equal(res$B,      model$model_par$B,      tolerance = 1e-3)
      expect_equal(res$dm1,    model$model_par$dm1,    tolerance = 1e-3)
      expect_equal(res$Omegaq, model$model_par$Omegaq, tolerance = 1e-3)
      expect_equal(res$mu,     model$posterior_par$mu, tolerance = 1e-3)
      expect_equal(gamma_cpp,  gamma_r,                tolerance = 1e-3)
      expect_equal(res$objective[-1], model$objective, tolerance = 1e-3)
    }
  }
})

test_that("ZINormalBlockUnknownClusters_fit matches ZINormalBlockUnknownClusters (diagonal/spherical, unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)

  for (nc in c("diagonal", "spherical")) {
    for (sparsity in c(0, 0.05)) {
      model <- ZINormalBlockUnknownClusters$new(data, q, sparsity = sparsity,
                                 control = NB_control(noise_covariance = nc, verbose = FALSE))
      init <- model$.__enclos_env__$private$EM_initialize()
      model$optimize(control = list(niter = niter, threshold = threshold))

      zi_cond_mean <- model$.__enclos_env__$private$ZI_cond_mean

      res <- ZINormalBlockUnknownClusters_fit(Y = data$Y, X = data$X,
                                      zeros_bar = data$zeros_bar, zi_cond_mean = zi_cond_mean,
                                      B0 = init$B, dm1_0 = init$dm1, Omegaq0 = init$Omegaq,
                                      C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
                                      sparsity = sparsity, sparsity_weights = model$sparsity_weights,
                                      noise_covariance = nc, fixed_tau = FALSE,
                                      niter = niter, threshold = threshold)

      expect_equal(res$B,      model$model_par$B,      tolerance = 1e-3)
      expect_equal(res$dm1,    model$model_par$dm1,    tolerance = 1e-3)
      expect_equal(res$Omegaq, model$model_par$Omegaq, tolerance = 1e-3)
      expect_equal(res$C,      model$var_par$tau,      tolerance = 1e-3)
      expect_equal(res$alpha,  model$model_par$alpha,  tolerance = 1e-3)
      expect_equal(res$M,      model$var_par$M,        tolerance = 1e-3)
      expect_equal(res$S,      model$var_par$S,        tolerance = 1e-3)
      expect_equal(res$objective[-1], model$objective, tolerance = 1e-3)
    }
  }
})
