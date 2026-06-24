###############################################################################
###############################################################################
## Use pre-save testdata (seed are hard to handle in testhat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

niter     <- 6
threshold <- -1 # never trigger early stopping: forces exactly `niter` iterations

###############################################################################
###############################################################################
## These tests check that the Rcpp/RcppArmadillo (V)EM core
## (NormalBlockKnownClusters_fit / NormalBlockUnknownClusters_fit, src/exports.cpp)
## reproduces *exactly* (up to numerical precision) the EM/VEM recursion
## already implemented in R (NormalBlockKnownClusters / NormalBlockUnknownClusters), starting from the
## same initial parameters, both unpenalized (sparsity = 0, plain inversion)
## and penalized (sparsity > 0, graphical lasso via glassoFast, called back
## from C++, see src/omega_estimation.h). The R6 heuristic initialization
## (private method EM_initialize) is reached into via R6's `.__enclos_env__`
## so that both implementations start from the very same point.

test_that("NormalBlockKnownClusters_fit matches NormalBlockKnownClusters (diagonal/spherical, unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)

  for (nc in c("diagonal", "spherical")) {
    for (sparsity in c(0, 0.05)) {
      model <- NormalBlockKnownClusters$new(data, C, sparsity = sparsity,
                                    control = NB_control(noise_covariance = nc, verbose = FALSE))
      init <- model$.__enclos_env__$private$EM_initialize()
      model$optimize(control = list(niter = niter, threshold = threshold))

      res <- NormalBlockKnownClusters_fit(Y = data$Y, X = data$X, C = C,
                                    B0 = init$B, dm1_0 = init$dm1, Omegaq0 = init$Omegaq,
                                    sparsity = sparsity, sparsity_weights = model$sparsity_weights,
                                    noise_covariance = nc, niter = niter, threshold = threshold)

      expect_equal(res$B,      model$model_par$B,      tolerance = 1e-8)
      expect_equal(res$dm1,    model$model_par$dm1,    tolerance = 1e-8)
      expect_equal(res$Omegaq, model$model_par$Omegaq, tolerance = 1e-8)
      expect_equal(res$gamma,  model$posterior_par$gamma, tolerance = 1e-8)
      expect_equal(res$mu,     model$posterior_par$mu,    tolerance = 1e-8)
      expect_equal(res$objective[-1], model$objective, tolerance = 1e-8)
    }
  }
})

test_that("NormalBlockUnknownClusters_fit matches NormalBlockUnknownClusters (diagonal/spherical, unpenalized/sparse)", {
  data <- NormalBlockData$new(Y, X)

  for (nc in c("diagonal", "spherical")) {
    for (sparsity in c(0, 0.05)) {
      model <- NormalBlockUnknownClusters$new(data, q, sparsity = sparsity,
                              control = NB_control(noise_covariance = nc, verbose = FALSE))
      init <- model$.__enclos_env__$private$EM_initialize()
      model$optimize(control = list(niter = niter, threshold = threshold))

      res <- NormalBlockUnknownClusters_fit(Y = data$Y, X = data$X,
                                      B0 = init$B, dm1_0 = init$dm1, Omegaq0 = init$Omegaq,
                                      C0 = init$C, alpha0 = init$alpha, M0 = init$M, S0 = init$S,
                                      sparsity = sparsity, sparsity_weights = model$sparsity_weights,
                                      noise_covariance = nc, fixed_tau = FALSE,
                                      niter = niter, threshold = threshold)

      expect_equal(res$B,      model$model_par$B,      tolerance = 1e-8)
      expect_equal(res$dm1,    model$model_par$dm1,    tolerance = 1e-8)
      expect_equal(res$Omegaq, model$model_par$Omegaq, tolerance = 1e-8)
      expect_equal(res$C,      model$var_par$tau,      tolerance = 1e-8)
      expect_equal(res$alpha,  model$model_par$alpha,  tolerance = 1e-8)
      expect_equal(res$M,      model$var_par$M,        tolerance = 1e-8)
      expect_equal(res$S,      model$var_par$S,        tolerance = 1e-8)
      expect_equal(res$objective[-1], model$objective, tolerance = 1e-8)
    }
  }
})
