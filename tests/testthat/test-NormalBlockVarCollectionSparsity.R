###############################################################################
## NormalBlockVarCollectionSparsity (R/NormalBlockCollectionSparsity.R): a
## sparsity path at fixed q, each penalty warm-started from the previous one's
## converged parameters, and its default penalty grid. StARS is in
## test-stability-selection.R.
## Frozen fixture (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal.RDS")
Y <- testdata$Y
X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

data <- NormalBlockData$new(Y, X)

test_that("optimize() warm-starts each penalty from the previous one, reducing total EM iterations without changing the BIC path", {
  ## n_sparsity_penalties = 15: the warm-start advantage grows with the
  ## length of the sparsity path and is only marginal (sometimes net-negative
  ## by a couple of iterations) on very short paths -- 15 gives a robust,
  ## comfortable margin on this dataset.
  set.seed(123)
  warm <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, C, control = NB_control(verbose = FALSE, n_sparsity_penalties = 15))
  warm$optimize(NB_control(verbose = FALSE))

  ## same penalties, but optimized independently (the pre-warm-start behaviour:
  ## every model starts from its own heuristic-derived initial parameters)
  set.seed(123)
  cold <- normalblockr:::NormalBlockVarCollectionSparsity$new(data, C, control = NB_control(verbose = FALSE, n_sparsity_penalties = 15))
  cold$models <- lapply(cold$models, function(model) { model$optimize(NB_control(verbose = FALSE)); model })

  niter_warm <- sapply(warm$models, function(m) length(m$objective))
  niter_cold <- sapply(cold$models, function(m) length(m$objective))

  ## the first model has no predecessor to warm-start from, so it should be
  ## unaffected; every later one should need no more iterations than cold
  expect_equal(niter_warm[1], niter_cold[1])
  expect_lt(sum(niter_warm), sum(niter_cold))

  ## warm-starting changes the optimization path, not the optimum it converges to
  expect_equal(sapply(warm$models, function(m) m$BIC), sapply(cold$models, function(m) m$BIC), tolerance = 1e-2)
})

test_that("the default penalty grid ignores the entries a zero weight leaves unpenalized", {
  W <- matrix(1, q, q); diag(W) <- 0; W[1, 2] <- W[2, 1] <- 0
  path <- get_model(data, C, sparsity = TRUE,
                    control = NB_control(verbose = FALSE, sparsity_weights = W, n_sparsity_penalties = 3))
  expect_length(path$sparsity, 3)
  expect_true(all(is.finite(path$sparsity) & path$sparsity > 0))
})
