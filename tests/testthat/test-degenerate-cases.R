###############################################################################
## Regression tests for two degenerate-but-realistic situations found while
## debugging an analysis script on real (heavily zero-inflated, many small
## design groups) data: an empty cluster in the inferred clustering, and a
## design group with zero non-zero-inflated observations for some variable.
set.seed(7)
ex <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 3)

test_that("plot_network() does not error when the inferred clustering has an empty cluster", {
  data  <- NormalBlockData$new(ex$Y, ex$X)
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 2, threshold = -1))

  ## Force an empty cluster: column 3 of the (soft) membership matrix never
  ## wins the argmax, regardless of the actual EM outcome above -- this is
  ## what naturally happens on real data when q is larger than the number of
  ## clusters the VEM actually uses.
  p <- model$p
  C_forced <- cbind(matrix(c(0.6, 0.4), p, 2, byrow = TRUE), 0)
  model$.__enclos_env__$private$C <- C_forced

  expect_length(tabulate(model$clustering, nbins = model$q), model$q)
  expect_no_error(model$plot_network(plot = FALSE))
})

test_that("ZI heuristic initialization tolerates a fully zero-inflated design group", {
  ## One-hot design with two groups; variable 1 is entirely zero (hence
  ## zero-inflated) within group B. This makes the per-column weighted normal
  ## equations exactly singular in the B-update heuristic -- regression test
  ## for the ginv() fix (was previously a hard `solve()` error).
  n <- 20
  p <- 3
  X <- cbind(A = rep(c(1, 0), each = n / 2), B = rep(c(0, 1), each = n / 2))
  Y <- matrix(rnorm(n * p, mean = 5, sd = 1), n, p)
  Y[X[, "B"] == 1, 1] <- 0

  data  <- NormalBlockData$new(Y, X)
  model <- ZINormalBlockVarKnownClusters$new(data, matrix(1, p, 1), control = NB_control(verbose = FALSE))

  expect_no_error(model$optimize(control = list(niter = 2, threshold = -1)))
  expect_false(anyNA(model$model_par$B))
})
