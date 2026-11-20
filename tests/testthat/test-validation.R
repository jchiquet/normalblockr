###############################################################################
## Exercises the stopifnot() input-validation checks throughout the package.
## These guard against silent misuse (mismatched dimensions, malformed
## clustering proposals, ...) and had zero regression coverage before -- two
## of them were previously broken (a `stopifnot(...) |> try()` that silently
## swallowed the error, and a dead `!q` condition) and went unnoticed.
set.seed(7)
ex <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 2)

test_that("NormalBlockData validates its inputs", {
  expect_error(NormalBlockData$new(as.data.frame(ex$Y), ex$X), "must be matrices")
  expect_error(NormalBlockData$new(ex$Y, ex$X[1:10, , drop = FALSE]), "same number of rows")
  expect_error(NormalBlockData$new(ex$Y, ex$X, formula = "X1 + X2"), "should start with")
})

test_that("NormalBlockVarKnownClusters validates the clustering matrix C", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockVarKnownClusters$new(data, C = as.data.frame(ex$parameters$C)), "must be a matrix")
  C_empty <- ex$parameters$C
  C_empty[, 2] <- 0
  expect_error(NormalBlockVarKnownClusters$new(data, C = C_empty), "empty clusters")
})

test_that("clustering_init proposals are validated for unknown-clusters models", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(
    NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(clustering_init = rep(1, 3))),
    "match the number of Y's columns"
  )
  expect_error(
    NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(clustering_init = rep(1:2, 4))),
    "number of clusters .* must be equal to q"
  )
})

test_that("NormalBlockVarCollectionClusters validates the list of candidate block counts", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockVarCollectionClusters$new(data, c(2, 2, 3)), "only be present once")
  expect_error(NormalBlockVarCollectionClusters$new(data, c(2, 3, 100)), "more blocks than there are entities")
})

test_that("NormalBlockVarCollectionSparsity validates blocks and sparsity_penalties", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockVarCollectionSparsity$new(data, blocks = c(2, 3)), "clustering matrix or a fixed number")
  expect_error(
    NormalBlockVarCollectionSparsity$new(data, blocks = 2, control = NB_control(sparsity_penalties = c(-1, 0.1))),
    "strictly positive"
  )
})

test_that("known-cluster models check that C has one row per variable", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  C_short <- ex$parameters$C[-1, ]
  for (cls in list(NormalBlockVarKnownClusters, ZINormalBlockVarKnownClusters,
                   NormalBlockMeanKnownClusters, ZINormalBlockMeanKnownClusters))
    expect_error(cls$new(data, C = C_short), "one row per variable")
})

test_that("normal_block() and NB_control() refuse malformed arguments", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(normal_block(data, blocks = c(1, 2, 1, 2, 1, 2, 1, 2)), "vector of cluster labels")
  expect_error(NB_control(noise_covariance = "block"), "noise_covariance must be one of")
  expect_error(NB_control(clustering_init = "random"), "must name a known heuristic")
  expect_error(normal_block(data, blocks = 2, model = "var",
                            control = NB_control(verbose = FALSE, noise_covariance = "full")),
               "only applies to mean-block models")
  expect_error(normal_block(data, blocks = 2:3, control = NB_control(clustering_init = list(rep(1:2, 4)))))
})

test_that("an initial clustering given as probabilities may have a nearly empty cluster", {
  set.seed(2)
  ex <- generate_normal_block_var_data(n = 40, p = 12, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  tau <- ex$parameters$C * (1 - 2e-3) + 1e-3
  tau[, 3] <- 1e-4; tau <- tau / rowSums(tau)   # cluster 3 carries less than one variable
  expect_lt(sum(tau[, 3]), 1)
  m <- NormalBlockVarUnknownClusters$new(d, 3, control = NB_control(clustering_init = tau, fixed_tau = TRUE, verbose = FALSE))
  expect_no_error(m$optimize(NB_control(fixed_tau = TRUE, verbose = FALSE)))
  ## an indicator matrix still needs a member in every cluster
  C <- ex$parameters$C; C[, 3] <- 0
  expect_error(NormalBlockVarUnknownClusters$new(d, 3, control = NB_control(clustering_init = C)), "empty clusters")
})
