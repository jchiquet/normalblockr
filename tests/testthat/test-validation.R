###############################################################################
## Exercises the stopifnot() input-validation checks throughout the package.
## These guard against silent misuse (mismatched dimensions, malformed
## clustering proposals, ...) and had zero regression coverage before -- two
## of them were previously broken (a `stopifnot(...) |> try()` that silently
## swallowed the error, and a dead `!q` condition) and went unnoticed.
set.seed(7)
ex <- generate_normal_block_data(n = 40, p = 8, d = 1, q = 2)

test_that("NormalBlockData validates its inputs", {
  expect_error(NormalBlockData$new(as.data.frame(ex$Y), ex$X), "must be matrices")
  expect_error(NormalBlockData$new(ex$Y, ex$X[1:10, , drop = FALSE]), "same number of rows")
  expect_error(NormalBlockData$new(ex$Y, ex$X, formula = "X1 + X2"), "should start with")
})

test_that("NormalBlockKnownClusters validates the clustering matrix C", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockKnownClusters$new(data, C = as.data.frame(ex$parameters$C)), "must be a matrix")
  C_empty <- ex$parameters$C
  C_empty[, 2] <- 0
  expect_error(NormalBlockKnownClusters$new(data, C = C_empty), "empty clusters")
})

test_that("clustering_init proposals are validated for unknown-clusters models", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(
    NormalBlockUnknownClusters$new(data, 2, control = NB_control(clustering_init = rep(1, 3))),
    "match the number of Y's columns"
  )
  expect_error(
    NormalBlockUnknownClusters$new(data, 3, control = NB_control(clustering_init = rep(1:2, 4))),
    "number of clusters .* must be equal to q"
  )
})

test_that("NormalBlockUnknownQ validates the list of candidate block counts", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockUnknownQ$new(data, c(2, 2, 3)), "only be present once")
  expect_error(NormalBlockUnknownQ$new(data, c(2, 3, 100)), "more blocks than there are entities")
})

test_that("NormalBlockChangingSparsity validates blocks and sparsity_penalties", {
  data <- NormalBlockData$new(ex$Y, ex$X)
  expect_error(NormalBlockChangingSparsity$new(data, blocks = c(2, 3)), "clustering matrix or a fixed number")
  expect_error(
    NormalBlockChangingSparsity$new(data, blocks = 2, control = NB_control(sparsity_penalties = c(-1, 0.1))),
    "strictly positive"
  )
})
