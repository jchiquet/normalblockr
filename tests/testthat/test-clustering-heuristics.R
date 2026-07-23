###############################################################################
## Regression tests for the clustering-heuristic registry (private$clustering_methods
## in NormalBlockVarBase.R): the 4 user-selectable algorithms (kmeans/ward2/sbm/
## spectral, via NB_control(clustering_init = ...)) and the ward2 fallback
## used whenever the chosen algorithm fails to produce q distinct clusters.
set.seed(123)
ex   <- generate_normal_block_var_data(n = 30, p = 10, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("every user-selectable clustering heuristic name produces a valid q-cluster init", {
  for (approx in c("kmeans", "ward2", "sbm", "spectral")) {
    model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_init = approx))
    expect_no_error(model$optimize(control = list(niter = 2, threshold = -1)))
    expect_equal(model$q, 3)
  }
})

test_that("heuristic_clustering() falls back to ward2 when the chosen method collapses clusters", {
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_init = "kmeans"))
  priv  <- model$.__enclos_env__$private

  ## Force the chosen heuristic to always return a single cluster, regardless of R
  priv$clustering_methods$kmeans <- function(R, q) rep(1, ncol(R))

  R <- matrix(rnorm(data$n * data$p), data$n, data$p)
  C <- suppressWarnings(priv$heuristic_clustering(R))

  expect_equal(dim(C), c(data$p, 3))
  expect_equal(length(unique(get_clusters(C))), 3)
})

test_that("ward2 does not error on a (near-)constant column (cor() = NA, e.g. a rare ZI variable)", {
  R <- matrix(rnorm(data$n * data$p), data$n, data$p)
  R[, 1] <- 0; R[1, 1] <- 5 # a single non-zero residual for variable 1
  priv <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))$.__enclos_env__$private
  expect_no_warning(cl <- priv$clustering_methods$ward2(R, 3))
  expect_length(cl, data$p)
  expect_equal(length(unique(cl)), 3)
})
