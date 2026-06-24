###############################################################################
## Regression tests for the clustering-heuristic registry (private$clustering_methods
## in NormalBlockBase.R): the 4 user-selectable algorithms (kmeans/ward2/sbm/
## kmeansvar, via NB_control(clustering_approx = ...)) and the hidden hclustvar
## fallback used whenever the chosen algorithm fails to produce q distinct
## clusters, previously dispersed across 4 separate private methods.
set.seed(123)
ex   <- generate_normal_block_data(n = 30, p = 10, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("every user-selectable clustering_approx produces a valid q-cluster init", {
  for (approx in c("kmeans", "ward2", "sbm", "kmeansvar")) {
    model <- NormalBlockUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_approx = approx))
    expect_no_error(model$optimize(control = list(niter = 2, threshold = -1)))
    expect_equal(model$q, 3)
  }
})

test_that("heuristic_clustering() falls back to hclustvar when the chosen method collapses clusters", {
  model <- NormalBlockUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE, clustering_approx = "kmeans"))
  priv  <- model$.__enclos_env__$private

  ## Force the chosen heuristic to always return a single cluster, regardless of R
  priv$clustering_methods$kmeans <- function(R, q) rep(1, ncol(R))

  R <- matrix(rnorm(data$n * data$p), data$n, data$p)
  C <- suppressWarnings(priv$heuristic_clustering(R))

  expect_equal(dim(C), c(data$p, 3))
  expect_equal(length(unique(get_clusters(C))), 3)
})
