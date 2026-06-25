###############################################################################
## Tests for the split/merge cluster-search machinery (NormalBlockBase$split(),
## $merge(), $candidates_split(), $candidates_merge()) and the SelectionNClusters
## class that drives them. This machinery only applies to NormalBlockUnknownClusters
## (the clustering -- and hence M/S/C -- is variational, unlike the fixed-C
## known-clusters models), and previously had no test coverage at all.
set.seed(42)
ex   <- generate_normal_block_data(n = 40, p = 8, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("split() increases q by one and keeps every parameter conformable", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  split_model <- model$split(1)

  expect_equal(split_model$q, model$q + 1)
  expect_equal(dim(split_model$memberships), c(model$p, split_model$q))
  expect_equal(dim(split_model$var_par$M), c(model$n, split_model$q))
  expect_length(split_model$var_par$S, split_model$q)
  expect_equal(dim(split_model$model_par$Omegaq), c(split_model$q, split_model$q))
  expect_equal(dim(split_model$sparsity_weights), c(split_model$q, split_model$q))
})

test_that("merge() decreases q by one and keeps every parameter conformable", {
  model <- NormalBlockUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  merged_model <- model$merge(c(1, 2))

  expect_equal(merged_model$q, model$q - 1)
  expect_equal(dim(merged_model$memberships), c(model$p, merged_model$q))
  expect_equal(dim(merged_model$var_par$M), c(model$n, merged_model$q))
  expect_length(merged_model$var_par$S, merged_model$q)
  expect_equal(dim(merged_model$model_par$Omegaq), c(merged_model$q, merged_model$q))
})

test_that("merge() from q = 2 down to q = 1 does not error (R's drop = TRUE default would collapse the relevant matrices to vectors/a scalar)", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  merged_model <- model$merge(c(1, 2))

  expect_equal(merged_model$q, 1)
  expect_equal(dim(merged_model$memberships), c(model$p, 1))
  expect_equal(dim(merged_model$var_par$M), c(model$n, 1))
  expect_length(merged_model$var_par$S, 1)
  expect_equal(dim(merged_model$model_par$Omegaq), c(1, 1))
  expect_no_error(merged_model$optimize(control = list(niter = 5, threshold = -1)))
})

test_that("split() then merge() round-trips back to the original q", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  back <- model$split(1)$merge(c(1, model$q + 1))
  expect_equal(back$q, model$q)
})

test_that("split()/merge() mark the result as warm-started, so EM_initialize() reuses their Omegaq/M/S/alpha instead of a fresh heuristic estimate", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  split_model <- model$split(1)
  priv <- split_model$.__enclos_env__$private
  expect_true(priv$warm_started)
  expect_equal(priv$alpha, colMeans(priv$C))
  init <- priv$EM_initialize()
  expect_identical(init$Omegaq, priv$Omegaq)
  expect_identical(init$M, priv$M)
  expect_identical(init$S, priv$S)

  model3 <- NormalBlockUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model3$optimize(control = list(niter = 5, threshold = -1))
  merged_model <- model3$merge(c(1, 2))
  priv2 <- merged_model$.__enclos_env__$private
  expect_true(priv2$warm_started)
  expect_equal(priv2$alpha, colMeans(priv2$C))
  init2 <- priv2$EM_initialize()
  expect_identical(init2$Omegaq, priv2$Omegaq)
})

test_that("split(in_place = TRUE) mutates self instead of returning a clone", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  out <- model$split(1, in_place = TRUE)
  expect_identical(out, model)
  expect_equal(model$q, 3)
})

test_that("candidates_split() returns optimized models with one extra cluster each", {
  model <- NormalBlockUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  candidates <- model$candidates_split()
  expect_true(length(candidates) > 0)
  for (cand in candidates) {
    expect_equal(cand$q, model$q + 1)
    expect_false(is.na(cand$ICL))
  }
})

test_that("candidates_merge() returns optimized models with one fewer cluster each", {
  model <- NormalBlockUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  candidates <- model$candidates_merge()
  expect_true(length(candidates) > 0)
  for (cand in candidates) {
    expect_equal(cand$q, model$q - 1)
    expect_false(is.na(cand$ICL))
  }
})

test_that("candidates_merge() requires at least two clusters", {
  model <- NormalBlockUnknownClusters$new(data, 1, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))
  expect_error(model$candidates_merge(), "at least two clusters")
})

test_that("SelectionNClusters explores the requested range and returns a consistent best model", {
  selection <- SelectionNClusters$new(data, n_clusters_range = c(2, 4),
                                      control = NB_control(verbose = FALSE))
  expect_no_error(selection$fit())

  expect_true(all(selection$ICL_explored$n_clusters >= 2))
  expect_true(all(selection$ICL_explored$n_clusters <= 4))

  best <- selection$best_model
  expect_true(isNB(best))
  expect_equal(best$ICL, min(selection$best_models$ICL))
})
