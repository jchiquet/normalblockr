###############################################################################
## Tests for the split/merge cluster-search machinery (NormalBlockVarBase$split(),
## $merge(), $candidates_split(), $candidates_merge()), which NormalBlockCollection
## Clusters$refine() drives. This machinery only applies to NormalBlockVarUnknownClusters
## (the clustering -- and hence M/S/C -- is variational, unlike the fixed-C
## known-clusters models), and previously had no test coverage at all.
set.seed(42)
ex   <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)

test_that("split() increases q by one and keeps every parameter conformable", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  split_model <- model$split(1)

  expect_equal(split_model$q, model$q + 1)
  expect_equal(dim(split_model$memberships), c(model$p, split_model$q))
  expect_equal(dim(split_model$var_par$M), c(model$n, split_model$q))
  expect_length(split_model$var_par$S, split_model$q)
  expect_equal(dim(split_model$model_par$Omega), c(split_model$q, split_model$q))
  expect_equal(dim(split_model$sparsity_weights), c(split_model$q, split_model$q))
})

test_that("split() does not error or scramble rows when p does not divide n evenly (M/S are individual-indexed, not variable-indexed)", {
  ## Regression test: split() used to index the n x q M (and, for ZI models,
  ## n x q S) matrices with split1/split2, masks computed over the p
  ## *variables*, not the n individuals. Whenever p happened to divide n
  ## evenly (e.g. this file's own n = 40, p = 8 fixture above), R silently
  ## recycled the mask without error -- scrambling values across unrelated
  ## individuals -- and errored outright otherwise (caught via the ZI case,
  ## where S is also n x q; see test-zi-weighted-fit.R-adjacent coverage
  ## below for that variant).
  set.seed(43)
  ex_uneven   <- generate_normal_block_var_data(n = 61, p = 17, d = 1, q = 3)
  data_uneven <- NormalBlockData$new(ex_uneven$Y, ex_uneven$X)
  model <- NormalBlockVarUnknownClusters$new(data_uneven, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  expect_no_error(split_model <- model$split(1))
  expect_equal(dim(split_model$var_par$M), c(model$n, split_model$q))
  ## the new column starts as an exact copy of its parent's
  expect_equal(split_model$var_par$M[, split_model$q], model$var_par$M[, 1])
})

test_that("split()/candidates_split()/refine() do not error on a zero-inflated model, where S is n x q rather than length q", {
  set.seed(44)
  exzi   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4, kappa = rep(0.3, 24))
  datazi <- NormalBlockData$new(exzi$Y, exzi$X, X0 = matrix(1, nrow(exzi$Y), 1))

  model <- ZINormalBlockVarUnknownClusters$new(datazi, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  expect_no_error(split_model <- model$split(1))
  expect_true(is.matrix(split_model$var_par$S))
  expect_equal(dim(split_model$var_par$S), c(model$n, split_model$q))
  expect_equal(split_model$var_par$S[, split_model$q], model$var_par$S[, 1])

  expect_no_error(model$candidates_split())

  coll <- NormalBlockVarCollectionClusters$new(datazi, 2:6, zero_inflation = TRUE, control = NB_control(verbose = FALSE))
  coll$optimize(control = list(niter = 50, threshold = 1e-4, verbose = FALSE))
  expect_no_error(coll$refine())
})

test_that("merge() decreases q by one and keeps every parameter conformable", {
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  merged_model <- model$merge(c(1, 2))

  expect_equal(merged_model$q, model$q - 1)
  expect_equal(dim(merged_model$memberships), c(model$p, merged_model$q))
  expect_equal(dim(merged_model$var_par$M), c(model$n, merged_model$q))
  expect_length(merged_model$var_par$S, merged_model$q)
  expect_equal(dim(merged_model$model_par$Omega), c(merged_model$q, merged_model$q))
})

test_that("merge() from q = 2 down to q = 1 does not error (R's drop = TRUE default would collapse the relevant matrices to vectors/a scalar)", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  merged_model <- model$merge(c(1, 2))

  expect_equal(merged_model$q, 1)
  expect_equal(dim(merged_model$memberships), c(model$p, 1))
  expect_equal(dim(merged_model$var_par$M), c(model$n, 1))
  expect_length(merged_model$var_par$S, 1)
  expect_equal(dim(merged_model$model_par$Omega), c(1, 1))
  expect_no_error(merged_model$optimize(control = list(niter = 5, threshold = -1)))
})

test_that("split() then merge() round-trips back to the original q", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  back <- model$split(1)$merge(c(1, model$q + 1))
  expect_equal(back$q, model$q)
})

test_that("split()/merge() mark the result as warm-started, so optim_initialize() reuses their Omegaq/M/S/alpha instead of a fresh heuristic estimate", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  split_model <- model$split(1)
  priv <- split_model$.__enclos_env__$private
  expect_true(priv$warm_started)
  expect_equal(priv$alpha, colMeans(priv$C))
  init <- priv$optim_initialize()
  expect_identical(init$Omega, priv$Omega)
  expect_identical(init$M, priv$M)
  expect_identical(init$S, priv$S)

  model3 <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model3$optimize(control = list(niter = 5, threshold = -1))
  merged_model <- model3$merge(c(1, 2))
  priv2 <- merged_model$.__enclos_env__$private
  expect_true(priv2$warm_started)
  expect_equal(priv2$alpha, colMeans(priv2$C))
  init2 <- priv2$optim_initialize()
  expect_identical(init2$Omega, priv2$Omega)
})

test_that("split(in_place = TRUE) mutates self instead of returning a clone", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  out <- model$split(1, in_place = TRUE)
  expect_identical(out, model)
  expect_equal(model$q, 3)
})

test_that("candidates_split() returns optimized models with one extra cluster each", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  candidates <- model$candidates_split()
  expect_true(length(candidates) > 0)
  for (cand in candidates) {
    expect_equal(cand$q, model$q + 1)
    expect_false(is.na(cand$ICL))
  }
})

test_that("candidates_merge() returns optimized models with one fewer cluster each", {
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))

  candidates <- model$candidates_merge()
  expect_true(length(candidates) > 0)
  for (cand in candidates) {
    expect_equal(cand$q, model$q - 1)
    expect_false(is.na(cand$ICL))
  }
})

test_that("candidates_merge() requires at least two clusters", {
  model <- NormalBlockVarUnknownClusters$new(data, 1, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))
  expect_error(model$candidates_merge(), "at least two clusters")
})

## Two well-separated groups of 3 variables, lumped into a single cluster:
## split(rule = "kmeans") must recover them, in both families.
two_groups <- function(family) {
  set.seed(7)
  n <- 60
  x <- rnorm(n)
  X <- cbind(1, x)
  if (family == "var") {
    W <- matrix(rnorm(2 * n), n, 2)
    signal <- W[, c(1, 1, 1, 2, 2, 2)]
  } else {
    signal <- cbind(X %*% c(0, 3), X %*% c(0, 3), X %*% c(0, 3),
                    X %*% c(0, -3), X %*% c(0, -3), X %*% c(0, -3))
  }
  Y <- signal + matrix(rnorm(6 * n, sd = 0.1), n, 6)
  NormalBlockData$new(Y, X)
}

test_that("split(rule = 'kmeans') separates two distinct groups lumped in one cluster", {
  model <- NormalBlockVarUnknownClusters$new(two_groups("var"), 1, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))
  set.seed(1)
  cl <- model$split(1, rule = "kmeans")$clustering
  expect_equal(length(unique(cl[1:3])), 1)
  expect_equal(length(unique(cl[4:6])), 1)
  expect_false(cl[1] == cl[4])

  model <- NormalBlockMeanUnknownClusters$new(two_groups("mean"), 1, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1, fixed_point_niter = 5))
  set.seed(1)
  cl <- model$split(1, rule = "kmeans")$clustering
  expect_equal(length(unique(cl[1:3])), 1)
  expect_equal(length(unique(cl[4:6])), 1)
  expect_false(cl[1] == cl[4])
})

test_that("split(rule = 'noise') keeps the previous median-of-noise-variance rule", {
  model <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))
  priv <- model$.__enclos_env__$private
  cl   <- model$clustering == 1
  var  <- 1 / priv$dm1
  expected_new <- which((var <= median(var[cl])) & cl)
  expect_setequal(which(model$split(1, rule = "noise")$clustering == model$q + 1), expected_new)
})

test_that("variance-block merge candidates are ranked by signed block correlation", {
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))
  model$optimize(control = list(niter = 5, threshold = -1))
  priv  <- model$.__enclos_env__$private
  Sigma <- matrix(c(1, -.8, .5,
                    -.8, 1, 0,
                    .5, 0, 1), 3, 3)
  priv$Omega <- solve(Sigma)
  pairs <- list(c(1, 2), c(1, 3), c(2, 3))
  expect_equal(priv$merge_score(pairs), c(-.8, .5, 0))
  ## the anticorrelated pair comes last: a merged block cannot hold
  ## negatively co-varying variables
  expect_equal(order(priv$merge_score(pairs), decreasing = TRUE), c(2, 3, 1))
})
