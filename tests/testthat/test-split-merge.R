###############################################################################
## The split/merge cluster search (split(), merge(), candidates_split(),
## candidates_merge()) that NormalBlockCollectionClusters$refine() drives,
## on variance-block models with an unknown clustering.
set.seed(42)
ex   <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
short <- list(niter = 5, threshold = -1)
set.seed(1)
m2 <- NormalBlockVarUnknownClusters$new(data, 2, control = NB_control(verbose = FALSE))$optimize(short)
set.seed(1)
m3 <- NormalBlockVarUnknownClusters$new(data, 3, control = NB_control(verbose = FALSE))$optimize(short)

test_that("split() increases q by one and keeps every parameter conformable", {
  split_model <- m2$split(1)
  expect_equal(split_model$q, m2$q + 1)
  expect_equal(dim(split_model$memberships), c(m2$p, split_model$q))
  expect_equal(dim(split_model$var_par$M), c(m2$n, split_model$q))
  expect_length(split_model$var_par$S, split_model$q)
  expect_equal(dim(split_model$model_par$Omega), c(split_model$q, split_model$q))
  expect_equal(dim(split_model$sparsity_weights), c(split_model$q, split_model$q))
  expect_equal(m2$q, 2)                       # a copy, not the model itself
})

test_that("split() copies the parent's column of M, whose rows are observations, also when p does not divide n", {
  set.seed(43)
  ex_uneven   <- generate_normal_block_var_data(n = 61, p = 17, d = 1, q = 3)
  data_uneven <- NormalBlockData$new(ex_uneven$Y, ex_uneven$X)
  model <- NormalBlockVarUnknownClusters$new(data_uneven, 2, control = NB_control(verbose = FALSE))
  model$optimize(short)

  expect_no_error(split_model <- model$split(1))
  expect_equal(dim(split_model$var_par$M), c(model$n, split_model$q))
  expect_equal(split_model$var_par$M[, split_model$q], model$var_par$M[, 1])
})

test_that("split()/candidates_split()/refine() work on a zero-inflated model, where S is n x q", {
  set.seed(44)
  exzi   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4, kappa = rep(0.3, 24))
  datazi <- NormalBlockData$new(exzi$Y, exzi$X, X0 = matrix(1, nrow(exzi$Y), 1))

  model <- ZINormalBlockVarUnknownClusters$new(datazi, 2, control = NB_control(verbose = FALSE))
  model$optimize(short)

  expect_no_error(split_model <- model$split(1))
  expect_true(is.matrix(split_model$var_par$S))
  expect_equal(dim(split_model$var_par$S), c(model$n, split_model$q))
  expect_equal(split_model$var_par$S[, split_model$q], model$var_par$S[, 1])
  expect_no_error(model$candidates_split())

  set.seed(1)
  coll <- normal_block(datazi, 2:4, zero_inflation = TRUE, control = NB_control(verbose = FALSE, refine = FALSE))
  expect_no_error(coll$refine(verbose = FALSE))
})

test_that("merge() decreases q by one and keeps every parameter conformable, down to q = 1", {
  for (model in list(m3, m2)) {
    merged_model <- model$merge(c(1, 2))
    expect_equal(merged_model$q, model$q - 1)
    expect_equal(dim(merged_model$memberships), c(model$p, merged_model$q))
    expect_equal(dim(merged_model$var_par$M), c(model$n, merged_model$q))
    expect_length(merged_model$var_par$S, merged_model$q)
    expect_equal(dim(merged_model$model_par$Omega), c(merged_model$q, merged_model$q))
    expect_no_error(merged_model$optimize(short))
  }
})

test_that("merging a split cluster back restores the clustering and the variational means", {
  back <- m2$split(1)$merge(c(1, m2$q + 1))
  expect_equal(back$q, m2$q)
  expect_equal(back$clustering, m2$clustering)
  expect_equal(back$memberships, m2$memberships, tolerance = 1e-8)
  expect_equal(back$var_par$M, m2$var_par$M)
  expect_equal(back$var_par$S, m2$var_par$S)
})

test_that("split()/merge() mark the result as warm-started, so optim_initialize() reuses their state", {
  for (new_model in list(m2$split(1), m3$merge(c(1, 2)))) {
    priv <- new_model$.__enclos_env__$private
    expect_true(priv$warm_started)
    expect_equal(priv$alpha, colMeans(priv$C))
    init <- priv$optim_initialize()
    expect_identical(init$Omega, priv$Omega)
    expect_identical(init$M, priv$M)
    expect_identical(init$S, priv$S)
  }
})

test_that("split(in_place = TRUE) mutates self instead of returning a clone", {
  model <- m2$clone(deep = TRUE)
  out <- model$split(1, in_place = TRUE)
  expect_identical(out, model)
  expect_equal(model$q, 3)
})

test_that("candidates_split()/candidates_merge() return trial-optimized models with one more/fewer cluster", {
  set.seed(1)
  for (cand in m2$candidates_split()) {
    expect_equal(cand$q, m2$q + 1)
    expect_true(is.finite(cand$ICL))
  }
  candidates <- m3$candidates_merge()
  expect_length(candidates, choose(3, 2))
  for (cand in candidates) {
    expect_equal(cand$q, m3$q - 1)
    expect_true(is.finite(cand$ICL))
  }
  set.seed(1)
  m1 <- NormalBlockVarUnknownClusters$new(data, 1, control = NB_control(verbose = FALSE))$optimize(short)
  expect_error(m1$candidates_merge(), "at least two clusters")
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
  model$optimize(short)
  set.seed(1)
  cl <- model$split(1, rule = "kmeans")$clustering
  expect_equal(length(unique(cl[1:3])), 1)
  expect_equal(length(unique(cl[4:6])), 1)
  expect_false(cl[1] == cl[4])

  model <- NormalBlockMeanUnknownClusters$new(two_groups("mean"), 1, control = NB_control(verbose = FALSE))
  model$optimize(c(short, fixed_point_niter = 5))
  set.seed(1)
  cl <- model$split(1, rule = "kmeans")$clustering
  expect_equal(length(unique(cl[1:3])), 1)
  expect_equal(length(unique(cl[4:6])), 1)
  expect_false(cl[1] == cl[4])
})

test_that("split(rule = 'noise') moves the members whose noise variance is at most the median", {
  priv <- m2$.__enclos_env__$private
  cl   <- m2$clustering == 1
  var  <- 1 / priv$dm1
  expected_new <- which((var <= median(var[cl])) & cl)
  expect_setequal(which(m2$split(1, rule = "noise")$clustering == m2$q + 1), expected_new)
})

test_that("variance-block merge candidates are ranked by signed block correlation", {
  model <- m3$clone(deep = TRUE)
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
