###############################################################################
## Work that a collection over q used to repeat once per model, even though it
## does not depend on q: the OLS fit, the zero-inflation logistic regressions,
## and the q-independent part of a clustering heuristic.
##
## These are pure removals of redundant computation -- the fits they feed must
## come out unchanged, which is what most of these tests check.
###############################################################################

test_that("ols_fit() is memoized and matches a fresh computation", {
  set.seed(1)
  ex <- generate_normal_block_var_data(n = 60, p = 20, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)

  first <- d$ols_fit()
  expect_named(first, c("B", "R", "Sigma"))
  expect_equal(first$B, d$XtXm1 %*% d$XtY)
  expect_equal(first$R, d$Y - d$X %*% first$B)
  expect_equal(first$Sigma, cov(first$R))

  ## same object back, not a recomputation
  expect_identical(d$ols_fit(), first)
})

test_that("zi_fit() is memoized and matches the per-variable logistic fits", {
  set.seed(2)
  ex <- generate_normal_block_var_data(n = 80, p = 12, d = 1, q = 2, kappa = rep(0.3, 12))
  d  <- NormalBlockData$new(ex$Y, ex$X)

  zi <- d$zi_fit()
  expect_named(zi, c("B0", "kappa", "ZI_cond_mean"))
  expect_identical(d$zi_fit(), zi)

  ## the same p independent logistic regressions, spelled out
  B0 <- t(sapply(lapply(1:d$p, function(j) {
    glm(zeros ~ 0 + ., family = binomial(link = "logit"),
        data = data.frame(zeros = d$zeros[, j], d$X0))$coefficients
  }), unlist))
  ## values only: B0's dimnames are a leaked deparse of whatever expression
  ## built the data.frame, which is a separate (pre-existing) wart
  expect_equal(zi$B0, B0, ignore_attr = TRUE)
  expect_equal(zi$ZI_cond_mean,
               sum(normalblockr:::xlogy(d$zeros, zi$kappa)) +
               sum(normalblockr:::xlogy(d$zeros_bar, 1 - zi$kappa)))
})

test_that("every model in a collection shares one data object, hence one fit", {
  set.seed(3)
  ex <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 2, kappa = rep(0.3, 15))
  d  <- NormalBlockData$new(ex$Y, ex$X)
  coll <- normal_block(d, blocks = 2:4, zero_inflation = TRUE,
                       control = NB_control(verbose = FALSE))

  ## R6 reference semantics: the models point at the caller's object, so the
  ## memoized zi_fit()/ols_fit() are computed once for the whole collection
  for (m in coll$models) expect_identical(m$data, d)
})

test_that("a q-range collection shares the q-independent part of its heuristic", {
  set.seed(4)
  ex <- generate_normal_block_var_data(n = 80, p = 25, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  R  <- normalblockr:::ols_residuals(d)
  q_list <- 2:5

  ## ward2: one tree, cut at each q -- identical to cutting a per-q tree
  path <- normalblockr:::clustering_path_for_collection(R, q_list, "ward2")
  expect_equal(length(path), length(q_list))
  for (i in seq_along(q_list)) {
    expect_equal(path[[i]], cutree(normalblockr:::ward2_tree(R), q_list[i]))
  }

  ## nothing is shareable for kmeans (it depends on q entirely) or for the
  ## model's own best_of_inits search
  expect_null(normalblockr:::clustering_path_for_collection(R, q_list, "kmeans"))
  expect_null(normalblockr:::clustering_path_for_collection(R, q_list, "best_of_inits"))
  expect_null(normalblockr:::clustering_path_for_collection(R, q_list, c("ward2", "kmeans")))
})

test_that("sharing the heuristic leaves the fitted collection unchanged", {
  set.seed(5)
  ex <- generate_normal_block_var_data(n = 80, p = 25, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  q_list <- 2:5

  shared <- normal_block(d, blocks = q_list, control = NB_control(verbose = FALSE))

  ## the same models, each cold-started from its own independently-built tree
  per_model <- lapply(q_list, function(q) {
    cl <- cutree(normalblockr:::ward2_tree(normalblockr:::ols_residuals(d)), q)
    normal_block(d, blocks = q, control = NB_control(verbose = FALSE, clustering_init = cl))
  })

  expect_equal(sapply(shared$models, function(m) m$loglik),
               sapply(per_model, function(m) m$loglik))
  expect_equal(lapply(shared$models, function(m) m$clustering),
               lapply(per_model, function(m) m$clustering))
})

test_that("a heuristic that collapses below q falls back to the model's own", {
  ## clustering_path_for_collection() must never hand down a clustering with
  ## the wrong number of groups: NormalBlockBase would reject it outright,
  ## where heuristic_clustering() has a documented fallback. Asking for more
  ## groups than there are distinct variables forces the degenerate case.
  set.seed(6)
  R <- matrix(rnorm(40), 20, 2)[, c(1, 1, 2, 2)] # only two distinct columns
  path <- normalblockr:::clustering_path_for_collection(R, c(2, 4), "ward2")
  ## whatever it returns per q is either exactly q groups, or the method name
  for (i in seq_along(path)) {
    cl <- path[[i]]
    expect_true(identical(cl, "ward2") || length(unique(cl)) == c(2, 4)[i])
  }
})
