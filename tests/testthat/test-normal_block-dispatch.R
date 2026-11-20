###############################################################################
## get_model()/normal_block() dispatch: every combination of the family,
## `blocks`, `sparsity` and `zero_inflation` reaches the class it is meant to
## build. The quality of the fits is not the point (niter is tiny).
set.seed(7)
ex   <- generate_normal_block_var_data(n = 40, p = 8, d = 1, q = 2)
data <- NormalBlockData$new(ex$Y, ex$X)
C    <- ex$parameters$C
fast <- NB_control(niter = 2, threshold = -1, verbose = FALSE, sparsity_penalties = c(0.1, 0.05))

test_that("normal_block() builds the class matching its arguments, in both families", {
  cases <- expand.grid(model = c("var", "mean"), blocks = c("C", "q", "q_list"),
                       sparsity = c("none", "fixed", "path"), zi = c(FALSE, TRUE),
                       stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cases))) {
    k <- cases[i, ]
    blocks <- switch(k$blocks, C = C, q = 2, q_list = 2:3)
    sparsity <- switch(k$sparsity, none = 0, fixed = 0.05, path = TRUE)
    family <- if (k$model == "var") "Var" else "Mean"
    if (k$model == "mean" && k$zi && k$sparsity != "none") {
      expect_error(normal_block(data, blocks, sparsity = sparsity, zero_inflation = k$zi,
                                model = k$model, control = fast), "sparsity")
      next
    }
    expected <- if (k$sparsity == "path" && k$blocks == "q_list") paste0("NormalBlock", family, "CollectionClustersSparsity")
                else if (k$sparsity == "path") paste0("NormalBlock", family, "CollectionSparsity")
                else if (k$blocks == "q_list") paste0("NormalBlock", family, "CollectionClusters")
                else paste0(if (k$zi) "ZI", "NormalBlock", family, if (k$blocks == "C") "Known" else "Unknown", "Clusters")
    leaf <- paste0(if (k$zi) "ZI", "NormalBlock", family, if (k$blocks == "C") "Known" else "Unknown", "Clusters")

    unfitted <- get_model(data, blocks, sparsity = sparsity, zero_inflation = k$zi,
                          model = k$model, control = fast)
    expect_s3_class(unfitted, expected)
    set.seed(1)
    fit <- normal_block(data, blocks, sparsity = sparsity, zero_inflation = k$zi,
                        model = k$model, control = fast)
    expect_s3_class(fit, expected)
    model <- fit
    while (inherits(model, "NormalBlockCollection")) model <- model$models[[1]]
    expect_s3_class(model, leaf)
    expect_true(isNB(model))
    expect_equal(isNB(fit), !inherits(fit, "NormalBlockCollection"))
    if (k$sparsity == "fixed") expect_equal(model$sparsity, 0.05)
    expect_true(is.finite(model$loglik))
  }
})

test_that("get_model() builds the model without fitting it", {
  model <- get_model(data, 2, control = fast)
  expect_true(anyNA(model$model_par$B))
  expect_false(isNB(data))
})

test_that("best_of_inits applies to unknown clusterings only, and not to sparsity paths", {
  boi <- NB_control(verbose = FALSE, clustering_init = "best_of_inits")
  ## an unknown clustering: the best of the heuristics after short trial runs
  set.seed(1)
  best <- normal_block(data, blocks = 2, control = boi)
  expect_s3_class(best, "NormalBlockVarUnknownClusters")
  loglik_by_init <- sapply(c("kmeans", "ward2", "spectral"), function(init) {
    set.seed(1)
    normal_block(data, 2, control = NB_control(verbose = FALSE, clustering_init = init))$loglik
  })
  expect_gte(best$loglik, min(loglik_by_init) - 1e-6)
  ## a known clustering is simply optimized
  set.seed(1); known <- normal_block(data, blocks = C, control = boi)
  set.seed(1); plain <- normal_block(data, blocks = C, control = NB_control(verbose = FALSE))
  expect_s3_class(known, "NormalBlockVarKnownClusters")
  expect_equal(known$loglik, plain$loglik)
  expect_error(plain$best_of_inits(), "only applies to an unknown clustering")
  expect_error(normal_block(data, blocks = 2, sparsity = TRUE, control = boi), "not supported together with sparsity")
  heuristic <- get_model(data, 2, control = NB_control(verbose = FALSE, heuristic = TRUE))
  expect_error(heuristic$best_of_inits(), "requires \\(V\\)EM inference")
})
