###############################################################################
## NB_control(heuristic = TRUE): the moment-based fit, without (V)EM, hence
## without likelihood, for every model of both families.
check_heuristic_fit <- function(model, Y) {
  expect_equal(model$inference_method, "heuristic")
  expect_true(is.na(model$loglik))
  expect_true(is.na(model$BIC))
  expect_equal(dim(model$fitted), dim(Y))
  expect_true(all(is.finite(model$fitted)))
  expect_length(model$clustering, ncol(Y))
}
heuristic <- NB_control(verbose = FALSE, heuristic = TRUE)

test_that("variance-block family: known or unknown clusters, zero-inflated or not", {
  for (zi in c(FALSE, TRUE)) {
    td <- readRDS(if (zi) "testdata/testdata_normal_zi.RDS" else "testdata/testdata_normal.RDS")
    data <- NormalBlockData$new(td$Y, td$X)
    C <- td$parameters$C
    for (blocks in list(C, ncol(C))) {
      set.seed(1)
      model <- normal_block(data, blocks, sparsity = 0.05, zero_inflation = zi, control = heuristic)
      check_heuristic_fit(model, td$Y)
      if (!is.matrix(blocks)) expect_equal(aricode::ARI(model$clustering, get_clusters(C)), 1)
    }
  }
  ## a sparsity path, with its StARS selection
  set.seed(1)
  path <- normal_block(data, C, sparsity = TRUE,
                       control = NB_control(verbose = FALSE, heuristic = TRUE, n_sparsity_penalties = 3,
                                            noise_covariance = "spherical"))
  path$stability_selection(n_subsamples = 2)
  check_heuristic_fit(path$get_best_model("StARS", stability = 0), td$Y)
})

test_that("mean-block family: known or unknown clusters, zero-inflated or not", {
  td <- readRDS("testdata/testdata_normal_mean_block_zi.RDS")
  data <- NormalBlockData$new(td$Yzi, td$X)
  truth <- get_clusters(td$C)
  for (zi in c(FALSE, TRUE)) for (blocks in list(td$C, 3)) {
    set.seed(1)
    model <- normal_block(data, blocks, model = "mean", zero_inflation = zi, control = heuristic)
    check_heuristic_fit(model, td$Yzi)
    expect_equal(model$criteria$niter, 0)
  }
})

test_that("a heuristic q-collection is fitted without refine(), which refuses it", {
  td <- readRDS("testdata/testdata_normal.RDS")
  set.seed(1)
  coll <- normal_block(NormalBlockData$new(td$Y, td$X), 2:4, control = heuristic)
  for (model in coll$models) check_heuristic_fit(model, td$Y)
  expect_error(coll$refine(), "heuristic")
  expect_error(coll$get_best_model("BIC"), "heuristic")
})
