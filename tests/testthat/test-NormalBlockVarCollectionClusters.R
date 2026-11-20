###############################################################################
## Collections over a range of q (R/NormalBlockCollectionClusters.R): a per-q
## clustering_init list, the defaults, and refine()'s split/merge search.
## Frozen fixture (seeds are unreliable inside testthat)
testdata <- readRDS("testdata/testdata_normal_zi.RDS")
Y <- testdata$Y ; X <- testdata$X
C <- testdata$parameters$C ; q <- ncol(C)

test_that("a collection starts each q from its own clustering_init", {
  clustering_init_3 <- get_clusters(C)
  clustering_init_2 <- clustering_init_3 ; clustering_init_2[clustering_init_2 == 3] <- 2
  clustering_init_1 <- clustering_init_2 ; clustering_init_1[clustering_init_1 == 2] <- 1
  clustering_init   <- list(clustering_init_1, clustering_init_2, clustering_init_3)

  data  <- NormalBlockData$new(Y, X)
  model <- NormalBlockVarCollectionClusters$new(data, 1:3,
             control = NB_control(clustering_init = clustering_init, verbose = FALSE))
  model$optimize(NB_control(verbose = FALSE))
  expect_lt(model$get_best_model("BIC")$BIC, 7455)
  expect_equal(model$get_best_model("BIC")$q, 3)
})

set.seed(7)
ex2 <- generate_normal_block_var_data(n = 60, p = 30, d = 1, q = 5)
data2 <- NormalBlockData$new(ex2$Y, ex2$X)
## a collection over 2:5, optimized without refine()
set.seed(1)
coll0 <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE, refine = FALSE))
coll0$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
deviances <- function(coll) sapply(coll$models, function(m) m$deviance)
## a copy of coll0 whose models can be modified
fresh <- function() {
  coll <- coll0$clone()
  coll$models <- lapply(coll0$models, function(m) m$clone())
  coll
}
## replaces a public method of every model of a collection
stub <- function(coll, method, fun) {
  for (m in coll$models) {
    unlockBinding(method, m)
    assign(method, fun, envir = m)
  }
}

test_that("NB_control()'s refine option defaults to TRUE, and kmeans is both families' initialization", {
  expect_true(NB_control()$refine)
  ## a single model resolves the heuristic itself...
  for (cls in list(NormalBlockVarUnknownClusters, NormalBlockMeanUnknownClusters))
    expect_equal(cls$new(data2, 3)$.__enclos_env__$private$clustering_approx, "kmeans")
  ## ... a collection through one shared path, the kmeans one by default
  for (family in c("var", "mean")) {
    set.seed(1); default <- normalblockr:::clustering_path_for_family(data2, 2:4, family)
    set.seed(1); km <- normalblockr:::clustering_path_for_family(data2, 2:4, family,
                                                                 control = NB_control(clustering_init = "kmeans"))
    expect_equal(default, km)
  }
  ## refine = TRUE leaves no deviance increase on these well-separated data
  set.seed(1)
  co <- normal_block(data2, blocks = 2:4, control = NB_control(verbose = FALSE))
  expect_length(co$deviance_violations, 0)
})

test_that("refine() never makes a model's deviance worse, and skips non-contiguous q gaps", {
  coll <- fresh()
  coll$models <- coll$models[-3]                  # q = 2, 3, 5
  stub(coll, "candidates_merge", function(...) list())
  ## q = 5 has no q = 4 neighbour to be split from
  tried <- c()
  for (m in coll$models[1:2]) {
    unlockBinding("candidates_split", m)
    local({ q <- m$q; f <- m$candidates_split
            assign("candidates_split", function(...) { tried <<- c(tried, q); f(...) }, envir = m) })
  }
  set.seed(1)
  expect_no_error(coll$refine(verbose = FALSE))
  expect_equal(tried, 2)                          # only q = 2 seeds a candidate (for q = 3)
  expect_true(all(deviances(coll) <= deviances(coll0)[-3] + 1e-8))
})

test_that("refine(directions = ) only builds the candidates it is asked for", {
  never <- function(...) stop("not this direction")
  coll <- fresh(); stub(coll, "candidates_merge", never)
  set.seed(1); expect_no_error(coll$refine(directions = "split", verbose = FALSE))
  coll <- fresh(); stub(coll, "candidates_split", never)
  set.seed(1); expect_no_error(coll$refine(directions = "merge", verbose = FALSE))
  expect_error(fresh()$refine(directions = "swap"))
})

test_that("refine(split_rule = ) reaches the split of every model", {
  coll <- fresh()
  rules <- c()
  stub(coll, "candidates_split", function(trial_niter = 5, split_rule = NULL) {
    rules <<- c(rules, split_rule); list()
  })
  coll$refine(directions = "split", split_rule = "noise", verbose = FALSE)
  expect_equal(rules, rep("noise", 3))
  expect_error(coll$refine(split_rule = "random"), "should be one of")
})

test_that("refine() with both directions is never worse than either alone", {
  coll_split <- fresh(); coll_merge <- fresh(); coll_both <- fresh()
  set.seed(1); coll_split$refine(directions = "split", verbose = FALSE)
  set.seed(1); coll_merge$refine(directions = "merge", verbose = FALSE)
  set.seed(1); coll_both$refine(directions = c("split", "merge"), verbose = FALSE)
  expect_true(all(deviances(coll_both) <= deviances(coll_split) + 1e-8))
  expect_true(all(deviances(coll_both) <= deviances(coll_merge) + 1e-8))
  expect_true(all(deviances(coll_both) <= deviances(coll0) + 1e-8))
})

test_that("NB_control(refine = TRUE) makes optimize() call refine() automatically", {
  set.seed(1)
  coll_auto <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE, refine = TRUE))
  coll_auto$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  set.seed(1)
  coll_manual <- NormalBlockVarCollectionClusters$new(data2, 2:5, control = NB_control(verbose = FALSE, refine = FALSE))
  coll_manual$optimize(control = list(niter = 100, threshold = 1e-4, verbose = FALSE))
  coll_manual$refine(verbose = FALSE)
  expect_equal(deviances(coll_auto), deviances(coll_manual))
})

test_that("a collection is sorted by q, with a per-q clustering_init list kept aligned", {
  data <- NormalBlockData$new(Y, X)
  cl3 <- normalblockr:::get_clusters(C)
  cl2 <- cl3; cl2[cl2 == 3] <- 2
  ctrl <- NB_control(verbose = FALSE, refine = FALSE, clustering_init = list(cl3, cl2))
  for (family in c("var", "mean")) {
    coll <- get_model(data, c(3, 2), model = family, control = ctrl)
    expect_equal(coll$q_list, c(2, 3))
    expect_equal(unname(coll$get_model(2)$clustering), cl2)
    expect_equal(unname(coll$get_model(3)$clustering), cl3)
  }
  path <- get_model(data, c(3, 2), sparsity = TRUE,
                    control = NB_control(verbose = FALSE, n_sparsity_penalties = 2))
  expect_equal(path$q_list, c(2, 3))
})

test_that("refine() discards a candidate whose trial fit fails instead of aborting", {
  data  <- NormalBlockData$new(Y, X)
  model <- normal_block(data, 2, control = NB_control(verbose = FALSE))
  good  <- model$split(1)
  bad   <- model$split(1)
  penv  <- bad$.__enclos_env__$private
  unlockBinding("optimizer", penv)
  penv$optimizer <- function(control) stop("degenerate")
  kept  <- model$.__enclos_env__$private$trial_fits(list(good, bad, good$clone()), 2)
  expect_length(kept, 2)
})
