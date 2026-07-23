###############################################################################
## Tests for sbm_clustering_path() (R/utils.R) and its use inside
## NormalBlockVarCollectionClusters / NormalBlockVarCollectionClustersSparsity: a single wide SBM
## exploration over the whole q_list range replaces one independent
## exploration per q whenever clustering_init = "sbm" is applied uniformly
## (not as an explicit per-q clustering). See inst/normal_block_models.qmd for the
## benchmark motivating this change (the per-q sbm heuristic was ~15-90x more
## expensive than this shared-path version, for the same or statistically
## equivalent quality).
set.seed(99)
ex   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4)
data <- NormalBlockData$new(ex$Y, ex$X)
R    <- ols_residuals(data)

test_that("sbm_clustering_path() returns one valid membership vector per q, named by q, never NULL", {
  q_list <- c(2, 3, 4, 5)
  path <- sbm_clustering_path(R, q_list)

  expect_equal(names(path), as.character(q_list))
  for (q in q_list) {
    membership <- path[[as.character(q)]]
    expect_false(is.null(membership))
    expect_length(membership, data$p)
    expect_equal(length(unique(membership)), q)
  }
})

test_that("sbm_clustering_path() falls back to a cheap ward2 clustering for q's the SBM exploration doesn't reach", {
  ## A q_list reaching far beyond what a tiny dataset's SBM exploration will
  ## ever cover forces the fallback path for the largest values -- regression
  ## test for the cost blow-up found on a real dataset where SBM's own model
  ## selection stopped early (re-running a dedicated SBM call per missing q
  ## reintroduced the cost this function exists to avoid).
  q_list <- 2:20
  path <- sbm_clustering_path(R, q_list)
  expect_true(all(!sapply(path, is.null)))
  for (q in q_list) expect_equal(length(unique(path[[as.character(q)]])), q)
})

test_that("NormalBlockVarCollectionClusters with clustering_init = 'sbm' eagerly resolves it from the shared path", {
  coll <- NormalBlockVarCollectionClusters$new(data, 2:5, control = NB_control(verbose = FALSE, clustering_init = "sbm"))

  for (model in coll$models) {
    ## every model converges to a valid q-cluster initialization once
    ## actually initialized, whether served by the shared path or (for q's it
    ## didn't cover) by the per-q heuristic fallback
    init <- model$.__enclos_env__$private$EM_initialize()
    expect_equal(length(unique(get_clusters(init$C))), model$q)
  }
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})

test_that("an explicit clustering_init is never overridden by the sbm path", {
  cl_init <- list(rep(1:2, length.out = data$p), rep(1:3, length.out = data$p))
  coll <- NormalBlockVarCollectionClusters$new(data, c(2, 3), control = NB_control(
    verbose = FALSE, clustering_init = cl_init
  ))

  for (i in seq_along(coll$models)) {
    C0 <- coll$models[[i]]$.__enclos_env__$private$C
    expect_equal(get_clusters(C0), cl_init[[i]])
  }
})

test_that("NormalBlockVarCollectionClustersSparsity also uses the shared sbm path", {
  coll <- NormalBlockVarCollectionClustersSparsity$new(data, c(2, 3), control = NB_control(
    verbose = FALSE, clustering_init = "sbm", n_sparsity_penalties = 3
  ))
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})

test_that("sbm_clustering_path()'s ward2 fallback does not error on a (near-)constant column", {
  R_const <- R
  R_const[, 1] <- 0; R_const[1, 1] <- 5 # a single non-zero residual for variable 1
  expect_no_warning(path <- sbm_clustering_path(R_const, 2:20))
  for (q in 2:20) expect_equal(length(unique(path[[as.character(q)]])), q)
})

test_that("zero-inflated collections also use the shared sbm path, clustering on zi_residuals() instead of ols_residuals()", {
  exzi   <- generate_normal_block_var_data(n = 60, p = 24, d = 1, q = 4, kappa = rep(0.3, 24))
  datazi <- NormalBlockData$new(exzi$Y, exzi$X, X0 = matrix(1, nrow(exzi$Y), 1))

  coll <- NormalBlockVarCollectionClusters$new(datazi, 2:3, zero_inflation = TRUE,
                                  control = NB_control(verbose = FALSE, clustering_init = "sbm"))
  ## clustering_init is now eagerly injected from the shared path, just like
  ## for non-ZI collections -- private$C is already a valid q-cluster
  ## indicator matrix before optimize()/EM_initialize() ever runs
  for (model in coll$models) {
    C0 <- model$.__enclos_env__$private$C
    expect_false(anyNA(C0))
    expect_equal(length(unique(get_clusters(C0))), model$q)
  }
  expect_no_error(coll$optimize(control = list(niter = 3, threshold = -1, verbose = FALSE)))
})
