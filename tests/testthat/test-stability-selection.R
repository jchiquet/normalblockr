###############################################################################
## StARS (stability_selection(), R/NormalBlockCollectionSparsity.R), shared by
## both families: each penalty path is refitted on subsamples with the
## clustering held fixed (known, or taken from the least penalized model), and
## the selection frequency of every edge is recorded at every penalty.

test_that("StARS works on both families, with known or unknown clusters", {
  set.seed(2)
  ex_var  <- generate_normal_block_var_data(n = 80, p = 12, d = 1, q = 3)
  ex_mean <- generate_normal_block_mean_data(n = 80, p = 12, d = 1, q = 3)
  ctrl <- NB_control(verbose = FALSE, n_sparsity_penalties = 4, threshold = 1e-3)
  for (family in c("var", "mean")) {
    ex <- if (family == "var") ex_var else ex_mean
    d  <- NormalBlockData$new(ex$Y, ex$X)
    n_nodes <- if (family == "var") 3 else 12
    for (blocks in list(ex$parameters$C, 3)) {
      set.seed(1)
      path <- normal_block(d, blocks = blocks, sparsity = TRUE, model = family, control = ctrl)
      expect_true(all(is.na(path$stability)))
      set.seed(3)
      path$stability_selection(n_subsamples = 3)
      sp <- path$stability_path
      expect_equal(nrow(sp), choose(n_nodes, 2) * 4)
      expect_true(all(sp$Prob >= 0 & sp$Prob <= 1))
      expect_true(all(sp$Prob * 3 == round(sp$Prob * 3)))   # a frequency over 3 subsamples
      expect_equal(sort(unique(sp$Penalty)), sort(path$sparsity))
      expect_match(sp$Edge[1], if (family == "var") "^Cluster_1\\|Cluster_2$" else "^Y1\\|Y2$")
      expect_true(all(path$stability >= 0 & path$stability <= 1))
      expect_equal(path$criteria$stability, path$stability)
      best <- path$get_best_model("StARS", stability = 0)
      expect_equal(best$sparsity, min(path$sparsity))
      expect_s3_class(best, if (family == "var") "NormalBlockVarBase" else "NormalBlockMeanBase")
      expect_message(path$get_best_model("StARS", stability = 1.1), "highest stability")
    }
  }
})

test_that("the subsample paths reuse the collection's penalties, so the stability path is labelled right", {
  set.seed(4)
  ex <- generate_normal_block_mean_data(n = 60, p = 10, d = 1, q = 2)
  path <- normal_block(NormalBlockData$new(ex$Y, ex$X), blocks = ex$parameters$C, sparsity = TRUE,
                       model = "mean", control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
  ## a subsample would otherwise derive its own penalty grid from its own
  ## covariance, and the probabilities would be filed under the wrong penalty
  ctrl <- path$control; ctrl$sparsity_penalties <- path$sparsity
  sub_path <- path$.__enclos_env__$private$subsample_path(path$data$subset(1:40), ctrl)
  expect_equal(sub_path$sparsity, path$sparsity)
})

test_that("the subsamples of a zero-inflated path keep X0 and the zeros mask", {
  set.seed(5)
  ex <- generate_normal_block_var_data(n = 60, p = 10, d = 1, q = 2, kappa = rep(0.2, 10))
  X0 <- cbind(1, rnorm(60))
  zeros <- 1 * (ex$Y == 0); zeros[1:3, 1] <- 1          # one more structural zero
  d <- NormalBlockData$new(ex$Y, ex$X, X0 = X0, zeros = zeros)
  set.seed(1)
  path <- normal_block(d, blocks = 2, sparsity = TRUE, zero_inflation = TRUE,
                       control = NB_control(verbose = FALSE, n_sparsity_penalties = 3))
  seen <- list()
  priv <- path$.__enclos_env__$private
  original <- priv$subsample_path
  unlockBinding("subsample_path", priv)
  priv$subsample_path <- function(data, control) {
    seen[[length(seen) + 1]] <<- data
    original(data, control)
  }
  subsamples <- list(1:45, 11:55)
  path$stability_selection(subsamples = subsamples)
  expect_length(seen, 2)
  for (i in 1:2) {
    expect_equal(seen[[i]]$X0, X0[subsamples[[i]], ])
    expect_equal(seen[[i]]$zeros, zeros[subsamples[[i]], ])
  }
  expect_false(anyNA(path$stability))
})
