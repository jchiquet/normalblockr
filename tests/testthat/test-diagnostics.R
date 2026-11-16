###############################################################################
## Diagnostics: within_block_correlations() on variance-block models (blocks
## must hold positively co-varying variables), and deviance_violations / the
## print() warning on collections over q (nested models: the deviance cannot
## increase with q at the global optimum).

## Block 1: three noisy copies of W1. Block 2: two copies of W2 and one of -W2,
## i.e. an anticorrelated variable the model cannot hold. Variable 7 alone.
set.seed(11)
n  <- 80
W  <- matrix(rnorm(2 * n), n, 2)
Y  <- cbind(W[, c(1, 1, 1, 2, 2)], -W[, 2], rnorm(n)) + matrix(rnorm(7 * n, sd = .3), n, 7)
X  <- matrix(1, n, 1)
C  <- model.matrix(~ 0 + factor(c(1, 1, 1, 2, 2, 2, 3)))
dat <- NormalBlockData$new(Y, X)

test_that("within_block_correlations() flags a block holding anticorrelated variables", {
  fit <- normal_block(dat, blocks = C, control = NB_control(verbose = FALSE))
  wb  <- fit$within_block_correlations()
  expect_equal(names(wb), c("block", "size", "mean_cor", "min_cor", "prop_negative", "model_cor"))
  expect_equal(wb$size, c(3L, 3L, 1L))
  ## block 1 is well specified
  expect_equal(wb$prop_negative[1], 0)
  expect_gt(wb$min_cor[1], 0.5)
  ## block 2: two of its three pairs are negatively correlated
  expect_equal(wb$prop_negative[2], 2 / 3)
  expect_lt(wb$min_cor[2], -0.5)
  ## the implied correlations are always in (0, 1)
  expect_true(all(wb$model_cor[1:2] > 0 & wb$model_cor[1:2] < 1))
  ## a singleton block has no pair
  expect_true(all(is.na(wb[3, c("mean_cor", "min_cor", "prop_negative", "model_cor")])))
})

test_that("within_block_correlations() works on zero-inflated models, ignoring the zeros", {
  Yzi <- Y; Yzi[sample(length(Yzi), 50)] <- 0
  fit <- normal_block(NormalBlockData$new(Yzi, X), blocks = C, zero_inflation = TRUE,
                      control = NB_control(verbose = FALSE))
  wb <- fit$within_block_correlations()
  expect_equal(nrow(wb), 3)
  expect_true(all(is.finite(wb$mean_cor[1:2])))
  expect_lt(wb$min_cor[2], 0)
})

test_that("deviance_violations lists the q whose deviance exceeds the previous one, and print() warns", {
  testdata <- readRDS("testdata/testdata_normal.RDS")
  d  <- NormalBlockData$new(testdata$Y, testdata$X)
  set.seed(1)
  co <- normal_block(d, blocks = 1:4, control = NB_control(verbose = FALSE))
  dev <- co$criteria$deviance[order(co$q_list)]
  expect_equal(co$deviance_violations, sort(co$q_list)[-1][diff(dev) > 1e-6 * pmax(1, abs(dev[-1]))])
  printed <- paste(capture.output(print(co)), collapse = "\n")
  if (length(co$deviance_violations) == 0) {
    expect_no_match(printed, "deviance increases")
  } else {
    expect_match(printed, "deviance increases")
  }

  ## replace the q = 3 model by a barely optimized one: its deviance now
  ## exceeds the q = 2 one, which nested models forbid at the optimum
  bad <- NormalBlockVarUnknownClusters$new(d, 3, control = NB_control(verbose = FALSE))
  bad$optimize(control = list(niter = 1, threshold = -1), warn = FALSE)
  i3  <- which(co$q_list == 3)
  co$models[[i3]] <- bad
  skip_if_not(bad$deviance > co$models[[which(co$q_list == 2)]]$deviance)
  expect_true(3 %in% co$deviance_violations)
  expect_match(paste(capture.output(print(co)), collapse = "\n"), "deviance increases with q at q = 3")
})
