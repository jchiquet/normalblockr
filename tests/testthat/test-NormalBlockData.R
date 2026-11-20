###############################################################################
## Regression tests for NormalBlockData's scale argument (R/NormalBlockData.R):
## with scale = TRUE each column of Y is divided by its own standard deviation
## (no centering -- the model's own intercept in X is expected to absorb the
## mean). The default, scale = NULL, defers to the model family: rescaled for
## variance-block models, unscaled for mean-block ones.
set.seed(1)
n <- 40; p <- 6
Y <- matrix(rnorm(n * p), n, p) * matrix(c(1, 10, 100, 0.1, 1000, 1), n, p, byrow = TRUE)
X <- matrix(1, n, 1)

test_that("scale = TRUE rescales each column to unit variance, without centering", {
  data <- NormalBlockData$new(Y, X, scale = TRUE)
  expect_equal(apply(data$Y, 2, sd), rep(1, p), tolerance = 1e-10)
  expect_equal(data$Y_scale, apply(Y, 2, sd))
  ## not centered: dividing by sd alone does not zero out the column means
  expect_false(isTRUE(all.equal(colMeans(data$Y), rep(0, p))))
  expect_equal(data$Y, Y / matrix(data$Y_scale, n, p, byrow = TRUE))
})

test_that("scale = FALSE leaves Y untouched", {
  data <- NormalBlockData$new(Y, X, scale = FALSE)
  expect_equal(data$Y, Y)
  expect_equal(data$Y_scale, rep(1, p))
})

test_that("scaling does not change which entries are zero (zeros/zeros_bar unaffected)", {
  Yz <- Y; Yz[1:5, 1] <- 0
  data_scaled   <- NormalBlockData$new(Yz, X, scale = TRUE)
  data_unscaled <- NormalBlockData$new(Yz, X, scale = FALSE)
  expect_equal(data_scaled$zeros_bar, data_unscaled$zeros_bar)
  expect_equal(data_scaled$nY, data_unscaled$nY)
})

test_that("constant columns of Y are refused, and named", {
  Yconst <- Y; Yconst[, c(1, 3)] <- 0
  expect_error(NormalBlockData$new(Yconst, X), "constant columns.*: 1, 3")
  colnames(Yconst) <- paste0("v", seq_len(ncol(Y)))
  expect_error(NormalBlockData$new(Yconst, X, scale = TRUE), "v1, v3")
  ## allowed for data only used to evaluate a fitted model
  expect_silent(NormalBlockData$new(Yconst, X, check_constant = FALSE))
})

test_that("scale = NULL (default) defers to the family: rescaled for var, unscaled for mean", {
  data <- NormalBlockData$new(Y, X)
  expect_null(data$scale)
  expect_equal(data$Y, Y)
  var_data <- data$for_family("var")
  expect_equal(var_data$Y, NormalBlockData$new(Y, X, scale = TRUE)$Y)
  expect_identical(data$for_family("var"), var_data)   # built once, shared
  expect_identical(data$for_family("mean"), data)
  ## by every model of a collection, hence one memoized OLS fit
  coll <- get_model(data, 2:3, control = NB_control(verbose = FALSE))
  for (m in coll$models) expect_identical(m$data, var_data)
  ## an explicit choice is kept by both families
  for (sc in c(TRUE, FALSE)) {
    d <- NormalBlockData$new(Y, X, scale = sc)
    expect_identical(d$for_family("var"), d)
    expect_identical(d$for_family("mean"), d)
  }
  expect_error(NormalBlockData$new(Y, X, scale = NA), "scale must be")
})

test_that("models fitted on default data use their family's scale, whatever the entry point", {
  C <- model.matrix(~ 0 + factor(rep(1:2, each = 3)))
  data <- NormalBlockData$new(Y, X)
  ctrl <- NB_control(verbose = FALSE)
  expect_equal(normal_block(data, blocks = C, control = ctrl)$data$Y_scale, apply(Y, 2, sd))
  expect_equal(normal_block(data, blocks = C, model = "mean", control = ctrl)$data$Y_scale, rep(1, p))
  expect_equal(NormalBlockVarKnownClusters$new(data, C)$data$Y_scale, apply(Y, 2, sd))
  expect_equal(NormalBlockMeanKnownClusters$new(data, C)$data$Y_scale, rep(1, p))
  ## a collection shares one data object across its models
  co <- normal_block(data, blocks = 1:2, control = ctrl)
  expect_identical(co$models[[1]]$data, co$models[[2]]$data)
})

test_that("subset() keeps X0, the zeros mask, the formula and the scale request; Y_original is Y as given", {
  Xf <- cbind(X1 = rnorm(n), X2 = rnorm(n))
  X0 <- cbind(Z1 = 1, Z2 = rnorm(n))
  zeros <- matrix(rbinom(n * p, 1, 0.2), n, p)
  d <- NormalBlockData$new(Y, Xf, X0 = X0, formula = ~ X1 | Z2, scale = TRUE, zeros = zeros)
  expect_equal(d$Y_original, Y)
  rows <- 3:30
  s <- d$subset(rows)
  expect_equal(s$Y_original, Y[rows, ])
  expect_true(isTRUE(s$scale))
  expect_equal(s$Y, NormalBlockData$new(Y[rows, ], Xf[rows, ], scale = TRUE)$Y)
  expect_equal(s$X, d$X[rows, ], ignore_attr = TRUE)
  expect_equal(s$X0, d$X0[rows, ], ignore_attr = TRUE)
  expect_equal(s$zeros, zeros[rows, ])
  expect_identical(s$formula, d$formula)
  expect_null(NormalBlockData$new(Y, X)$subset(rows)$scale)
})

test_that("an explicit zeros mask overrides the default Y == 0", {
  mask <- matrix(0, n, p); mask[1:5, 1:2] <- 1
  default  <- NormalBlockData$new(Y, X)
  explicit <- NormalBlockData$new(Y, X, zeros = mask)
  expect_equal(sum(default$zeros), 0)
  expect_equal(explicit$zeros, mask)
  expect_equal(explicit$zeros_bar, 1 - mask)
  expect_equal(explicit$npY, n * p - 10)
  expect_equal(explicit$nY, colSums(1 - mask))
  expect_error(NormalBlockData$new(Y, X, zeros = mask[, -1]), "same dimensions")
})
