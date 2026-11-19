set.seed(11)
ex_var  <- generate_normal_block_var_data(n = 60, p = 12, d = 2, q = 3)
ex_mean <- generate_normal_block_mean_data(n = 60, p = 12, d = 2, q = 3)
name_data <- function(ex) {
  Y <- ex$Y; X <- ex$X
  colnames(Y) <- paste0("var", seq_len(ncol(Y))); rownames(Y) <- paste0("obs", seq_len(nrow(Y)))
  colnames(X) <- paste0("cov", seq_len(ncol(X)))
  NormalBlockData$new(Y, X)
}

test_that("the extractors keep the names of the variables, covariates and clusters", {
  data <- name_data(ex_var)
  m <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
  expect_equal(dimnames(coef(m)), list(paste0("cov", 1:2), paste0("var", 1:12)))
  expect_equal(dimnames(sigma(m)), list(paste0("Cluster_", 1:3), paste0("Cluster_", 1:3)))
  expect_equal(dimnames(fitted(m)), list(paste0("obs", 1:60), paste0("var", 1:12)))
  new_X <- data$X[1:2, , drop = FALSE]
  rownames(new_X) <- c("a", "b")
  expect_equal(dimnames(predict(m, new_X)), list(c("a", "b"), paste0("var", 1:12)))

  data <- name_data(ex_mean)
  m <- normal_block(data, blocks = 3, model = "mean", control = NB_control(verbose = FALSE))
  expect_equal(dimnames(coef(m)), list(paste0("cov", 1:2), paste0("Cluster_", 1:3)))
  expect_equal(dimnames(sigma(m)), list(paste0("var", 1:12), paste0("var", 1:12)))
  expect_equal(dim(predict(m, data$X[1:2, , drop = FALSE])), c(2, 12))
})

test_that("the extractors work on unnamed data", {
  data <- NormalBlockData$new(ex_var$Y, ex_var$X)
  m <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
  expect_equal(dim(coef(m)), c(2, 12))
  expect_null(colnames(coef(m)))
  expect_equal(colnames(sigma(m)), paste0("Cluster_", 1:3))
})

test_that("print() and summary() end the model type with a period, without a space", {
  data <- NormalBlockData$new(ex_var$Y, ex_var$X)
  m <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
  expect_false(any(grepl(" \\.$", capture.output(print(m)))))
  expect_false(any(grepl(" \\.$", capture.output(print(summary(m))))))
})

test_that("coef() of a variance-block model is on the original scale of Y", {
  Y <- ex_var$Y * matrix(c(1, 10), nrow(ex_var$Y), ncol(ex_var$Y), byrow = TRUE)
  data <- NormalBlockData$new(Y, ex_var$X, scale = TRUE)
  m <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
  expect_equal(predict(m, ex_var$X), ex_var$X %*% coef(m), ignore_attr = TRUE)
})
