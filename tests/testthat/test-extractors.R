###############################################################################
## The extractors (coef(), sigma(), fitted(), predict(), latent_network()):
## names of the variables, covariates and clusters, and Y's original scale.
set.seed(11)
ex_var  <- generate_normal_block_var_data(n = 60, p = 12, d = 2, q = 3)
ex_mean <- generate_normal_block_mean_data(n = 60, p = 12, d = 2, q = 3)
name_data <- function(ex) {
  Y <- ex$Y; X <- ex$X
  colnames(Y) <- paste0("var", seq_len(ncol(Y))); rownames(Y) <- paste0("obs", seq_len(nrow(Y)))
  colnames(X) <- paste0("cov", seq_len(ncol(X)))
  NormalBlockData$new(Y, X)
}
ctrl <- NB_control(verbose = FALSE)
set.seed(1)
m_var <- normal_block(name_data(ex_var), blocks = 3, control = ctrl)

test_that("the extractors keep the names of the variables, covariates and clusters", {
  expect_equal(dimnames(coef(m_var)), list(paste0("cov", 1:2), paste0("var", 1:12)))
  expect_equal(dimnames(sigma(m_var)), list(paste0("Cluster_", 1:3), paste0("Cluster_", 1:3)))
  expect_equal(dimnames(m_var$latent_network()), dimnames(sigma(m_var)))
  expect_equal(dimnames(fitted(m_var)), list(paste0("obs", 1:60), paste0("var", 1:12)))
  new_X <- m_var$data$X[1:2, , drop = FALSE]
  rownames(new_X) <- c("a", "b")
  expect_equal(dimnames(predict(m_var, new_X)), list(c("a", "b"), paste0("var", 1:12)))

  data <- name_data(ex_mean)
  set.seed(1)
  m <- normal_block(data, blocks = 3, model = "mean", control = ctrl)
  expect_equal(dimnames(coef(m)), list(paste0("cov", 1:2), paste0("Cluster_", 1:3)))
  expect_equal(dimnames(sigma(m)), list(paste0("var", 1:12), paste0("var", 1:12)))
  expect_equal(dimnames(m$latent_network()), dimnames(sigma(m)))
  ## predict() maps the cluster-level predictor back to the variables
  expect_equal(dim(predict(m, data$X[1:2, , drop = FALSE])), c(2, 12))
  expect_equal(predict(m, data$X), fitted(m), ignore_attr = TRUE)
})

test_that("the extractors work on unnamed data", {
  set.seed(1)
  m <- normal_block(NormalBlockData$new(ex_var$Y, ex_var$X), blocks = 3, control = ctrl)
  expect_equal(dim(coef(m)), c(2, 12))
  expect_null(colnames(coef(m)))
  expect_equal(colnames(sigma(m)), paste0("Cluster_", 1:3))
})

test_that("the variance-block extractors undo the rescaling of Y", {
  Y <- ex_var$Y * matrix(c(1, 10), nrow(ex_var$Y), ncol(ex_var$Y), byrow = TRUE)
  data <- NormalBlockData$new(Y, ex_var$X, scale = TRUE)
  set.seed(1)
  m <- normal_block(data, blocks = 3, control = ctrl)
  s <- data$Y_scale
  expect_false(all(s == 1))
  B <- m$model_par$B
  expect_equal(m$B_original, B * matrix(s, nrow(B), ncol(B), byrow = TRUE))
  expect_equal(coef(m), m$B_original, ignore_attr = TRUE)
  expect_equal(m$dm1_original, m$model_par$dm1 / s^2)
  expect_equal(predict(m, ex_var$X), ex_var$X %*% coef(m), ignore_attr = TRUE)
  ## fitted adds the posterior mean of the latent signals, on the same scale
  expect_equal(m$fitted, ex_var$X %*% m$B_original +
                 tcrossprod(m$var_par$M, m$var_par$tau) * matrix(s, nrow(Y), ncol(Y), byrow = TRUE))
  ## the heuristic fit has no latent signal: fitted = X B_original = X coef
  h <- normal_block(data, blocks = ex_var$parameters$C, control = NB_control(verbose = FALSE, heuristic = TRUE))
  expect_equal(h$fitted, data$X %*% h$B_original)
  expect_equal(h$fitted, data$X %*% coef(h), ignore_attr = TRUE)
  ## without rescaling, both scales coincide
  u <- normal_block(NormalBlockData$new(Y, ex_var$X, scale = FALSE), blocks = ex_var$parameters$C,
                    control = NB_control(verbose = FALSE, heuristic = TRUE))
  expect_equal(u$B_original, u$model_par$B)
  expect_equal(u$dm1_original, u$model_par$dm1)
})

test_that("the mean-block extractors undo the rescaling of Y", {
  Y <- ex_mean$Y * matrix(c(1, 10), nrow(ex_mean$Y), ncol(ex_mean$Y), byrow = TRUE)
  data <- NormalBlockData$new(Y, ex_mean$X, scale = TRUE)
  set.seed(1)
  m <- normal_block(data, blocks = 3, model = "mean",
                    control = NB_control(verbose = FALSE, noise_covariance = "full"))
  expect_equal(predict(m, ex_mean$X), fitted(m), ignore_attr = TRUE)
  ## the p x p residual covariance is mapped back to Y's units
  s <- data$Y_scale
  expect_equal(sigma(m), chol2inv(chol(m$model_par$Omega)) * tcrossprod(s), ignore_attr = TRUE)
})

test_that("warm_start_from() copies the parameters on the scale of the fit", {
  data <- NormalBlockData$new(ex_var$Y, ex_var$X)
  model <- NormalBlockVarUnknownClusters$new(data, 3, control = ctrl)
  model$warm_start_from(m_var)
  expect_equal(model$.__enclos_env__$private$B, m_var$model_par$B)
  expect_true(model$.__enclos_env__$private$warm_started)
  expect_error(model$warm_start_from(normal_block(data, 2, control = ctrl)), "same q")
})
