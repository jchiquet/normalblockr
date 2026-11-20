###############################################################################
## Model-selection criteria (R/NormalBlockBase.R active bindings): the BIC
## penalty charges each group of parameters the log of the number of
## observations informing it (Daudin, Picard & Robin 2008; Keribin et al.
## 2015), the ICL completes only the discrete partition (Biernacki, Celeux &
## Govaert 2000), and the EBIC uses the dimension of the network. These are
## structural checks: they do not depend on the quality of the fit.
ctrl <- NB_control(verbose = FALSE, niter = 20, threshold = -1)

set.seed(1)
ex_var  <- generate_normal_block_var_data(n = 60, p = 15, d = 2, q = 3)
dat_var <- NormalBlockData$new(ex_var$Y, ex_var$X)
ex_mean  <- generate_normal_block_mean_data(n = 60, p = 15, d = 2, q = 3)
dat_mean <- NormalBlockData$new(ex_mean$Y, ex_mean$X)
C_known <- ex_var$parameters$C
set.seed(3)
fit_unknown <- normal_block(dat_var, blocks = 3, control = ctrl)

check_bic <- function(fit) {
  pen <- fit$penalty_terms
  expect_equal(fit$nb_param, sum(pen$nu))
  expect_equal(fit$penalty, sum(pen$nu * log(pen$N)))
  expect_equal(fit$BIC, fit$deviance + fit$penalty)
}

test_that("variance-block, unknown clusters: alpha is charged log(p), the rest log(n)", {
  fit <- fit_unknown
  check_bic(fit)
  pen <- fit$penalty_terms
  n <- fit$n; p <- fit$p; q <- fit$q; d <- fit$d
  expect_equal(pen$group, c("B", "Omega", "noise", "alpha"))
  expect_equal(pen$nu, c(p * d, q + fit$n_edges, p, q - 1))
  expect_equal(pen$N, c(n, n, n, p))
})

test_that("a spherical noise variance is charged log(n * p)", {
  fit <- normal_block(dat_var, blocks = 3,
                      control = NB_control(verbose = FALSE, niter = 20, threshold = -1, noise_covariance = "spherical"))
  pen <- fit$penalty_terms
  expect_equal(pen[pen$group == "noise", c("nu", "N")],
               data.frame(nu = 1L, N = fit$n * fit$p), ignore_attr = TRUE)
})

test_that("ICL completes the partition only: ICL = BIC + 2 H(tau), and ICL = BIC for known clusters", {
  fit <- fit_unknown
  tau <- fit$var_par$tau
  expect_equal(fit$entropy, -sum(ifelse(tau < .Machine$double.eps, 0, tau * log(tau))))
  expect_equal(fit$ICL, fit$BIC + 2 * fit$entropy)

  known_var <- normal_block(dat_var, blocks = C_known, control = ctrl)
  expect_equal(known_var$ICL, known_var$BIC)
  known_mean <- normal_block(dat_mean, blocks = ex_mean$parameters$C, model = "mean", control = ctrl)
  expect_equal(known_mean$ICL, known_mean$BIC)
})

test_that("mean-block: cluster profiles are charged log(n * p), alpha log(p)", {
  set.seed(4)
  fit <- normal_block(dat_mean, blocks = 3, model = "mean", control = ctrl)
  check_bic(fit)
  pen <- fit$penalty_terms
  n <- fit$n; p <- fit$p; q <- fit$q; d <- fit$d
  expect_equal(pen$group, c("B", "noise", "alpha"))
  expect_equal(pen$nu, c(q * d, p, q - 1))
  expect_equal(pen$N, c(n * p, n, p))
})

test_that("EBIC uses the dimension of the network: q for variance-block, p for mean-block", {
  fit_var <- normal_block(dat_var, blocks = 3, sparsity = 0.05, control = ctrl)
  expect_equal(fit_var$EBIC - fit_var$BIC, 2 * fit_var$n_edges * log(fit_var$q))

  fit_mean <- normal_block(dat_mean, blocks = 3, sparsity = 0.05, model = "mean", control = ctrl)
  expect_gt(fit_mean$n_edges, 0)
  expect_equal(fit_mean$EBIC - fit_mean$BIC, 2 * fit_mean$n_edges * log(fit_mean$p))
})

test_that("zero-inflated models add p * d0 coefficients charged log(n)", {
  set.seed(2)
  ex  <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3, kappa = rep(0.3, 15))
  dat <- NormalBlockData$new(ex$Y, ex$X)
  fit <- normal_block(dat, blocks = 3, zero_inflation = TRUE, control = NB_control(verbose = FALSE))
  check_bic(fit)
  pen <- fit$penalty_terms
  expect_equal(pen[pen$group == "kappa", c("nu", "N")],
               data.frame(nu = fit$p * fit$data$d0, N = fit$n), ignore_attr = TRUE)
})
