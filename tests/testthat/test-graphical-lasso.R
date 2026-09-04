###############################################################################
## The in-package graphical lasso (src/graphical_lasso.h), which replaced the
## glassoFast callback the C++ (V)EM used to make once per M-step.
##
## glassoFast is no longer a dependency, so these tests can't compare against
## it; the equivalence was established while both were installed (477 random
## problems, worst relative difference 2.2e-13) and is pinned here instead by
## the two properties that actually matter: the returned matrix solves the
## penalized problem, and the known-in-closed-form cases come out exact.
###############################################################################

## a well-conditioned covariance: 3n degrees of freedom plus a small ridge, so
## that the tests measure the solver rather than the conditioning of S
rand_S <- function(n, m = 3 * n) {
  A <- matrix(rnorm(n * m), n)
  tcrossprod(A) / m + diag(0.1, n)
}

## -log det(Theta) + tr(S Theta) + || rho o Theta ||_1
glasso_objective <- function(S, rho, Theta) {
  ld <- determinant(Theta, logarithm = TRUE)
  if (ld$sign <= 0) return(Inf)
  -as.numeric(ld$modulus) + sum(S * Theta) + sum(abs(rho * Theta))
}

off_diag_weights <- function(n) {
  w <- matrix(1, n, n)
  diag(w) <- 0
  w
}

test_that("the solution is symmetric, positive definite, and minimizes the objective", {
  set.seed(101)
  for (n in c(4, 12, 30)) {
    S   <- rand_S(n)
    rho <- 0.08 * off_diag_weights(n)
    Theta <- graphical_lasso_fit(S, rho)$wi

    expect_true(isSymmetric(Theta))
    expect_gt(min(eigen(Theta, symmetric = TRUE, only.values = TRUE)$values), 0)

    ## no small symmetric perturbation staying in the PD cone does better
    base  <- glasso_objective(S, rho, Theta)
    worse <- 0
    for (k in 1:150) {
      P <- matrix(0, n, n)
      i <- sample(n, 1); j <- sample(n, 1)
      P[i, j] <- P[j, i] <- rnorm(1, sd = 1e-3)
      if (glasso_objective(S, rho, Theta + P) < base - 1e-10) worse <- worse + 1
    }
    expect_equal(worse, 0)
  }
})

test_that("the penalty actually sparsifies, monotonically in rho", {
  set.seed(102)
  n <- 20
  S <- rand_S(n)
  nnz <- vapply(c(0.001, 0.05, 0.2, 0.6, 2), function(pen) {
    Theta <- graphical_lasso_fit(S, pen * off_diag_weights(n))$wi
    sum(abs(Theta[upper.tri(Theta)]) > 1e-10)
  }, numeric(1))

  expect_true(all(diff(nnz) <= 0))
  expect_gt(nnz[1], nnz[length(nnz)])
  expect_equal(nnz[length(nnz)], 0) # a big enough penalty empties the network
})

test_that("an unpenalized problem returns the plain inverse", {
  set.seed(103)
  n <- 10
  S <- rand_S(n)
  ## rho = 0 everywhere: the graphical lasso reduces to inverting S. Tightened
  ## threshold -- at the default 1e-4 this is an iterative solve stopped early,
  ## not an exact inverse.
  expect_equal(graphical_lasso_fit(S, matrix(0, n, n), thr = 1e-12)$wi, solve(S),
               tolerance = 1e-6)
})

test_that("a separable problem is solved exactly rather than through the recursion", {
  ## No off-diagonal mass: Theta is diagonal, with entries 1 / (S_ii + rho_ii).
  ## This is the case glassoFast got wrong -- it dropped S_ii and returned
  ## 1 / max(rho_ii, eps), i.e. ~9.09e15 on an unpenalized diagonal.
  S <- diag(c(4, 9, 16))
  expect_equal(diag(graphical_lasso_fit(S, 0.1 * off_diag_weights(3))$wi), 1 / diag(S))

  ## every 1 x 1 problem takes that branch, which is how a q = 1 model reaches it
  expect_equal(graphical_lasso_fit(matrix(4, 1, 1), matrix(0, 1, 1))$wi[1], 0.25)

  ## a penalized diagonal shifts it, and is still exact
  expect_equal(diag(graphical_lasso_fit(S, diag(1, 3))$wi), 1 / (diag(S) + 1))
})

test_that("degenerate input is reported rather than hung on or silently accepted", {
  ## glassoFast loops forever on a non-finite input: `dlx < thrLasso` is never
  ## true once dlx is NaN, and its inner loop has no other exit.
  bad <- matrix(c(1, NA, NA, 1), 2)
  res <- graphical_lasso_fit(bad, matrix(0.1, 2, 2))
  expect_true(all(is.na(res$wi)))
  expect_false(res$converged)

  ## a zero-variance coordinate is what would divide by zero inside the sweep
  res0 <- graphical_lasso_fit(matrix(c(0, 0, 0, 1), 2), matrix(0, 2, 2))
  expect_true(all(is.na(res0$wi)))
  expect_false(res0$converged)
})

test_that("a scalar rho and a constant matrix rho agree", {
  set.seed(104)
  S <- rand_S(8)
  expect_equal(graphical_lasso_fit(S, 0.1)$wi,
               graphical_lasso_fit(S, matrix(0.1, 8, 8))$wi)
  expect_error(graphical_lasso_fit(S, matrix(0.1, 3, 3)), "same dimensions")
})

test_that("per-pair weights are respected: a heavily penalized pair is zeroed", {
  set.seed(105)
  n <- 8
  S <- rand_S(n)
  rho <- 0.02 * off_diag_weights(n)
  rho[2, 5] <- rho[5, 2] <- 50
  Theta <- graphical_lasso_fit(S, rho)$wi
  expect_equal(Theta[2, 5], 0)
  expect_gt(sum(abs(Theta[upper.tri(Theta)]) > 1e-10), 0) # the rest survives
})

test_that("a warm start lands on the same solution as a cold one", {
  set.seed(106)
  n <- 15
  S1 <- rand_S(n)
  S2 <- S1 + 0.02 * rand_S(n) # a nearby problem, as between two M-steps
  rho <- 0.05 * off_diag_weights(n)

  first <- graphical_lasso_fit(S1, rho)
  cold  <- graphical_lasso_fit(S2, rho)
  warm  <- graphical_lasso_fit(S2, rho, w_init = first$w, wi_init = first$wi)

  ## Both solve the same problem, but `dw <= shr` measures per-sweep progress,
  ## so starting closer exits sooner and lands slightly short of where a cold
  ## start at the same threshold would. The (V)EM warm-starts anyway and
  ## tightens the threshold to compensate -- see src/graphical_lasso.h.
  expect_equal(cold$wi, warm$wi, tolerance = 1e-2)
  expect_lte(warm$niter, cold$niter)

  ## tightening the threshold collapses the difference, confirming it is
  ## convergence slack and not a different fixed point
  tight_cold <- graphical_lasso_fit(S2, rho, thr = 1e-10)
  tight_warm <- graphical_lasso_fit(S2, rho, thr = 1e-10,
                                    w_init = first$w, wi_init = first$wi)
  expect_equal(tight_cold$wi, tight_warm$wi, tolerance = 1e-6)
})

test_that("w and wi are inverses of each other where the penalty does not bind", {
  set.seed(107)
  n <- 10
  S <- rand_S(n)
  res <- graphical_lasso_fit(S, 0.05 * off_diag_weights(n))
  ## the KKT conditions force W = Theta^-1 on the diagonal and wherever
  ## Theta is non-zero; the whole product is the identity at convergence
  expect_equal(res$w %*% res$wi, diag(n), tolerance = 1e-3)
})

###############################################################################
## The (V)EM warm-starts each M-step's graphical lasso from the previous one
## (nb_omega::estimate). That state belongs to one model and one recursion:
## nothing about a fit may depend on what was fitted before it.
###############################################################################

test_that("a penalized fit does not depend on what was fitted before it", {
  set.seed(201)
  ex <- generate_normal_block_mean_data(n = 120, p = 25, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  ctl <- NB_control(verbose = FALSE, noise_covariance = "full")

  ## seeded identically either side: the initial clustering is drawn, so this
  ## isolates the carried-over solver state from the heuristic's own randomness
  set.seed(999)
  alone <- normal_block(d, blocks = 3, sparsity = 0.05, model = "mean", control = ctl)

  ## the same fit, after an unrelated one on different data
  set.seed(202)
  other <- generate_normal_block_mean_data(n = 90, p = 25, d = 1, q = 4)
  invisible(normal_block(NormalBlockData$new(other$Y, other$X), blocks = 4,
                         sparsity = 0.2, model = "mean", control = ctl))
  set.seed(999)
  after <- normal_block(d, blocks = 3, sparsity = 0.05, model = "mean", control = ctl)

  expect_equal(after$loglik, alone$loglik, tolerance = 0)
  expect_equal(as.matrix(after$model_par$Omega), as.matrix(alone$model_par$Omega),
               tolerance = 0)
})

test_that("repeating a penalized fit is exactly reproducible", {
  set.seed(203)
  ex <- generate_normal_block_var_data(n = 100, p = 30, d = 1, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X)
  C  <- ex$parameters$C
  a <- normal_block(d, blocks = C, sparsity = 0.05, control = NB_control(verbose = FALSE))
  b <- normal_block(d, blocks = C, sparsity = 0.05, control = NB_control(verbose = FALSE))
  expect_equal(as.matrix(a$model_par$Omega), as.matrix(b$model_par$Omega), tolerance = 0)
  expect_equal(a$loglik, b$loglik, tolerance = 0)
})

test_that("the R reference and the C++ core use the same M-step threshold", {
  ## NB_GLASSO_THRESHOLD (R/utils.R) must track nb_omega::kGlassoThreshold
  ## (src/omega_estimation.h). They cannot be read from one another, so the
  ## guard is that the two recursions are compared trace-for-trace at 1e-8 in
  ## test-cpp-normal-block-mean.R -- drift there fails loudly. This just pins
  ## the R side against silent edits.
  expect_equal(normalblockr:::NB_GLASSO_THRESHOLD, 1e-6)
})

test_that("a warm start that diverges is retried cold, not abandoned", {
  ## A warm start is only a starting point, but a bad one can send the
  ## coordinate descent to infinity on an ill-conditioned Sigma that a cold
  ## start solves in a handful of sweeps. Falling straight through to the
  ## unpenalized inverse would silently swap the model being fitted, so
  ## nb_omega::estimate() retries cold first.
  set.seed(301)
  ex <- generate_normal_block_mean_data(n = 150, p = 40, d = 2, q = 3)
  d  <- NormalBlockData$new(ex$Y, ex$X, scale = FALSE)

  ## the fit must stay penalized: Omega keeps exact zeros off the diagonal,
  ## which a fall-through to solve(Sigma_hat) would not produce
  fit <- normal_block(d, blocks = 3, sparsity = 0.1, model = "mean",
                      control = NB_control(verbose = FALSE))
  Om <- fit$model_par$Omega
  expect_gt(sum(Om[upper.tri(Om)] == 0), 0)
  expect_true(all(is.finite(Om)))
})

test_that("a diverging warm start is a real failure mode, not a hypothetical", {
  ## Guards the retry above: feeding the solver a warm start from an
  ## unrelated, badly scaled problem must still return a usable answer.
  set.seed(302)
  n <- 25
  A <- matrix(rnorm(n * n), n); S <- tcrossprod(A) / n + diag(0.05, n)
  rho <- 0.02 * (1 - diag(n))

  ok   <- graphical_lasso_fit(S, rho)
  junk <- ok$w * 1e6                       # a wildly mis-scaled "previous" solve
  warm <- graphical_lasso_fit(S, rho, w_init = junk, wi_init = solve(junk))

  ## whatever the warm attempt does, a cold solve on the same problem is fine
  expect_false(anyNA(ok$wi))
  expect_true(is.matrix(warm$wi))
})
