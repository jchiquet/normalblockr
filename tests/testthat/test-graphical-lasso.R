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
  ## so starting closer exits sooner and lands slightly short. That trade is
  ## why the (V)EM does *not* warm-start between M-steps -- see the note in
  ## src/graphical_lasso.h.
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
