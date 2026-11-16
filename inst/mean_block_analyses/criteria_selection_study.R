## ---------------------------------------------------------------------------
## Model-selection criteria: do the BIC and ICL of PR #55 select the true q?
##
## PR #55 changed the criteria of both families (see NEWS, 0.4.0):
##   * penalty: one log(N_g) per group of parameters (alpha: log p; mean-block
##     cluster profiles and spherical noise: log(np); the rest: log n) instead
##     of log(n) for every parameter;
##   * ICL: BIC + 2 H(tau) only, as in Biernacki et al. (2000), instead of also
##     adding the entropy of the Gaussian latent W in the variance-block family.
## This study compares the old and the new definitions on simulated data where
## the true q is known, in both families, and splits the two changes apart:
##
##   criterion      penalty     entropy added to the BIC
##   BIC_old        log n       none
##   BIC_new        per group   none
##   ICL_old        log n       H(tau) + H(W) (var) / H(tau) (mean)
##   ICL_new        per group   H(tau)
##   ICL_oldpen     log n       H(tau)               (new entropy, old penalty)
##   ICL_exact      per group   none, but the exact log-likelihood of the MAP
##                              partition (known-clusters fit) + log p(C)
##   BIC/ICL_betalogn            as BIC/ICL_new, with the mean-block cluster
##                              profiles charged log n again (= *_new for var)
##
## All of them are recomputed from the same fits (deviance, nb_param,
## penalty, entropies), so the comparison is free of optimization noise.
##
## Optimization noise is controlled separately. Each dataset is fitted as a
## collection over a range of q around the truth, then refine()d; at the true
## q, a fit started from the true partition ("oracle start") replaces the
## collection's model when its deviance is lower. A criterion that misses the
## truth on the refined collection but finds it once the oracle start is in
## is failing because of the optimizer, not because of its penalty.
##
## Data are simulated from the models themselves, so they are fitted with
## scale = FALSE (well specified: column rescaling breaks the shared latent
## signal / shared mean profile within a block). scale = TRUE, the package
## default, is kept as a secondary arm.
##
## Kept in inst/mean_block_analyses/ for convenience (excluded from the build),
## although it covers both families. Run from the package root.
## ---------------------------------------------------------------------------

suppressMessages({
  if (!"normalblockr" %in% loadedNamespaces()) devtools::load_all(quiet = TRUE)
  library(parallel)
})

OUT     <- "inst/mean_block_analyses/criteria_selection_study.rds"
N_REP   <- as.integer(Sys.getenv("N_REP", 50))
N_CORES <- as.integer(Sys.getenv("N_CORES", 20))

## --- data-generating processes ----------------------------------------------

## variance-block: the package generator; `noise` scales the range of the
## idiosyncratic variances D relative to the default c(0.5, 1.5)
simulate_var <- function(n, p, q, noise) {
  ex <- generate_normal_block_var_data(n = n, p = p, d = 1, q = q, range_D = as.numeric(noise) * c(0.5, 1.5))
  list(Y = ex$Y, X = ex$X, cl = as.vector(ex$parameters$C %*% seq_len(q)))
}

## mean-block: the package generator's parameters, with the residual
## covariance either as generated (full, Erdos-Renyi precision: the default
## diagonal fit is then misspecified) or reduced to its diagonal
simulate_mean <- function(n, p, q, SNR, sigma) {
  X   <- matrix(runif(n * 2, 0, 10), n, 2)
  par <- normalblockr:::generate_normal_block_mean_param(X, p, q, SNR = SNR)
  Sigma <- if (sigma == "diagonal") diag(diag(par$Sigma)) else par$Sigma
  Y <- X %*% par$B %*% t(par$C) + MASS::mvrnorm(n, rep(0, p), Sigma)
  list(Y = Y, X = X, cl = as.vector(par$C %*% seq_len(q)))
}

## --- criteria components of one fit -----------------------------------------

components <- function(fit, truth_cl, family, dat) {
  tau <- fit$memberships
  H_tau <- -sum(xlogx(tau))
  H_W <- if (family == "var") {
    S <- fit$var_par$S
    .5 * fit$n * fit$q * log(2 * pi * exp(1)) + .5 * fit$n * sum(log(S))
  } else 0
  ## exact complete-data criterion at the MAP partition: W integrated out
  ## exactly by the known-clusters fit, C completed (Biernacki's ICL proper).
  ## Empty clusters of the MAP partition are dropped, so the alpha penalty
  ## counts the non-empty ones only.
  cl    <- fit$clustering
  known <- suppressWarnings(normal_block(dat, blocks = model.matrix(~ 0 + factor(cl)), model = family,
                                         control = NB_control(verbose = FALSE)))
  prop  <- table(cl) / length(cl)
  ICL_exact <- -2 * (known$loglik + sum(log(prop[as.character(cl)]))) +
    known$penalty + (length(prop) - 1) * log(fit$p)
  data.frame(q = fit$q, deviance = fit$deviance, nb_param = fit$nb_param,
             penalty = fit$penalty, H_tau = H_tau, H_W = H_W, ICL_exact = ICL_exact,
             ARI = aricode::ARI(cl, truth_cl))
}

criteria_from <- function(df, n) {
  within(df, {
    BIC_old    <- deviance + log(n) * nb_param
    BIC_new    <- deviance + penalty
    ICL_old    <- BIC_old + 2 * (H_tau + H_W)
    ICL_new    <- BIC_new + 2 * H_tau
    ICL_oldpen <- BIC_old + 2 * H_tau
  })
}

## --- one replicate ------------------------------------------------------------

one_run <- function(cfg, rep) {
  set.seed(1000 * cfg$id + rep)
  sim <- if (cfg$family == "var") simulate_var(cfg$n, cfg$p, cfg$q, cfg$noise)
         else simulate_mean(cfg$n, cfg$p, cfg$q, as.numeric(cfg$SNR), cfg$sigma)
  dat <- NormalBlockData$new(sim$Y, sim$X, scale = cfg$scale)
  q_range <- max(2, cfg$q - 3):(cfg$q + 4)
  ctrl <- NB_control(verbose = FALSE)
  co <- suppressWarnings(normal_block(dat, blocks = q_range, model = cfg$family, control = ctrl))
  raw <- do.call(rbind, lapply(co$models, components, sim$cl, cfg$family, dat))
  suppressWarnings(co$refine(verbose = FALSE))
  ref <- do.call(rbind, lapply(co$models, components, sim$cl, cfg$family, dat))
  oracle_fit <- suppressWarnings(normal_block(dat, blocks = cfg$q, model = cfg$family,
                                              control = NB_control(verbose = FALSE, clustering_init = sim$cl)))
  orc <- ref
  i <- which(orc$q == cfg$q)
  if (oracle_fit$deviance < orc$deviance[i]) orc[i, ] <- components(oracle_fit, sim$cl, cfg$family, dat)
  out <- rbind(cbind(stage = "raw", raw), cbind(stage = "refined", ref), cbind(stage = "oracle", orc))
  cbind(cfg[c("id", "family", "scale", "n", "p", "q_true", "noise", "SNR", "sigma")], rep = rep,
        criteria_from(out, cfg$n))
}

## --- design ---------------------------------------------------------------------
## "-" marks a factor that does not apply to the family (NA would be dropped by split())

grid_var  <- expand.grid(family = "var", n = c(30, 60, 120), p = 60, q = c(4, 8),
                         noise = c(1, 4), SNR = "-", sigma = "-", scale = c(FALSE, TRUE), stringsAsFactors = FALSE)
grid_mean <- expand.grid(family = "mean", n = c(30, 60, 120, 240), p = 60, q = c(3, 6),
                         noise = "-", SNR = c(1, 5, 25), sigma = c("diagonal", "full"), scale = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
design <- rbind(grid_var, grid_mean)
design$id <- seq_len(nrow(design))
design$q_true <- design$q

## --- run ------------------------------------------------------------------------
## OUT keeps one row per (dataset, stage): the q selected by each criterion and
## the ARI at the true q (the per-q table weighs ~5 MB; set OUT_FULL to keep it).

CRITS <- c("BIC_old", "BIC_new", "BIC_betalogn", "ICL_old", "ICL_oldpen", "ICL_new", "ICL_betalogn", "ICL_exact")

select_q <- function(res) {
  ## mean-block variant: cluster profiles charged log(n) again (alpha kept at
  ## log p), i.e. PR #55 without its log(np) for beta_k; d = 2 covariates
  m <- res$family == "mean"
  res$BIC_betalogn <- res$BIC_new - m * res$q * 2 * log(res$p)
  res$ICL_betalogn <- res$ICL_new - m * res$q * 2 * log(res$p)
  do.call(rbind, lapply(split(res, list(res$id, res$rep, res$stage), drop = TRUE), function(d) {
    sel <- sapply(CRITS, function(cr) d$q[which.min(d[[cr]])])
    cbind(d[1, c("id", "family", "scale", "n", "p", "q_true", "noise", "SNR", "sigma", "rep", "stage")],
          t(sel), ARI_true_q = d$ARI[d$q == d$q_true[1]])
  }))
}

if (!file.exists(OUT) || nzchar(Sys.getenv("RERUN"))) {
  jobs <- expand.grid(id = design$id, rep = seq_len(N_REP))
  t0 <- Sys.time()
  res <- mclapply(seq_len(nrow(jobs)), function(j) {
    tryCatch(one_run(as.list(design[jobs$id[j], ]), jobs$rep[j]),
             error = function(e) { message("job ", j, ": ", conditionMessage(e)); NULL })
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  res <- do.call(rbind, res)
  if (nzchar(Sys.getenv("OUT_FULL"))) saveRDS(res, Sys.getenv("OUT_FULL"))
  selected <- select_q(res)
  attr(selected, "elapsed") <- Sys.time() - t0
  saveRDS(selected, OUT)
}
selected <- readRDS(OUT)
cat("elapsed:", format(attr(selected, "elapsed")), "-- datasets:", nrow(selected) / 3, "\n")

## --- summaries ------------------------------------------------------------------

summarise <- function(sel) {
  key <- c("family", "scale", "n", "q_true", "noise", "SNR", "sigma", "stage")
  do.call(rbind, lapply(split(sel, sel[key], drop = TRUE), function(d) {
    hit  <- sapply(CRITS, function(cr) mean(d[[cr]] == d$q_true))
    bias <- sapply(CRITS, function(cr) mean(d[[cr]] - d$q_true))
    cbind(d[1, key], n_rep = nrow(d), ARI = round(mean(d$ARI_true_q), 2),
          t(round(100 * hit)), setNames(as.data.frame(t(round(bias, 2))), paste0("bias_", CRITS)))
  }))
}

tab <- summarise(selected)
tab <- tab[order(tab$family, tab$scale, tab$stage, tab$q_true, tab$noise, tab$SNR, tab$sigma, tab$n), ]
options(width = 250)
for (st in c("refined", "oracle", "raw")) {
  cat("\n==== stage:", st, " -- % of replicates selecting the true q (then mean selected q - true q)\n")
  print(tab[tab$stage == st, setdiff(names(tab), "stage")], row.names = FALSE)
}
cat("\n==== overall % exact, by family and stage\n")
print(aggregate(selected[CRITS] == selected$q_true, selected[c("family", "scale", "stage")], function(x) round(100 * mean(x))))

## by regime: is the true partition recoverable at the true q at all?
ref <- selected[selected$stage == "refined", ]
ref$regime <- ifelse(ref$ARI_true_q >= 0.9, "separable", "weak")
for (what in c("exact", "over", "under")) {
  cat("\n==== refined, % of datasets where the selected q is", what, "\n")
  f <- switch(what, exact = function(v) round(100 * mean(v == 0)),
              over = function(v) round(100 * mean(v > 0)), under = function(v) round(100 * mean(v < 0)))
  print(aggregate(ref[CRITS] - ref$q_true, ref[c("family", "scale", "regime")], f))
}
cat("\ndatasets per cell:\n"); print(with(ref, table(paste(family, scale), regime)))
