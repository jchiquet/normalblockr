## ---------------------------------------------------------------------------
## Should variance-block models rescale Y by default? A penalty-free check.
##
## On data simulated from the variance-block model, rescaling the columns of
## Y can only hurt: the generator gives every variable of a block the same
## loading on the block's latent signal, which rescaling by sqrt(Sigma_kk +
## d_j) destroys, and every criterion then over-selects q (see
## criteria_selection_study.R and the general vignette). On real data the
## question is empirical: variables of one block may carry the signal with
## different amplitudes, which the raw fit forbids and rescaling evens out.
## This compares the two fits by their cross-validated log-likelihood
## (normal_block_cv(), exact, on Y's original scale with the Jacobian of the
## rescaling included), paired over the observations, on the three real
## datasets of the initialization benchmark. Same folds for both scalings.
## Variance-block family only (despite the folder name); ~5 min on 6 cores.
## ---------------------------------------------------------------------------
suppressMessages({devtools::load_all(quiet = TRUE); library(parallel)})
data(brca_rppa); data(onema); data(university)
sets <- list(
  brca  = list(Y = as.matrix(brca_rppa$expr), X = model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates), q = 2:30),
  onema = list(Y = log1p(as.matrix(onema$biomass)), X = model.matrix(~ 1, data = onema$covariates), q = 2:15),
  web   = list(Y = log1p(university$frequencies[, university$terms]), X = matrix(1, nrow(university$frequencies), 1), q = 2:15))
jobs <- expand.grid(set = names(sets), scale = c(TRUE, FALSE), stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(j) {
  s <- sets[[jobs$set[j]]]
  set.seed(1)  # same folds for both scalings
  folds <- sample(rep_len(1:5, nrow(s$Y)))
  t0 <- Sys.time()
  cv <- normal_block_cv(NormalBlockData$new(s$Y, s$X, scale = jobs$scale[j]), blocks = s$q, folds = folds)
  list(cv = cv, time = as.numeric(Sys.time() - t0, units = "mins"))
}, mc.cores = 6, mc.preschedule = FALSE)
saveRDS(list(jobs = jobs, res = res), file.path(tempdir(), "var_scale_cv_real_data.rds"))
for (nm in names(sets)) {
  a <- res[[which(jobs$set == nm & jobs$scale)]]$cv; b <- res[[which(jobs$set == nm & !jobs$scale)]]$cv
  pa <- attr(a, "per_observation"); pb <- attr(b, "per_observation")
  ia <- which.max(a$loglik); ib <- which.max(b$loglik)
  ja <- match(attr(a, "best"), a$q); jb <- match(attr(b, "best"), b$q)
  d_max <- pa[, ia] - pb[, ib]; d_1se <- pa[, ja] - pb[, jb]
  cat(sprintf("%-6s scaled: max %.0f (q=%d), 1SE q=%d | raw: max %.0f (q=%d), 1SE q=%d\n", nm,
              a$loglik[ia], a$q[ia], attr(a, "best"), b$loglik[ib], b$q[ib], attr(b, "best")))
  cat(sprintf("       scaled - raw, paired over the %d observations: at the max %+.0f (se %.0f), at the 1SE choice %+.0f (se %.0f); obs. favouring scaled: %.0f%%\n",
              length(d_max), sum(d_max), sd(d_max) * sqrt(length(d_max)), sum(d_1se), sd(d_1se) * sqrt(length(d_1se)), 100 * mean(d_1se > 0)))
}
cat("minutes:", round(sapply(res, `[[`, "time"), 1), "\n")
## Results (2026-11-17), scaled - raw at the one-SE choice of each:
##   brca  +2,080 (se 124), 88% of observations favour rescaling
##   onema   +169 (se  65), 67%
##   web   -1,255 (se 754), 71% -- not significant, a few documents dominate
## Rescaling stays the variance-block default.
