## Diagnostic script for the mean-block model on onema (real data, 2nd
## dataset -- cross-checks the brca_rppa findings against a different
## n/p/covariate structure before writing the vignette).
suppressPackageStartupMessages({
  library(normalblockr)
  library(aricode)
})

data(onema)
Y <- log1p(as.matrix(onema$biomass))
X <- model.matrix(~ scale(alt) + scale(slope) + scale(temperature_med) + scale(flow_med),
                  data = onema$covariates)
data_env <- NormalBlockData$new(Y, X)
cat("n =", nrow(Y), " p =", ncol(Y), " d =", ncol(X), "\n")
cat("rank(X) =", qr(X)$rank, "(full column rank:", qr(X)$rank == ncol(X), ")\n\n")

## ---------------------------------------------------------------------
## 1. Initialization comparison (q = 8, as for brca_rppa)
## ---------------------------------------------------------------------
cat("=== 1. clustering_init comparison (q = 8) ===\n")
for (h in c("kmeans", "ward2", "spectral")) {
  t0 <- Sys.time()
  m <- normal_block(data_env, blocks = 8, model = "mean",
                    control = NB_control(verbose = FALSE, clustering_init = h))
  cat(sprintf("%-9s ELBO: %10.2f  niter: %3d  time: %.1fs\n",
              h, m$loglik, m$criteria$niter, as.numeric(Sys.time() - t0, units = "secs")))
}
t0 <- Sys.time()
m_boi <- normal_block(data_env, blocks = 8, model = "mean",
                      control = NB_control(verbose = FALSE, clustering_init = "best_of_inits"))
cat(sprintf("%-9s ELBO: %10.2f  niter: %3d  time: %.1fs\n",
            "best_of_inits", m_boi$loglik, m_boi$criteria$niter, as.numeric(Sys.time() - t0, units = "secs")))

## ---------------------------------------------------------------------
## 2. Collection over q (default init) + refine()
## ---------------------------------------------------------------------
cat("\n=== 2. collection over q = 2:15 (default init) ===\n")
t0 <- Sys.time()
coll <- normal_block(data_env, blocks = 2:15, model = "mean", control = NB_control(verbose = FALSE))
cat("total time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
print(coll$criteria[, c("q", "loglik", "deviance", "BIC", "ICL", "niter")])
best_before <- coll$get_best_model("ICL")
cat("ICL selects q =", best_before$q, "\n")

cat("\n=== 3. refine() ===\n")
t0 <- Sys.time()
coll$refine(verbose = TRUE)
cat("refine() time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
print(coll$criteria[, c("q", "loglik", "deviance", "BIC", "ICL")])
best_after <- coll$get_best_model("ICL")
cat("ICL selects q =", best_after$q, "(was", best_before$q, "before refine())\n")

## ---------------------------------------------------------------------
## 4. Sanity: does the selected clustering track slope/altitude at all?
## (no species-level external annotation available for onema, unlike
## brca_rppa's GO terms -- station-level covariates are the closest proxy,
## but they are already regressed out via X, so a strong signal here would
## point at a modeling issue rather than validate the clustering)
## ---------------------------------------------------------------------
cat("\n=== 4. ELBO monotonicity on the selected model ===\n")
ll <- best_after$objective
cat("niter:", length(ll), " monotone:",
    if (length(ll) > 1) all(diff(ll) >= -1e-6) else NA, "\n")
cat("cluster sizes:", paste(table(best_after$clustering), collapse = ", "), "\n")
