## Diagnostic script for the mean-block model on brca_rppa (real data).
## Not a vignette: exercises initialization, convergence, model selection
## and refine() to surface problems/limits before writing one.
suppressPackageStartupMessages({
  library(normalblockr)
  library(aricode)
})

data(brca_rppa)
Y <- as.matrix(brca_rppa$expr)
X <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
data_subtype <- NormalBlockData$new(Y, X)
cat("n =", nrow(Y), " p =", ncol(Y), " d =", ncol(X), "\n")
cat("rank(X) =", qr(X)$rank, "(full column rank:", qr(X)$rank == ncol(X), ")\n\n")

## ---------------------------------------------------------------------
## 1. Initialization: does kmeans beat ward2 here too (as on simulated data)?
## ---------------------------------------------------------------------
cat("=== 1. clustering_init comparison (q = 8) ===\n")
for (h in c("kmeans", "ward2", "spectral")) {
  t0 <- Sys.time()
  m <- normal_block(data_subtype, blocks = 8, model = "mean",
                    control = NB_control(verbose = FALSE, clustering_init = h))
  cat(sprintf("%-9s ELBO: %10.2f  niter: %3d  time: %.1fs\n",
              h, m$loglik, m$criteria$niter, as.numeric(Sys.time() - t0, units = "secs")))
}
t0 <- Sys.time()
m_boi <- normal_block(data_subtype, blocks = 8, model = "mean",
                      control = NB_control(verbose = FALSE, clustering_init = "best_of_inits"))
cat(sprintf("%-9s ELBO: %10.2f  niter: %3d  time: %.1fs\n",
            "best_of_inits", m_boi$loglik, m_boi$criteria$niter, as.numeric(Sys.time() - t0, units = "secs")))

## ---------------------------------------------------------------------
## 2. Convergence diagnostics across a q range (default: kmeans init)
## ---------------------------------------------------------------------
cat("\n=== 2. collection over q = 2:20 (default init) ===\n")
t0 <- Sys.time()
coll <- normal_block(data_subtype, blocks = 2:20, model = "mean",
                     control = NB_control(verbose = FALSE))
cat("total time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
print(coll$criteria[, c("q", "nb_param", "loglik", "deviance", "BIC", "ICL", "niter")])

best_before <- coll$get_best_model("ICL")
cat("\nICL selects q =", best_before$q, "\n")

## ---------------------------------------------------------------------
## 3. refine(): does the neighbor-seeded search change the outcome?
## ---------------------------------------------------------------------
cat("\n=== 3. refine() ===\n")
t0 <- Sys.time()
coll$refine(verbose = TRUE)
cat("refine() time:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
print(coll$criteria[, c("q", "loglik", "deviance", "BIC", "ICL")])
best_after <- coll$get_best_model("ICL")
cat("\nICL selects q =", best_after$q, "(was", best_before$q, "before refine())\n")

## ---------------------------------------------------------------------
## 4. External validation: inferred clustering vs. GO annotation
## ---------------------------------------------------------------------
cat("\n=== 4. clustering vs. go_bp_term ===\n")
go <- brca_rppa$gene_annotation$go_bp_term
cat("distinct GO terms:", length(unique(go)), " (top 5 groups cover",
    sum(sort(table(go), decreasing = TRUE)[1:5]), "/", length(go), "proteins)\n")
go <- as.integer(factor(go)) # ARI()/sort_pairs() as.integer()'s character labels directly
cat("ARI(mean-block clustering, GO term):", round(ARI(best_after$clustering, go), 4), "\n")

## same comparison for the variance-block model, as a reference point
NB_var <- normal_block(data_subtype, blocks = 2:20, control = NB_control(verbose = FALSE))
NB_var$refine()
best_var <- NB_var$get_best_model("ICL")
cat("(reference) variance-block: q =", best_var$q,
    " ARI vs GO term:", round(ARI(best_var$clustering, go), 4), "\n")
cat("ARI(mean-block clustering, variance-block clustering):",
    round(ARI(best_after$clustering, best_var$clustering), 4), "\n")

## ---------------------------------------------------------------------
## 5. ELBO monotonicity check on the selected model (sanity check for the
## Gauss-Seidel tau update and the M-step, on real, non-simulated data)
## ---------------------------------------------------------------------
cat("\n=== 5. ELBO monotonicity on the selected model ===\n")
ll <- best_after$objective
cat("niter:", length(ll), " monotone:",
    if (length(ll) > 1) all(diff(ll) >= -1e-6) else NA, "\n")
