## Diagnostic script for the mean-block model at larger scale (simulated):
## does refine()'s cost/benefit and the kmeans-over-ward2 initialization
## finding (brca_rppa, onema) hold up as n/p/q grow, and where does the
## wall-clock time actually go?
suppressPackageStartupMessages({
  library(normalblockr)
  library(aricode)
})

set.seed(1)
n <- 300; p <- 150; d <- 5; q_true <- 10
ex <- generate_normal_block_mean_data(n = n, p = p, d = d, q = q_true, SNR = 2)
data_sim <- NormalBlockData$new(ex$Y, ex$X, scale = FALSE)
vrai <- get_clusters(ex$parameters$C)
cat("n =", n, " p =", p, " d =", d, " q_true =", q_true, "\n\n")

## ---------------------------------------------------------------------
## 1. Initialization comparison (q = q_true)
## ---------------------------------------------------------------------
cat("=== 1. clustering_init comparison (q =", q_true, ") ===\n")
for (h in c("kmeans", "ward2", "spectral")) {
  t0 <- Sys.time()
  m <- normal_block(data_sim, blocks = q_true, model = "mean",
                    control = NB_control(verbose = FALSE, clustering_init = h))
  cat(sprintf("%-9s ELBO: %10.2f  ARI: %.3f  niter: %3d  time: %.1fs\n",
              h, m$loglik, ARI(m$clustering, vrai), m$criteria$niter,
              as.numeric(Sys.time() - t0, units = "secs")))
}
t0 <- Sys.time()
m_boi <- normal_block(data_sim, blocks = q_true, model = "mean",
                      control = NB_control(verbose = FALSE, clustering_init = "best_of_inits"))
cat(sprintf("%-9s ELBO: %10.2f  ARI: %.3f  niter: %3d  time: %.1fs\n",
            "best_of_inits", m_boi$loglik, ARI(m_boi$clustering, vrai), m_boi$criteria$niter,
            as.numeric(Sys.time() - t0, units = "secs")))

## ---------------------------------------------------------------------
## 2. Collection over q + refine(), with timing broken down
## ---------------------------------------------------------------------
cat("\n=== 2. collection over q = 5:15 (default init) ===\n")
t0 <- Sys.time()
coll <- normal_block(data_sim, blocks = 5:15, model = "mean", control = NB_control(verbose = FALSE))
t_coll <- as.numeric(Sys.time() - t0, units = "secs")
cat("collection fit time:", round(t_coll, 1), "s\n")
best_before <- coll$get_best_model("ICL")
cat("ICL selects q =", best_before$q, " ARI vs truth:", round(ARI(best_before$clustering, vrai), 3), "\n")

cat("\n=== 3. refine() ===\n")
t0 <- Sys.time()
coll$refine(verbose = FALSE)
t_refine <- as.numeric(Sys.time() - t0, units = "secs")
cat("refine() time:", round(t_refine, 1), "s (", round(t_refine / t_coll, 1), "x the collection fit)\n")
best_after <- coll$get_best_model("ICL")
cat("ICL selects q =", best_after$q, " ARI vs truth:", round(ARI(best_after$clustering, vrai), 3),
    "(was q =", best_before$q, "ARI", round(ARI(best_before$clustering, vrai), 3), "before refine())\n")

## ---------------------------------------------------------------------
## 4. Where does the time go? profile a single fit at q_true
## ---------------------------------------------------------------------
cat("\n=== 4. profiling a single fit (q =", q_true, ", kmeans init) ===\n")
Rprof(tmp <- tempfile())
invisible(normal_block(data_sim, blocks = q_true, model = "mean", control = NB_control(verbose = FALSE)))
Rprof(NULL)
prof <- summaryRprof(tmp)$by.self
print(head(prof[order(-prof$self.time), ], 10))
