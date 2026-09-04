suppressMessages({library(igraph); library(MASS); library(dplyr); devtools::load_all(quiet=TRUE)})
source("inst/normalblockmean/code_nestor/01_generation_donnees.R")
set.seed(42)
res <- list()
for (snr in c(0.3, 1.5)) for (qv in c(3, 5)) for (n in c(80, 150, 220, 300)) for (rep in 1:10) {
  dd <- generer_donnees(n = n, p = 60, d = 2, q = qv, snr = snr)
  data <- NormalBlockData$new(dd$Z, dd$X, scale = FALSE)
  f <- suppressWarnings(normal_block(data, blocks = 2:6, model = "mean",
        control = NB_control(verbose = FALSE, noise_covariance = "full")))
  res[[length(res)+1]] <- data.frame(snr = snr, q = qv, n = n,
    q_BIC = f$get_best_model("BIC")$q, q_ICL = f$get_best_model("ICL")$q)
}
r <- bind_rows(res); saveRDS(r, "/tmp/claude-1001/-home-jchiquet-Repos-jeanne-normalblockr/fd1c5035-aedc-4039-9568-ad44954786fd/scratchpad/qsel_nopen.rds")
cat("=== selection de q, Sigma pleine NON penalisee (nu identique a celui de Nestor) ===\n")
print(as.data.frame(r |> group_by(snr, q) |> summarise(
  pct_BIC = 100*mean(q_BIC == q), pct_ICL = 100*mean(q_ICL == q),
  sous = 100*mean(q_BIC < q), sur = 100*mean(q_BIC > q), n = n(), .groups="drop")), digits = 3)
