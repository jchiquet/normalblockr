suppressMessages(devtools::load_all(".", quiet = TRUE))
library(aricode); RhpcBLASctl::blas_set_num_threads(1)
ctl <- function(...) NB_control(verbose = FALSE, ...)

## Etape 1 : clustering dans la moyenne (Sigma diagonale). Etape 2 : clustering
## dans la covariance des residus (modele Var, X reduit a l'intercept).
sequential <- function(Y, X, qs_mean, qs_var) {
  m  <- normal_block(NormalBlockData$new(Y, X), blocks = qs_mean, model = "mean", control = ctl())
  m1 <- if (inherits(m, "NormalBlockCollection")) m$get_best_model("ICL") else m
  R  <- Y - fitted(m1)
  v  <- normal_block(NormalBlockData$new(R, model.matrix(~ 1, data.frame(row.names = seq_len(nrow(Y))))),
                     blocks = qs_var, control = ctl())
  v1 <- if (inherits(v, "NormalBlockCollection")) v$get_best_model("ICL") else v
  list(mean = m1, var = v1)
}

cat("===== CONTROLE POSITIF : deux structures, partitions differentes =====\n")
set.seed(42)
n <- 300; p <- 100; d <- 3; qm <- 5; qv <- 4
X   <- cbind(1, matrix(rnorm(n * (d - 1)), n, d - 1))
Cm  <- normalblockr:::as_indicator(sample(rep(1:qm, length.out = p)))   # partition moyenne
Cv  <- normalblockr:::as_indicator(sample(rep(1:qv, length.out = p)))   # partition covariance
B   <- matrix(rnorm(d * qm, sd = 2), d, qm)
W   <- matrix(rnorm(n * qv, sd = 1.5), n, qv)                          # facteurs latents par bloc
Y   <- X %*% B %*% t(Cm) + W %*% t(Cv) + matrix(rnorm(n * p, sd = 0.5), n, p)
cat("ARI entre les deux vraies partitions :", round(ARI(get_clusters(Cm), get_clusters(Cv)), 3), "(par construction ~0)\n")
fit <- sequential(Y, X, seq(3, 15, by = 2), 1:8)
cat(sprintf("etape 1 : q_mu = %2d | ARI vs vraie partition moyenne     : %.3f\n",
            fit$mean$q, ARI(fit$mean$clustering, get_clusters(Cm))))
cat(sprintf("etape 2 : q_Sigma = %2d | ARI vs vraie partition covariance : %.3f\n",
            fit$var$q, ARI(fit$var$clustering, get_clusters(Cv))))
cat(sprintf("ARI entre les deux partitions estimees : %.3f\n\n", ARI(fit$mean$clustering, fit$var$clustering)))
flush.console()

cat("===== BRCA =====\n")
data(brca_rppa)
Yb <- as.matrix(brca_rppa$expr)
Xb <- model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
fb <- sequential(Yb, Xb, seq(5, 90, by = 5), 1:20)
direct <- normal_block(NormalBlockData$new(Yb, Xb), blocks = 1:20, control = ctl())$get_best_model("ICL")
cat(sprintf("etape 1 : q_mu = %d | etape 2 : q_Sigma = %d | Var direct : q = %d\n",
            fb$mean$q, fb$var$q, direct$q))
cat(sprintf("ARI(moyenne, covariance)              : %.3f\n", ARI(fb$mean$clustering, fb$var$clustering)))
cat(sprintf("ARI(covariance sequentielle, directe) : %.3f\n", ARI(fb$var$clustering, direct$clustering)))
