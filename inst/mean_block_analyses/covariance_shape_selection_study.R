suppressMessages(devtools::load_all(".", quiet = TRUE))
library(aricode); RhpcBLASctl::blas_set_num_threads(1)
q_true <- 5; n <- 200; qs <- 1:10; nseed <- 12
res <- NULL
for (p in c(60, 150)) for (gen in c("Sigma pleine", "Sigma diagonale")) for (s in seq_len(nseed)) {
  set.seed(1000 * p + s)
  ex <- generate_normal_block_mean_data(n = n, p = p, d = 3, q = q_true, SNR = 2)
  Y  <- if (gen == "Sigma pleine") ex$Y else {
    ## meme moyenne, mais bruit independant entre variables (variances conservees)
    mu <- ex$X %*% ex$parameters$B %*% t(ex$parameters$C)
    mu + matrix(rnorm(n * p, sd = rep(sqrt(diag(ex$parameters$Sigma)), each = n)), n, p)
  }
  dat  <- NormalBlockData$new(Y, ex$X, scale = FALSE)
  vrai <- get_clusters(ex$parameters$C)
  for (nc in c("full", "diagonal")) {
    co <- tryCatch(normal_block(dat, blocks = qs, model = "mean",
                     control = NB_control(verbose = FALSE, noise_covariance = nc)),
                   error = function(e) NULL)
    if (is.null(co)) next
    crit <- co$criteria
    res <- rbind(res, data.frame(
      p = p, gen = gen, seed = s, forme = nc,
      q_BIC = crit$q[which.min(crit$BIC)], q_ICL = crit$q[which.min(crit$ICL)],
      ARI_qtrue = ARI(co$get_model(q_true)$clustering, vrai),
      violations = sum(diff(crit$deviance) > 1e-6)))
  }
  cat("."); flush.console()
}
cat("\n\n")
saveRDS(res, "/tmp/claude-1001/-home-jchiquet-Repos-jeanne-normalblockr/fd1c5035-aedc-4039-9568-ad44954786fd/scratchpad/study.rds")
summ <- do.call(rbind, lapply(split(res, list(res$p, res$gen, res$forme)), function(z) data.frame(
  p = z$p[1], generatif = z$gen[1], ajuste = z$forme[1],
  BIC_exact = sprintf("%d/%d", sum(z$q_BIC == q_true), nrow(z)),
  ICL_exact = sprintf("%d/%d", sum(z$q_ICL == q_true), nrow(z)),
  ecart_BIC = round(mean(abs(z$q_BIC - q_true)), 2),
  ecart_ICL = round(mean(abs(z$q_ICL - q_true)), 2),
  ARI = round(mean(z$ARI_qtrue), 3),
  viol = round(mean(z$violations), 1))))
print(summ[order(summ$p, summ$generatif, summ$ajuste), ], row.names = FALSE)
