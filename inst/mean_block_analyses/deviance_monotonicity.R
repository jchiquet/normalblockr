suppressMessages(devtools::load_all(".", quiet = TRUE))
data(brca_rppa)
d <- NormalBlockData$new(as.matrix(brca_rppa$expr),
                         model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates))
diag_report <- function(dev, q) {
  inc <- diff(dev)                      # should be <= 0 everywhere
  viol <- which(inc > 1e-6)
  sprintf("%2d/%2d violations | pire remontee %8.1f | deviance q=1 -> q=%d : %.0f -> %.0f",
          length(viol), length(inc), if (length(viol)) max(inc[viol]) else 0,
          max(q), dev[1], dev[length(dev)])
}
for (nc in c("full", "diagonal", "spherical")) {
  set.seed(1)
  co <- normal_block(d, blocks = 1:15, model = "mean",
                     control = NB_control(verbose = FALSE, noise_covariance = nc))
  cat(sprintf("%-10s avant refine : %s\n", nc, diag_report(co$criteria$deviance, co$q_list)))
  co$refine(verbose = FALSE)
  cat(sprintf("%-10s apres refine : %s\n", nc, diag_report(co$criteria$deviance, co$q_list)))
  flush.console()
}
