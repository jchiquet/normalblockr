## Which split() rule should refine() use by default for the mean-block family?
## "noise" (median cut of the noise variances) vs "kmeans" (2-means on each
## variable's least-squares trajectory X B_ols). For each seed, the same
## initial collection is refined with each rule (serialize/unserialize copy),
## so the two rules start from the very same models.
## Source from the package root after devtools::load_all().
copy <- function(x) unserialize(serialize(x, NULL))
one_run <- function(label, d, blocks, seed, zi = FALSE, truth = NULL, q_true = NA) {
  set.seed(seed)
  co <- suppressWarnings(normal_block(d, blocks = blocks, model = "mean", zero_inflation = zi,
                                      control = NB_control(verbose = FALSE)))
  dev0 <- sum(co$criteria$deviance)
  do.call(rbind, lapply(c("noise", "kmeans"), function(rule) {
    cr <- copy(co); set.seed(seed + 1000)
    suppressWarnings(cr$refine(verbose = FALSE, split_rule = rule))
    crit <- cr$criteria
    data.frame(data = label, seed = seed, rule = rule,
               gain  = dev0 - sum(crit$deviance),          # deviance removed by refine(), summed over q
               viol  = sum(diff(crit$deviance) > 1e-6),     # monotonicity violations left
               q_icl = cr$q_list[which.min(crit$ICL)],
               ari_qtrue = if (is.null(truth)) NA else aricode::ARI(cr$get_model(q_true)$clustering, truth))
  }))
}

data(brca_rppa)
d_brca <- NormalBlockData$new(as.matrix(brca_rppa$expr),
                              model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates))
data(onema)
d_onema <- NormalBlockData$new(log(1 + onema$biomass),
                               model.matrix(~ 1 + temperature_med, data = onema$covariates))
res <- rbind(
  do.call(rbind, lapply(1:5, function(s) one_run("brca", d_brca, 1:40, s))),
  do.call(rbind, lapply(1:3, function(s) one_run("onema (ZI)", d_onema, 1:15, s, zi = TRUE))),
  do.call(rbind, lapply(1:10, function(s) {
    set.seed(500 + s)
    ex <- generate_normal_block_mean_data(n = 120, p = 60, d = 3, q = 6)
    one_run("simulation", NormalBlockData$new(ex$Y, ex$X), 1:12, s,
            truth = apply(ex$parameters$C, 1, which.max), q_true = 6)
  }))
)
print(res, digits = 4)
## Results (2026-11-17, unscaled data -- the mean-block default since #60):
## "noise" removes more deviance on every brca seed (2-15% more) and every
## onema seed (about 8%); the simulations tie 7/10, "kmeans" is ahead 2/10
## and "noise" 1/10. Selected q and ARI at the true q are identical in every
## run. Hence "noise" stays the mean-block default.
## The previous run (2026-11-16, on rescaled data) reached the same conclusion
## with smaller margins (~1.7% on brca, <1% on onema).
