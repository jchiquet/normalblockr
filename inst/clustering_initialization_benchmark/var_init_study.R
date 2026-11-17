## ---------------------------------------------------------------------------
## Which clustering_init should the variance-block family use by default?
##
## The mean-block default is "kmeans" (Nestor's protocol and the brca/onema
## diagnostics, inst/mean_block_analyses/); the variance-block default is
## "ward2", chosen on cold-started collections before refine() was reworked
## (#56) and before the criteria of #55. The re-run of the best_of_inits +
## refine() benchmark (#62) found kmeans + refine() ahead. This study decides
## whether the variance-block default can become "kmeans" too, so that both
## families share one.
##
## Part 1, simulations from the variance-block model (the package generator),
## with known truth: n in 30/60/120, p in 60/120, q in 4/8, two noise levels,
## scale = FALSE (well specified) and TRUE (the variance-block default), 20
## datasets per cell. On each, a collection over max(2, q-3):(q+4) per
## heuristic, recorded before and after refine() (refine() is not called by
## default, NB_control(refine = FALSE), so both matter).
## Part 2, real data: brca_rppa, onema (zero-inflated) and university, the
## datasets of clustering_initialization_benchmark.qmd, 10 seeds per
## heuristic (kmeans, spectral and sbm draw random starts; ward2 does not).
##
## Every job runs in its own process with one BLAS thread: run with
## OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1 from the package root.
## ---------------------------------------------------------------------------

suppressMessages({
  if (!"normalblockr" %in% loadedNamespaces()) devtools::load_all(quiet = TRUE)
  library(parallel)
})

DIR     <- Sys.getenv("OUT_DIR", "inst/clustering_initialization_benchmark/results")
OUT_SIM <- file.path(DIR, "var_init_study_simulations.rds")
OUT_REAL <- file.path(DIR, "var_init_study_real_data.rds")
N_REP   <- as.integer(Sys.getenv("N_REP", 20))
N_CORES <- as.integer(Sys.getenv("N_CORES", 20))
METHODS <- c("kmeans", "ward2", "spectral", "sbm")

## one collection per heuristic, before and after refine()
fit_collection <- function(dat, q_range, method, zi = FALSE, truth = NULL) {
  t0 <- Sys.time()
  co <- suppressWarnings(normal_block(dat, blocks = q_range, zero_inflation = zi,
                                      control = NB_control(verbose = FALSE, clustering_init = method)))
  t_fit <- as.numeric(Sys.time() - t0, units = "secs")
  grab <- function(stage, t) {
    cr <- co$criteria
    data.frame(method = method, stage = stage, q = cr$q, deviance = cr$deviance,
               BIC = cr$BIC, ICL = cr$ICL, time = t,
               n_viol = length(co$deviance_violations),
               ARI = if (is.null(truth)) NA_real_ else
                 sapply(co$models, function(m) aricode::ARI(m$clustering, truth)))
  }
  raw <- grab("raw", t_fit)
  t1 <- Sys.time()
  suppressWarnings(co$refine(verbose = FALSE))
  rbind(raw, grab("refined", t_fit + as.numeric(Sys.time() - t1, units = "secs")))
}

## --- part 1: simulations ------------------------------------------------------

design <- expand.grid(n = c(30, 60, 120), p = c(60, 120), q = c(4, 8), noise = c(1, 4),
                      scale = c(FALSE, TRUE))
design$id <- seq_len(nrow(design))

one_sim <- function(cfg, rep) {
  set.seed(1000 * cfg$id + rep)
  ex <- generate_normal_block_var_data(n = cfg$n, p = cfg$p, d = 1, q = cfg$q,
                                       range_D = cfg$noise * c(0.5, 1.5))
  truth <- as.vector(ex$parameters$C %*% seq_len(cfg$q))
  dat <- NormalBlockData$new(ex$Y, ex$X, scale = cfg$scale)
  q_range <- max(2, cfg$q - 3):(cfg$q + 4)
  res <- do.call(rbind, lapply(METHODS, function(m) {
    set.seed(rep)
    fit_collection(dat, q_range, m, truth = truth)
  }))
  cbind(cfg[c("id", "n", "p", "noise", "scale")], q_true = cfg$q, rep = rep, res)
}

if (!file.exists(OUT_SIM) || nzchar(Sys.getenv("RERUN"))) {
  jobs <- expand.grid(id = design$id, rep = seq_len(N_REP))
  t0 <- Sys.time()
  sim <- mclapply(seq_len(nrow(jobs)), function(j)
    tryCatch(one_sim(design[jobs$id[j], ], jobs$rep[j]),
             error = function(e) { message("job ", j, ": ", conditionMessage(e)); NULL }),
    mc.cores = N_CORES, mc.preschedule = FALSE)
  sim <- do.call(rbind, sim)
  attr(sim, "elapsed") <- Sys.time() - t0
  saveRDS(sim, OUT_SIM)
}

## --- part 2: real data ----------------------------------------------------------

if (!file.exists(OUT_REAL) || nzchar(Sys.getenv("RERUN"))) {
  data(brca_rppa); data(onema); data(university)
  datasets <- list(
    brca  = list(data = NormalBlockData$new(brca_rppa$expr,
                   model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)),
                 q = 2:30, zi = FALSE),
    onema = list(data = NormalBlockData$new(log(1 + onema$biomass),
                   model.matrix(~ 1, data = onema$covariates)),
                 q = 2:8, zi = TRUE),
    web   = list(data = NormalBlockData$new(log(1 + university$frequencies[, university$terms]),
                   X = matrix(1, nrow(university$frequencies), 1)),
                 q = 2:15, zi = FALSE))
  jobs <- expand.grid(seed = 1:10, method = METHODS, dataset = names(datasets),
                      stringsAsFactors = FALSE)
  jobs <- jobs[jobs$method != "ward2" | jobs$seed == 1, ]   # ward2 is deterministic
  jobs <- jobs[jobs$seed <= as.integer(Sys.getenv("N_SEEDS", 10)), ]
  t0 <- Sys.time()
  real <- mclapply(seq_len(nrow(jobs)), function(j) {
    ds <- datasets[[jobs$dataset[j]]]
    set.seed(jobs$seed[j])
    cbind(dataset = jobs$dataset[j], seed = jobs$seed[j],
          fit_collection(ds$data, ds$q, jobs$method[j], zi = ds$zi))
  }, mc.cores = N_CORES, mc.preschedule = FALSE)
  real <- do.call(rbind, real)
  attr(real, "elapsed") <- Sys.time() - t0
  saveRDS(real, OUT_REAL)
}

## --- summaries ------------------------------------------------------------------

sim  <- readRDS(OUT_SIM)
real <- readRDS(OUT_REAL)
options(width = 200)
cat("simulations:", format(attr(sim, "elapsed")), "-- real data:", format(attr(real, "elapsed")), "\n")

## per dataset and stage: rank of each method's deviance at every q (1 = best),
## ARI at the true q, selected q (ICL) and its ARI
per_dataset <- do.call(rbind, lapply(split(sim, list(sim$id, sim$rep, sim$stage), drop = TRUE), function(d) {
  d$rank <- ave(d$deviance, d$q, FUN = function(x) rank(x, ties.method = "min"))
  do.call(rbind, lapply(split(d, d$method), function(m) {
    i_sel <- which.min(m$ICL); i_true <- which(m$q == m$q_true[1])
    data.frame(m[1, c("id", "n", "p", "noise", "scale", "q_true", "rep", "stage", "method")],
               mean_rank = mean(m$rank), ARI_true = m$ARI[i_true],
               exact_q = m$q[i_sel] == m$q_true[1], ARI_sel = m$ARI[i_sel],
               n_viol = m$n_viol[1], time = m$time[1])
  }))
}))

agg <- function(d, by) aggregate(cbind(mean_rank, ARI_true, exact_q, ARI_sel, n_viol, time) ~ .,
                                 data = d[c(by, "mean_rank", "ARI_true", "exact_q", "ARI_sel", "n_viol", "time")],
                                 FUN = function(x) round(mean(x), 3))
cat("\n==== simulations, by scale, stage and method\n")
print(agg(per_dataset, c("scale", "stage", "method")), row.names = FALSE)
cat("\n==== simulations, scale = TRUE (default), by q, noise and method\n")
print(agg(subset(per_dataset, scale), c("stage", "q_true", "noise", "method")), row.names = FALSE)

cat("\n==== real data: deviance rank per q, total deviance, violations, ICL choice, time (mean over seeds)\n")
real$key <- paste(real$dataset, real$stage, real$q)
best <- tapply(real$deviance, real$key, min)
real$gap <- real$deviance - best[real$key]
rs <- do.call(rbind, lapply(split(real, list(real$dataset, real$method, real$stage), drop = TRUE), function(d) {
  per_seed <- do.call(rbind, lapply(split(d, d$seed), function(s)
    data.frame(tot_gap = sum(s$gap), n_viol = s$n_viol[1], q_ICL = s$q[which.min(s$ICL)], time = s$time[1])))
  data.frame(dataset = d$dataset[1], stage = d$stage[1], method = d$method[1], seeds = nrow(per_seed),
             gap_to_best = round(mean(per_seed$tot_gap)), gap_sd = round(sd(per_seed$tot_gap)),
             n_viol = round(mean(per_seed$n_viol), 1),
             q_ICL = paste(sort(unique(per_seed$q_ICL)), collapse = "/"),
             time_s = round(mean(per_seed$time), 1))
}))
print(rs[order(rs$dataset, rs$stage, rs$gap_to_best), ], row.names = FALSE)
