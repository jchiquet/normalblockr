## Generates results/best_of_inits_refine_benchmark.rds, used by the
## "best_of_inits() and refine(), collection-level" section of
## clustering_initialization_benchmark.qmd. Unlike the original per-(q,
## method) protocol (independent dummy fits, no refine()), this one builds a
## full collection per (dataset, method) via normal_block() itself -- refine()
## only applies to a whole collection (it needs adjacent q's) -- and records
## criteria both before and after a single refine() call, so the cost of
## building the collection is paid once, not duplicated across the
## with/without-refine comparison.
devtools::load_all(".", quiet = TRUE)

run_method <- function(data, q_values, method, zero_inflation, dataset_name) {
  ci <- if (method == "best_of_inits") "best_of_inits" else method
  t0 <- Sys.time()
  model <- normal_block(data, blocks = q_values, zero_inflation = zero_inflation,
                         control = NB_control(verbose = FALSE, clustering_init = ci))
  t_build <- as.numeric(Sys.time() - t0, units = "secs")
  before <- model$criteria; before$stage <- "before_refine"
  t1 <- Sys.time()
  model$refine(verbose = FALSE)
  t_refine <- as.numeric(Sys.time() - t1, units = "secs")
  after <- model$criteria; after$stage <- "after_refine"
  res <- rbind(before, after)
  res$method <- method; res$dataset <- dataset_name
  res$t_build <- t_build; res$t_refine <- t_refine
  res
}

data(brca_rppa)
data_brca <- NormalBlockData$new(
  brca_rppa$expr,
  model.matrix(~ 0 + PAM50_SUBTYPE, data = brca_rppa$covariates)
)
data(onema)
data_onema <- NormalBlockData$new(
  log(1 + onema$biomass),
  model.matrix(~ 1, data = onema$covariates)
)
data(university)
Y_web    <- log(1 + university$frequencies[, university$terms])
data_web <- NormalBlockData$new(Y_web, X = matrix(1, nrow(Y_web), 1))

q_lists  <- list(brca = 2:30, onema = 2:8, web = 2:15)
methods  <- c("kmeans", "ward2", "sbm", "spectral", "best_of_inits")
datasets <- list(brca = data_brca, onema = data_onema, web = data_web)
zi_flags <- c(brca = FALSE, onema = TRUE, web = FALSE)

benchmark <- do.call(rbind, lapply(names(datasets), function(nm) {
  do.call(rbind, lapply(methods, function(m) {
    cat(nm, m, "\n")
    run_method(datasets[[nm]], q_lists[[nm]], m, zi_flags[[nm]], nm)
  }))
}))

saveRDS(benchmark, "inst/clustering_initialization_benchmark/results/best_of_inits_refine_benchmark.rds")
