library(parallel)
library(tibble)
library(normalblockr)

source("data_generation.R")
source("omega_generation.R")
source("parameters_generation.R")
source("measures.R")
source("utils.R")


#' @description runs one simulation for a given set of parameters
#' @returns a dataframe with the measures for each model: integrated inference
#' with fixed blocks, fixed q and baseline inference procedures
#' @param simu rank of the simulation among the set of simulations it is run in
#' @param n number of rows in the abundance matrix
#' @param p number of columns in the abundance matrix
#' @param q number of clusters in the model
#' @param d number of covariates
#' @param kappa vector of zero-inflation probabilities
#' @param omega_structure graph structure for the precision matrix
#' @param total_simu optional, number of simulations in the set this one simulation
#' is part of
#' @param saving_file optional, path to saving file to save the results from that
#' one simulation
#' @param verbose whether information should be given throughout the simulation
one_simu_sparse_fixed_q <- function(simu, n, p, q, d, kappa,
                                    omega_structure = c("erdos_renyi",
                                                        "preferential_attachment",
                                                        "community"),
                                    total_simu = NA, saving_file = NULL,
                                    verbose = 1){

  omega_structure <- match.arg(omega_structure)
  ## Show progression
  if (verbose >= 1) cat("Simu", simu, "of", total_simu)

  ## Draw parameters and data according to the Normal-Block model
  X <- generate_X(n, d)
  true_param <- generate_parameters(X, p, q, kappa, omega_structure)
  data <- NormalBlockData$new(Y = normal_block_data(true_param, X), X = X)

  ## Arguments common to all function call
  args_common <- list(data,
                      sparsity = TRUE,
                      zero_inflation = (min(kappa) > 0),
                      control = NB_control(verbose = FALSE))

  ## --------------------------------------------------------------
  ## NORMAL BLOCKS WITH KNOWN FIXED BLOCKS
  ##
  if (verbose > 1) cat("\n V-EM Known fixed blocks")

  args_fixed_blocks <- c(args_common, list(blocks = true_param$C))

  t0 <- Sys.time()
  myNB_fixed_blocks <- do.call(normal_block, args_fixed_blocks)
  t  <- Sys.time() - t0
  ref_penalty <- myNB_fixed_blocks$get_best_model("BIC")$sparsity


  res_fixed_blocks_BIC   <- get_measures(myNB_fixed_blocks, true_param, "BIC",
                                         fixed_blocks = TRUE)
  res_fixed_blocks_BIC   <- c(res_fixed_blocks_BIC, "time" = as.numeric(t))

  res_fixed_blocks_StARS   <- get_measures(myNB_fixed_blocks, true_param, "StARS",
                                         fixed_blocks = TRUE)
  res_fixed_blocks_StARS   <- c(res_fixed_blocks_StARS, "time" = as.numeric(t))


  ## --------------------------------------------------------------
  ## NORMAL BLOCKS WITH UNKNOWN BLOCKS BUT KNOWN NUMBER OF BLOCKS
  ##
  if (verbose > 1) cat("\n V_EM Unknown blocks")
  args_fixed_q <- c(args_common, list(blocks = q))
  t0 <- Sys.time()
  myNB_fixed_q <- do.call(normal_block, args_fixed_q)
  t  <- Sys.time() - t0

  res_fixed_q_BIC   <- get_measures(myNB_fixed_q, true_param, "BIC")
  res_fixed_q_BIC   <- c(res_fixed_q_BIC, "time" = as.numeric(t))

  res_fixed_q_StARS   <- get_measures(myNB_fixed_q, true_param, "StARS")
  res_fixed_q_StARS   <- c(res_fixed_q_StARS, "time" = as.numeric(t))

  ## --------------------------------------------------------------
  ## G-LASSO ON RESIDUALS OF NORMAL BLOCKS WITH KNOWN BLOCKS
  ##
  if (verbose > 1) cat("\n 2-step Known blocks")
  args_baseline_fixed_blocks <- args_fixed_blocks
  args_baseline_fixed_blocks$control <-  NB_control(verbose = FALSE,
                                                    heuristic = TRUE)
  t0 <- Sys.time()
  baseline_fixed_blocks <- do.call(normal_block, args_baseline_fixed_blocks)
  t  <- Sys.time() - t0

  res_baseline_fixed_blocks_ref_penalty <- get_measures(baseline_fixed_blocks, true_param,
                                                        ref_penalty, fixed_blocks = TRUE)
  res_baseline_fixed_blocks_ref_penalty   <- c(res_baseline_fixed_blocks_ref_penalty, "time" = as.numeric(t))

  ## --------------------------------------------------------------
  ## G-LASSO ON RESIDUALS OF NORMAL BLOCKS WITH CLUSTERING ON SIGMA
  ##
  if(q <= 10){ # otherwise the sbm takes too long
    if (verbose > 1) cat("\n 2-step Unknown blocks (cluster the variance)")
    args_baseline_sigma_clustering <- args_fixed_q
    args_baseline_sigma_clustering$control <- NB_control(verbose = FALSE,
                                                         heuristic = TRUE,
                                                         clustering_init = "sbm")
    t0 <- Sys.time()
    baseline_sigma_clustering <- do.call(normal_block, args_baseline_sigma_clustering)
    t  <- Sys.time() - t0

    res_baseline_sigma_clustering_ref_penalty   <- get_measures(baseline_sigma_clustering, true_param,
                                                                ref_penalty)
    res_baseline_sigma_clustering_ref_penalty   <- c(res_baseline_sigma_clustering_ref_penalty, "time" = as.numeric(t))
  }


  ## --------------------------------------------------------------
  ## G-LASSO ON RESIDUALS OF NORMAL BLOCKS WITH CLUSTERING ON RESIDUALS
  ##
  if (verbose > 1) cat("\n 2-step Unknown blocks (cluster the variance)")
  args_baseline_residuals_clustering         <- args_fixed_q
  args_baseline_residuals_clustering$control <- NB_control(verbose = FALSE,
                                                           heuristic = TRUE,
                                                           clustering_init = "kmeans")

  t0 <- Sys.time()
  baseline_residuals_clustering <- do.call(normal_block, args_baseline_residuals_clustering)
  t  <- Sys.time() - t0

  res_baseline_residuals_clustering_ref_penalty <- get_measures(baseline_residuals_clustering,
                                                                true_param, ref_penalty)
  res_baseline_residuals_clustering_ref_penalty   <- c(res_baseline_residuals_clustering_ref_penalty,
                                                       "time" = as.numeric(t))

  if(q <= 10){
    res <- as.data.frame(
      cbind(simu = simu, n = n, p = p, q = q, zero_inflation = mean(kappa), omega_structure = omega_structure,
            rbind(c(method = "NB", res_fixed_blocks_BIC),
                  c(method = "NB", res_fixed_blocks_StARS),
                  c(method = "NB", res_fixed_q_BIC),
                  c(method = "NB", res_fixed_q_StARS),
                  c(method = "BL1", res_baseline_fixed_blocks_ref_penalty),
                  c(method = "BL2", res_baseline_sigma_clustering_ref_penalty),
                  c(method = "BL3", res_baseline_residuals_clustering_ref_penalty)
            )
      )
    )
  }else{
    res <- as.data.frame(
      cbind(simu = simu, n = n, p = p, q = q, zero_inflation = mean(kappa), omega_structure = omega_structure,
            rbind(c(method = "NB", res_fixed_blocks_BIC),
                  c(method = "NB", res_fixed_blocks_StARS),
                  c(method = "NB", res_fixed_q_BIC),
                  c(method = "NB", res_fixed_q_StARS),
                  c(method = "BL1", res_baseline_fixed_blocks_ref_penalty),
                  c(method = "BL3", res_baseline_residuals_clustering_ref_penalty)
            )
      )
    )
  }

  if (!is.null(saving_file)) write.csv(res, saving_file)
  if (verbose > 1) cat("\ndone \r")

  res
}


#' @description runs multiple simulations with one set of parameters
#' @returns a dataframe with the measures for each simulation
#' @param n_simu number of simulations to run
#' @param n number of rows in the abundance matrix
#' @param p number of columns in the abundance matrix
#' @param q number of clusters in the model
#' @param d number of covariates
#' @param kappa vector of zero-inflation probabilities
#' @param saving_file_prefix  optional, prefix to saving file if the results should
#' be saved
#' @param mc.cores number of cores to run the simulations on
multiple_simulations <-
  function(n_simu, n, p, q, d, kappa,
           omega_structure = c("erdos_renyi", "preferential_attachment","community"),
           saving_file_prefix = NULL, mc.cores = 1){

    cat("Settings: (n, p, q, d, kappa) = (",n, p, q, d, mean(kappa), ")\n")
    omega_structure <- match.arg(omega_structure)
    multiple_res <- mclapply(1:n_simu,
                             one_simu_sparse_fixed_q,
                             n, p, q, d, kappa, omega_structure,
                             mc.cores = mc.cores, total_simu = n_simu)
    res <- do.call(rbind, multiple_res)
    res
  }

#' @description runs multiple simulations for multiple sets of parameters
#' @returns a dataframe with the measures for each model: integrated inference
#' with fixed blocks, fixed q and baseline inference procedures
#' @param n_simu number of simulations to run for each set of parameters
#' @param n_list list of values for n, number of rows in the abundance matrix
#' @param p_list list of values for p, number of columns in the abundance matrix
#' @param q_list list of values for q number of clusters in the model
#' @param d_list list of values for d, number of covariates
#' @param kappa_list p_list list of values for kappa, vector of zero-inflation probabilities
#' @param saving_folder  optional, folder to save the results in
#' @param mc.cores number of cores to run the simulations on
grid_simulations <- function(n_simu, n_list, p_list, q_list, d_list, kappa_list,
                             omega_structure_list, saving_folder = NULL,
                             mc.cores = 1){

  settings <- expand.grid(n = n_list,
                          p = p_list,
                          d = d_list,
                          q = q_list,
                          kappa = kappa_list,
                          omega_structure = omega_structure_list,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>% as_tibble()

  settings$n_simu <- n_simu

  final_res <- purrr::pmap(settings, multiple_simulations, mc.cores = mc.cores)
  final_res <- do.call(rbind, final_res) %>% as_tibble()
  final_res
}

