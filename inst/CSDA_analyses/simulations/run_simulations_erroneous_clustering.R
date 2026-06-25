############### Loading useful libraries and external files ####################
library(parallel)
library(tibble)
library(normalblockr)

source("data_generation.R")
source("omega_generation.R")
source("parameters_generation.R")
source("measures.R")
source("utils.R")

############### Functions to assess Normal-Block with an erroneous clustering ##

#' @description introduces errors in a one-hot encoded clustering matrix
#' @returns a one-hot encoded clustering matrix with errors compared to the input
#' @param C one-hot-encoded clustering matrix
#' @param error_rate percentage of errors to introduce in the clustering
introduce_clustering_errors <- function(C, error_rate = 0.15){
  q <- ncol(C) ; p = nrow(C)
  rows_to_modify <- sample.int(p, round(error_rate * p))
  new_cols <- unlist(lapply(rows_to_modify, f <- function(x){
    cl <- which(C[x,] == 1)
    new_cl <- sample(setdiff(1:q, c(cl)), 1)
  }))
  old_cols <- max.col(C[rows_to_modify, ])
  C[cbind(rows_to_modify, old_cols)] <- 0
  C[cbind(rows_to_modify, new_cols)] <- 1
  return(C)
}

#' @description runs one simulation for a given set of parameters with an
#' erroneous clustering
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
#' @param error_rate percentage of errors to introduce in the clustering
one_simu_erroneous_clustering <- function(simu, n, p, q, d, kappa,
                                          omega_structure, total_simu = NA,
                                          saving_file = NULL, verbose = 1,
                                          error_rate = 0.15){
  ## Show progression
  if (verbose >= 1) cat("Simu", simu, "of", total_simu)

  ## Draw parameters and data according to the Normal-Block model
  X <- generate_X(n, d)
  true_param <- generate_parameters(X, p, q, kappa, omega_structure)
  data <- NormalBlockData$new(Y = normal_block_data(true_param, X), X = X)

  ## --------------------------------------------------------------
  ## NORMAL BLOCKS WITH FIXED BLOCKS
  ##

  if (verbose > 1) cat("\n V-EM Known fixed blocks")

  args_fixed_blocks <- list(data,
                            sparsity = TRUE,
                            zero_inflation = (min(kappa) > 0),
                            control = NB_control(verbose = FALSE),
                            blocks = true_param$C)

  t0 <- Sys.time()
  myNB_fixed_blocks <- do.call(normal_block, args_fixed_blocks)
  t  <- Sys.time() - t0
  ref_penalty <- myNB_fixed_blocks$get_best_model("BIC")$sparsity


  res_fixed_blocks_BIC   <- get_measures(myNB_fixed_blocks, true_param, "BIC",
                                         fixed_blocks = TRUE)
  res_fixed_blocks_BIC   <- c(res_fixed_blocks_BIC, "time" = as.numeric(t))


  ## --------------------------------------------------------------
  ## NORMAL BLOCKS WITH ERRONEOUS FIXED BLOCKS
  ##
  C_with_errors <- introduce_clustering_errors(true_param$C, error_rate)
  args_fixed_blocks_erroneous_clustering <- list(data,
                                             sparsity = TRUE,
                                             zero_inflation = (min(kappa) > 0),
                                             control = NB_control(verbose = FALSE),
                                             blocks = C_with_errors)
  if (verbose > 1) cat("\n V-EM Known fixed blocks with erroneous clustering")
  t0 <- Sys.time()
  myNB_fixed_blocks_erroneous_clustering <- do.call(normal_block,
                                                args_fixed_blocks_erroneous_clustering)
  t  <- Sys.time() - t0
  ref_penalty <- myNB_fixed_blocks_erroneous_clustering$get_best_model("BIC")$sparsity


  res_fixed_blocks_erroneous_clustering_BIC   <- get_measures(myNB_fixed_blocks_erroneous_clustering,
                                                          true_param, "BIC",
                                                          fixed_blocks = TRUE)
  res_fixed_blocks_erroneous_clustering_BIC   <- c(res_fixed_blocks_erroneous_clustering_BIC,
                                               "time" = as.numeric(t))

  res <- as.data.frame(cbind(simu = simu, n = n, p = p, q = q,
                             zero_inflation = mean(kappa), error_rate = error_rate,
                             omega_structure = omega_structure,
                             rbind(c(method = "NB", res_fixed_blocks_BIC),
                                   c(method = "NB_erroneous_clustering", res_fixed_blocks_erroneous_clustering_BIC))))
  if (verbose > 1) cat("\ndone \r")

  res
}


#' @description runs multiple simulations with an erroneous clustering, for one
#' set of parameters
#' @returns a dataframe with the measures for each simulation
#' @param n_simu number of simulations to run
#' @param n number of rows in the abundance matrix
#' @param p number of columns in the abundance matrix
#' @param q number of clusters in the model
#' @param d number of covariates
#' @param kappa vector of zero-inflation probabilities
#' @param error_rate percentage of errors to introduce in the clustering
#' @param saving_file_prefix  optional, prefix to saving file if the results should
#' be saved
#' @param mc.cores number of cores to run the simulations on
multiple_simulations_erroneous_clustering <-
  function(n_simu, n, p, q, d, kappa, error_rate,
           omega_structure = c("erdos_renyi", "preferential_attachment","community"),
           saving_file_prefix = NULL, mc.cores = max(1, detectCores() - 2)){

    cat("Settings: (n, p, q, d, kappa) = (",n, p, q, d, mean(kappa), ")\n")
    omega_structure <- match.arg(omega_structure)
    multiple_res <- mclapply(1:n_simu,
                             one_simu_erroneous_clustering,
                             n, p, q, d, kappa, omega_structure,
                             error_rate = error_rate,
                             mc.cores = mc.cores, total_simu = n_simu)
    res <- do.call(rbind, multiple_res)
    res
  }


#' @description runs multiple simulations with erreoneous clustering for multiple
#' sets of parameters
#' @returns a dataframe with the measures for each model: integrated inference
#' with fixed blocks, fixed q and baseline inference procedures
#' @param n_simu number of simulations to run for each set of parameters
#' @param n_list list of values for n, number of rows in the abundance matrix
#' @param p_list list of values for p, number of columns in the abundance matrix
#' @param q_list list of values for q number of clusters in the model
#' @param d_list list of values for d, number of covariates
#' @param kappa_list p_list list of values for kappa, vector of zero-inflation probabilities
#' @param error_rate_list list of error rates for the clustering
#' @param saving_folder  optional, folder to save the results in
#' @param mc.cores number of cores to run the simulations on
grid_simulations_erroneous_clustering <- function(n_simu, n_list, p_list, q_list,
                                              d_list, kappa_list, omega_structure_list,
                                              error_rate_list,
                                              saving_folder = NULL,
                                              mc.cores = max(1, detectCores() - 2)){

  settings <- expand.grid(n = n_list,
                          p = p_list,
                          d = d_list,
                          q = q_list,
                          kappa = kappa_list,
                          omega_structure = omega_structure_list,
                          error_rate = error_rate_list,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>% as_tibble()
  settings$n_simu <- n_simu

  final_res <- purrr::pmap(settings, multiple_simulations_erroneous_clustering,
                           mc.cores = mc.cores)
  final_res <- do.call(rbind, final_res) %>% as_tibble()
  final_res
}


############### Run simulations ################################################

n_simu = 2
n_list = c(50)
p_list = c(100)
q_list = c(5)
d_list = c(1)
kappa_list =  list(rep(0, 100))
error_rate_list = c(0.05, 0.1, 0.15)
omega_structure_list = c("erdos_renyi")

res <- grid_simulations_erroneous_clustering(n_simu, n_list, p_list, q_list,
                                         d_list, kappa_list, omega_structure_list,
                                         error_rate_list, saving_folder = NULL,
                                         mc.cores = max(1, detectCores() - 2))

