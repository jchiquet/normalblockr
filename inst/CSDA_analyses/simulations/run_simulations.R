# library(normalblockr)

source("simulations.R")

n_simu     <- 2
n_list     <- c(100)
p_list     <- c(50, 100)
Q_list     <- c(3)
d_list     <- c(1)
kappa_list <- list(rep(0, 50))#, rep(0.2, 50))

omega_structure_list <- c("erdos_renyi", "preferential_attachment", "community")

res <- grid_simulations(n_simu, n_list, p_list, Q_list, d_list, kappa_list,
                        omega_structure_list, mc.cores = 1)

