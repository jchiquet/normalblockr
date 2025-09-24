#' @description generates a diagonal variance matrix
#' @returns diagonal variance matrix D
#' @param p dimension required for the matrix
#' @param min_D minimum value for D values
#' @param max_D maximum value for D values
generate_D <- function(p, min_D = 0.5, max_D = 1.5){
  D <- matrix(rep(0, p*p), nrow = p)
  diag(D) <- runif(p, min_D, max_D)
  return(D)
}

#' @description generates a covariate matrix
#' @returns an (n, d) covariates matrix X
#' @param n number of rows for X
#' @param d number of columns for X
#' @param min_X minimum value for X (either a single value for the entire matrix
#' or a list of length d to give one minimal value per dimension)
#' @param max_X maximum value for X (either a single value for the entire matrix
#' or a list of length d to give one minimal value per dimension)
generate_X <- function(n, d, min_X = 0, max_X = 10){
  X = matrix(rep(1, n * d), nrow=n)
  for(dim in 1:d){X[,dim] = runif(n, min=min_X[[dim]], max = max_X[[dim]])}
  return(X)
}

#' @description generates a regression matrix
#' @returns a (d, p) regression matrix B
#' @param p number of columns for B
#' @param X covariates matrix of dimensions (n, d)
#' @param Sigma variance-covariance matrix used in the simulation
#' @param SNR signal-to-noise ratio to define the relative weights of Sigma and XB
generate_B <- function(p, X, Sigma, SNR = 0.75){
  d <- ncol(X)
  B <- matrix(rep(1, d*p), nrow=d)
  for(dim in 1:d){B[dim,] = runif(p, min=0, max = 1)}
  correcting_factor <- SNR * var(as.vector(Sigma)) / (var(as.vector(X %*% B)))
  B <- sqrt(correcting_factor) * B
  return(B)
}

#' @description generates a clustering matrix
#' @returns a (p, q) clustering matrix C (one-hot encoding of the clustering)
#' @param p number of rows for C
#' @param q number of clusters = number of columns for C
#' @param alpha list of probabilities for each cluster, default is equiprobable
#' distribution
generate_blocks <- function(p, q, alpha = NULL){
  if(is.null(alpha)) alpha <- rep(1/q, q)
  C <- matrix(rep(0, p * q), nrow = p)
  while(0 %in% colSums(C)){
    groups = sample(1 : q, size = p, replace = TRUE)
    for(dim in 1:p){C[dim, groups[[dim]]] = 1}
  }
  return(C)
}

#' @description generates a list of parameters required to generate data under
#' the normal-block model
#' @returns a named list with all parameters
#' @param X covaroates matrix
#' @param p number of columns in the final abundance matrix
#' @param q number of clusters
#' @param kappa vector of zero-inflation probabilities
#' @param omega_structure network type of structure for omega
#' @param SNR signal-to-noise ratio, defines the relative weight of Sigma and XB
#' in the simulation
generate_parameters <- function(X, p, q, kappa,
                                omega_structure = c("erdos_renyi", "community",
                                                    "preferential_attachment"),
                                SNR = 0.75) {
  omega_structure <- match.arg(omega_structure)
  Omega <- generate_omega(q, omega_structure)
  Sigma <- chol2inv(chol(Omega))
  list(
    B = generate_B(p, X, Sigma, SNR),
    C = generate_blocks(p, q),
    D = generate_D(p),
    Omega = Omega,
    Sigma = Sigma,
    kappa = kappa
  )
}
