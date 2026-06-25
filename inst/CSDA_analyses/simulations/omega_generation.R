library(igraph)

#'@description generates an Erdos-Renyi structure
#'@returns an igraph object Erdos-Renyi graph
#'@param q number of nodes in the graph
#'@param p probability of having an edge between two nodes
erdos_renyi_graph <- function(q, p = 0.5){
  as_adjacency_matrix(sample_gnp(q, p))
}

#'@description generates a preferential attachment structure
#'@returns an igraph object with a preferential attachment structure
#'@param q number of nodes in the graph
preferential_attachment_graph <- function(q){
  as_adjacency_matrix(sample_pa(q, directed = FALSE))
}

#'@description generates a community structure
#'@returns an igraph object with a community structure
#'@param q number of nodes in the graph
#'@param prob vector of probabilities describing the probability of belonging to
#'each community, the number of which is defined by the length of prob
#'@param p_in probability of having an edge between 2 nodes of one community
#'@param p_out probability of having an edge between 2 nodes of different communities
community_graph <- function(q, prob = c(1/2,1/4,1/4), p_in = 0.5, p_out = 0.1) {
  pref_mat <- matrix(p_out, length(prob), length(prob))
  diag(pref_mat) <- p_in
  graph_mat <- as_adjacency_matrix(sample_sbm(q,
                                              pref.matrix = pref_mat,
                                              block.sizes = c(rmultinom(1, q, prob)) ))
  graph_mat
}

#'@description generates a precision matrix
#'@returns a precision matrix of dimensions (q, q)
#'@param q dimension of the precision matrix
#'@param omega_structure graph structure for the precision matrix
#'@param v parameter used to compute omega from a graph
#'@param u parameter used to compute omega from a graph
generate_omega <- function(q, omega_structure = c("erdos_renyi",
                                                  "preferential_attachment",
                                                  "community"),
                           v = 0.3, u = 0.1){
  omega_structure <- match.arg(omega_structure)
  cond <- FALSE
  while(!cond){
    if(omega_structure == "erdos_renyi") G <- erdos_renyi_graph(q)
    if(omega_structure == "preferential_attachment") G <- preferential_attachment_graph(q)
    if(omega_structure == "community") G <- community_graph(q)

    # Ensuring that the network is not empty for AUC to make sense
    if(max(G) == 0){
      off_diag_indices <- which(row(matrix(1:q, q, q)) != col(matrix(1:q, q, q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      G[selected_index[["row"]], selected_index[["col"]]] <- 1
      G[selected_index[["col"]], selected_index[["row"]]] <- 1
    }
    omega_tilde <- G * v
    omega <- omega_tilde + diag(abs(min(eigen(omega_tilde)$values)) + u, q, q)
    # Ensuring that the network is not full for AUC to make sense
    if(min(omega) > 0){ # Ensuring that the network has 0s for AUC to make sense
      off_diag_indices <- which(row(matrix(1:q, q, q)) != col(matrix(1:q, q, q)), arr.ind = TRUE)
      selected_index <- off_diag_indices[sample(nrow(off_diag_indices), 1), ]
      omega[selected_index[["row"]], selected_index[["col"]]] <- 0
      omega[selected_index[["col"]], selected_index[["row"]]] <- 0
    }
    cond <- all(eigen(omega)$values > 0)
  }
  as.matrix(omega)
}
