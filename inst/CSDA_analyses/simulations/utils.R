library(sbm)
library(ClustOfVar)

#' @description for a list of edges, give corresponding (node1, node2) list
edge_to_node <- function(x, n = max(x)) {
  x <- x - 1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1 + 8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  return(data.frame(node1 = i + 1, node2 = j + 1))
}


#' @description turns a clustering given as a list of labels into a one-hot-encoded
#' matrix
#' @returns a matrix of dimensions (p, q) with p the number of clustered elements
#' and q the number of clusters that describes the clustering in a one-hot-encoded fashion
#' @param clustering clustering given as a list of labels
as_indicator <- function(clustering) {
  q <- max(clustering)
  N <- length(clustering)
  Z <- matrix(0, N, q)
  Z[cbind(seq.int(N), clustering)] <- 1
  Z
}


#' @description turns a clustering probability matrix (gives the probability for
#' each element of belonging to each cluster) into a list of clustering labels
#' giving the cluster each element has the max probability of belonging to
#' @returns a list of cluster labels
#' @param tau probability matrix of dimensions (p, q) with p the number of clustered elements
#' and q the number of clusters, tau_j,k is the probability that element k belongs to cluster k
get_clusters <- function(tau) {
  return(apply(tau, 1, which.max))
}
