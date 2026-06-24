# turns a list clustering of q cluster labels for N elements into a matrix of
# dimensions (p, q) with a one-hot encoding of the clustering
#
# @param clustering a list of labels
as_indicator <- function(clustering) {
  q <- max(clustering)
  p <- length(clustering)
  Z <- matrix(0, p, q)
  Z[cbind(seq.int(p), clustering)] <- 1
  Z
}

# removes machine's 0 to elements equal to  1 in x
check_one_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x >= 1 - zero] <- 1 - zero
  x
}

# adds machine's 0 to elements equal to 0 in x
check_zero_boundary <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x < zero]  <- zero
  x
}

# computes xlogx, setting it to 0 if x = 0
xlogx <- function(x) ifelse(x < .Machine$double.eps, 0, x * log(x))

# computes xlogy, setting it to 0 if x = 0
xlogy <- function(x,y) ifelse(x < .Machine$double.eps, 0, x * log(y))

# computes softmax
softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

# gets cluster labels from probability matrix
get_clusters <- function(tau) {
  apply(tau, 1, which.max)
}

# for a list of edges, give corresponding (node1, node2) list.
edge_to_node <- function(x, n = max(x)) {
  x <- x - 1 ## easier for arithmetic to number edges starting from 0
  n.node <- round((1 + sqrt(1 + 8*n)) / 2) ## n.node * (n.node -1) / 2 = n (if integer)
  j.grid <- cumsum(0:n.node)
  j <- findInterval(x, vec = j.grid)
  i <- x - j.grid[j]
  ## Renumber i and j starting from 1 to stick with R convention
  data.frame(node1 = i + 1, node2 = j + 1)
}

sigmoid <- function(x){
  return(1 / (1 + exp(-x)))
}

# OLS residuals of Y on X (used to seed the clustering heuristics). Factored
# out of NormalBlockBase$multivariate_normal_inference() so that collections
# fitting several q values (NormalBlockUnknownQ, NormalBlockUnknownQChangingSparsity)
# can get the same residual once, instead of recomputing it once per model
# (see sbm_clustering_path()).
ols_residuals <- function(data) {
  B <- data$XtXm1 %*% data$XtY
  data$Y - data$X %*% B
}

# Clusters the (n x p) residual matrix R into every q in q_list using a
# SINGLE sbm::estimateSimpleSBM exploration over [min(q_list), max(q_list)],
# instead of one independent exploration per q. A wide SBM exploration
# already visits every intermediate block count on its way to fitting any
# single one (its cost grows steeply with the explored range, but the whole
# path comes "for free" with the most expensive single point) -- one
# exploration plus near-instant setModel() calls reproduces the same
# clusterings as calling estimateSimpleSBM separately for each q, at a
# fraction of the cost (see inst/normal_block_models.qmd, "Implementation
# notes": benchmarked ~15-90x cheaper depending on the q range, same or
# statistically equivalent ICL after the (V)EM polish).
# estimateSimpleSBM's own model selection is adaptive and can stop exploring
# before exploreMax if its ICL stops improving, so not every q in q_list is
# necessarily covered by a single exploration; and even a stored model for q
# blocks can have an empty block, just like the VEM's own argmax clustering
# (see plot_network()/cluster_sizes' tabulate() fix). Re-running a dedicated,
# increasingly expensive SBM exploration for just those q's would silently
# reintroduce the very cost (and, empirically, no quality benefit -- see
# inst/normal_block_models.qmd) this function exists to avoid, so those
# entries fall back to a single shared, cheap ward2 clustering instead.
# Returns a list of membership vectors (one per q, never NULL), named by q.
sbm_clustering_path <- function(R, q_list) {
  options <- list(verbosity = 0, exploreMin = min(q_list), exploreMax = max(q_list),
                  plot = FALSE, nbCores = 1)
  mySBM <- sbm::estimateSimpleSBM(stats::cov(R), "gaussian", estimOptions = options)
  explored <- mySBM$storedModels$nbBlocks
  fallback_tree <- stats::hclust(stats::dist(1 - stats::cor(R)), method = "ward.D2")

  stats::setNames(
    lapply(q_list, function(q) {
      if (q %in% explored) {
        mySBM$setModel(q)
        memberships <- mySBM$memberships
        if (length(unique(memberships)) == q) return(memberships)
      }
      stats::cutree(fallback_tree, q)
    }),
    q_list
  )
}

