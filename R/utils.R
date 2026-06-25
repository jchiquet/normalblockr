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

# equivalent of check_zero_boundary(check_one_boundary(x)), used to keep
# variational probabilities (tau) away from the {0, 1} boundaries (mirrors
# clip_probabilities() in src/utils_arma.h). Previously spelled out at both
# call sites (NormalBlockUnknownClusters.R/ZINormalBlockUnknownClusters.R).
clip_probabilities <- function(x, zero = .Machine$double.eps) {
  check_zero_boundary(check_one_boundary(x, zero), zero)
}

# Projects a symmetric matrix onto the positive-definite cone by flooring its
# eigenvalues, leaving an already-PD matrix unchanged (up to symmetrization).
# Used by split()/merge() (NormalBlockBase.R): their new_Omegaq is built by
# hand-editing a handful of entries of an existing precision matrix (halving/
# averaging diagonal entries, zero-filling the new row/column), which has no
# general guarantee of staying PD -- and an indefinite Omegaq handed directly
# to the (V)EM solver fails on log_det_sympd(). A plain numeric safeguard
# rather than a principled re-derivation, since the entries being edited are
# themselves already a crude split/merge heuristic, not an exact update.
ensure_pd <- function(M, floor = 1e-6) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  eig$vectors %*% diag(pmax(eig$values, floor), nrow(M)) %*% t(eig$vectors)
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
# fitting several q values (NormalBlockCollectionClusters, NormalBlockCollectionClustersSparsity)
# can get the same residual once, instead of recomputing it once per model
# (see sbm_clustering_path()).
ols_residuals <- function(data) {
  B <- data$XtXm1 %*% data$XtY
  data$Y - data$X %*% B
}

# Iteratively reweighted least squares fit of B under a zero-inflation mask
# (weights = zeros_bar, dm1 re-estimated between iterates), operating on
# `data` alone. Shared by zi_residuals() (collection-level, no model instance
# available yet) and NormalBlockBase$private$zi_diag_normal_inference()
# (per-model, additionally carries along the already-computed kappa) -- same
# math, two different call sites, so the fit itself lives here once.
# The zero-inflation mask makes the weight matrix vary by both row and
# column, so each column of B is solved independently (mirrors
# nb_optim::solve_wls in src/zi_closed_form_solvers.h -- no iterative
# optimizer is needed since the objective is exactly quadratic in B). Uses a
# pseudo-inverse rather than solve() because, e.g. with a one-hot design and
# enough zero inflation, a whole design level can be all-zero for some
# variable, making XtWX exactly singular; ginv() falls back to the
# minimum-norm solution instead of erroring (mirroring arma::solve()'s
# automatic pinv fallback on the C++ side).
# The weighted residual sum of squares (ssq) is floored away from exact zero:
# a variable with very few non-zero observations relative to X's degrees of
# freedom (e.g. a rare species seen at a single station) can be fit exactly
# by B on one iterate, driving its residual -- and so dm1 = nY/ssq -- to
# Inf; on the next iterate, that Inf weight propagates into XtWX, and
# ginv()'s svd() errors on the non-finite input. Flooring ssq keeps dm1 large
# but finite instead.
zi_weighted_fit <- function(data) {
  ssq <- function(B) pmax(colSums(data$zeros_bar * (data$Y - data$X %*% B)^2), .Machine$double.eps)

  B   <- data$XtXm1 %*% data$XtY
  dm1 <- data$nY / ssq(B)
  for (i in 1:3) { # a couple of iterates is enough
    DM1 <- matrix(dm1, data$n, data$p, byrow = TRUE) * data$zeros_bar
    for (j in seq_len(data$p)) {
      w <- DM1[, j]
      XtWX <- crossprod(data$X, data$X * w)
      XtWy <- crossprod(data$X, data$Y[, j] * w)
      B[, j] <- MASS::ginv(XtWX) %*% XtWy
    }
    dm1 <- data$nY / ssq(B)
  }
  list(B = B, dm1 = dm1, R = data$zeros_bar * (data$Y - data$X %*% B))
}

# Zero-inflation analogue of ols_residuals(): the zero-inflation-masked
# residuals used to seed clustering heuristics for ZI models. Zero-inflation
# probabilities (kappa) are deliberately not computed here: the residual only
# depends on the weighted fit of B, not on kappa, so the p logistic
# regressions that would produce it are skipped entirely. Like ols_residuals(),
# factored out so collections fitting several q values can get the same
# residual once instead of recomputing it once per model (see
# sbm_clustering_path()).
zi_residuals <- function(data) zi_weighted_fit(data)$R

# Hierarchical (Ward.D2) clustering tree of the p columns of R by their
# pairwise correlation distance (1 - cor). Shared by
# NormalBlockBase$private$clustering_methods$ward2 (the "ward2"
# clustering_init heuristic, and the fallback whenever any chosen heuristic
# collapses to fewer than q clusters) and sbm_clustering_path()'s own
# fallback below -- same computation, two call sites.
# cor() is NA for any pair involving a (near-)constant column (e.g. a rare
# ZI variable with a single non-zero residual): treated as "uncorrelated"
# (cor = 0, i.e. the neutral, maximal 1 - cor distance) rather than letting
# dist()/hclust() fail on NA input.
ward2_tree <- function(R) {
  cor_R <- suppressWarnings(stats::cor(R))
  cor_R[is.na(cor_R)] <- 0
  stats::hclust(stats::dist(1 - cor_R), method = "ward.D2")
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
  fallback_tree <- ward2_tree(R)

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

# Builds the shared sbm_clustering_path() for a collection over q_list
# (NormalBlockCollectionClusters/NormalBlockCollectionClustersSparsity), or returns NULL
# when the optimization doesn't apply: clustering_init isn't the heuristic
# name "sbm" applied uniformly (it names a different heuristic, or is already
# an explicit clustering/list of clusterings -- in which case every model
# falls back to its own per-q heuristic_clustering() call, as before). ZI
# collections cluster on zi_residuals() (the zero-inflation-masked residual a
# ZI model would otherwise derive itself in zi_diag_normal_inference())
# instead of the plain ols_residuals(). Previously duplicated verbatim in
# both collection classes' initialize().
sbm_path_for_collection <- function(mydata, q_list, zero_inflation, control) {
  if (!identical(control$clustering_init, "sbm")) return(NULL)
  R <- if (zero_inflation) zi_residuals(mydata) else ols_residuals(mydata)
  sbm_clustering_path(R, q_list)
}

