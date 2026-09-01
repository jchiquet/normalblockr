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
# call sites (NormalBlockVarUnknownClusters.R/ZINormalBlockVarUnknownClusters.R).
clip_probabilities <- function(x, zero = .Machine$double.eps) {
  check_zero_boundary(check_one_boundary(x, zero), zero)
}

# Projects a symmetric matrix onto the PD cone by flooring its eigenvalues.
# Used by split()/merge()'s new_Omega, hand-edited from an existing
# precision matrix with no general guarantee of staying PD.
ensure_pd <- function(M, floor = 1e-6) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  eig$vectors %*% diag(pmax(eig$values, floor), nrow(M)) %*% t(eig$vectors)
}

# Graphical-lasso estimate of a precision matrix from its covariance estimate,
# with a plain inversion as the fallback when the solver can't produce a finite
# answer. Wraps the in-package solver (src/graphical_lasso.h) so that the R
# reference implementations and the C++ (V)EM cores run the very same code --
# tests/testthat/test-cpp-*.R compare the two at 1e-8, which only holds if they
# share this step exactly.
glasso_omega <- function(Sigma, rho) {
  glasso_out <- graphical_lasso_fit(Sigma, rho)
  if (anyNA(glasso_out$wi)) {
    warning("GLasso fails, the penalty is probably too small and the system badly ",
            "conditionned \n reciprocal condition number =", rcond(Sigma),
            "\n We send back the original matrix and its inverse (unpenalized).",
            call. = FALSE)
    return(chol2inv(chol(Sigma)))
  }
  Matrix::symmpart(glasso_out$wi)
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

# OLS residuals of Y on X, used to seed the clustering heuristics; shared
# across a collection's models instead of recomputing per q.
ols_residuals <- function(data) {
  B <- data$XtXm1 %*% data$XtY
  data$Y - data$X %*% B
}

# Iteratively reweighted least squares fit of B under a zero-inflation mask
# (weights = zeros_bar, dm1 re-estimated between iterates). Shared by
# zi_residuals() and NormalBlockVarBase's zi_diag_normal_inference(). Each
# column of B is solved independently (the mask varies by row and column);
# ginv() rather than solve() since a design level can be all-zero for some
# variable, making XtWX singular. ssq is floored away from 0 to keep dm1
# finite when a variable has very few non-zero observations.
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

# Zero-inflation analogue of ols_residuals(): kappa isn't computed here since
# the residual only depends on the weighted fit of B, not on kappa.
zi_residuals <- function(data) zi_weighted_fit(data)$R

# Ward.D2 clustering tree of the p columns of R by pairwise correlation
# distance (1 - cor); shared by the "ward2" heuristic, its fallback for any
# heuristic that collapses to fewer than q clusters, and sbm_clustering_path().
# cor() is NA for a (near-)constant column: treated as uncorrelated (cor = 0)
# rather than letting dist()/hclust() fail on NA input.
ward2_tree <- function(R) {
  cor_R <- suppressWarnings(stats::cor(R))
  cor_R[is.na(cor_R)] <- 0
  stats::hclust(stats::dist(1 - cor_R), method = "ward.D2")
}

# Clusters R into every q in q_list from a SINGLE sbm::estimateSimpleSBM
# exploration over [min(q_list), max(q_list)] rather than one per q (see
# inst/normal_block_models.qmd, "Implementation notes", for the cost/quality
# comparison). The exploration is adaptive and may stop before exploreMax, or
# return an empty block for some q; either falls back to a shared ward2
# clustering for that q. Returns a list of membership vectors, named by q.
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

# Builds the shared sbm_clustering_path() for a collection over q_list, or
# NULL if clustering_init isn't "sbm" applied uniformly (every model then
# falls back to its own per-q heuristic_clustering()).
sbm_path_for_collection <- function(mydata, q_list, zero_inflation, control) {
  if (!identical(control$clustering_init, "sbm")) return(NULL)
  R <- if (zero_inflation) zi_residuals(mydata) else ols_residuals(mydata)
  sbm_clustering_path(R, q_list)
}

# Whether NB_control(clustering_init = "best_of_inits") was requested (see
# best_of_inits() in NormalBlockVarBase.R). identical() keeps this safe when
# clustering_init is a list or an explicit clustering.
uses_best_of_inits <- function(control) identical(control$clustering_init, "best_of_inits")

