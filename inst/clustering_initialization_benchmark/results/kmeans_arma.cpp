// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Lloyd's k-means with multiple random restarts (mirrors stats::kmeans(X, k,
// nstart = ...) -- same multi-start strategy, but plain Lloyd's reassignment
// instead of R's default Hartigan-Wong algorithm). Used to test whether a
// hand-rolled Armadillo implementation is worth adding to the package's
// clustering_methods registry (R/NormalBlockBase.R) in place of/alongside
// stats::kmeans().

namespace {

double assign_clusters(const arma::mat& X, const arma::mat& centers, arma::uvec& assignment) {
  arma::uword n = X.n_rows, k = centers.n_rows;
  double wcss = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    double best_d = arma::datum::inf;
    arma::uword best_k = 0;
    for (arma::uword c = 0; c < k; ++c) {
      double d = arma::accu(arma::square(X.row(i) - centers.row(c)));
      if (d < best_d) { best_d = d; best_k = c; }
    }
    assignment(i) = best_k;
    wcss += best_d;
  }
  return wcss;
}

arma::uvec kmeans_one_run(const arma::mat& X, arma::uword k, arma::uword max_iter,
                          const arma::uvec& init_idx, double& out_wcss) {
  arma::uword n = X.n_rows, p = X.n_cols;
  arma::mat centers = X.rows(init_idx);
  arma::uvec assignment(n, arma::fill::zeros);
  arma::uvec prev_assignment(n);
  prev_assignment.fill(n); // sentinel: no real cluster id equals n
  double wcss = 0.0;

  for (arma::uword iter = 0; iter < max_iter; ++iter) {
    wcss = assign_clusters(X, centers, assignment);
    if (arma::all(assignment == prev_assignment)) break;
    prev_assignment = assignment;

    arma::mat new_centers(k, p, arma::fill::zeros);
    arma::uvec counts(k, arma::fill::zeros);
    for (arma::uword i = 0; i < n; ++i) {
      new_centers.row(assignment(i)) += X.row(i);
      counts(assignment(i))++;
    }
    for (arma::uword c = 0; c < k; ++c) {
      if (counts(c) > 0) {
        new_centers.row(c) /= static_cast<double>(counts(c));
      } else {
        // re-seed an empty cluster on a random point, mirroring how most
        // k-means implementations handle this degenerate case
        arma::uword ridx = arma::as_scalar(arma::randi<arma::uvec>(1, arma::distr_param(0, n - 1)));
        new_centers.row(c) = X.row(ridx);
      }
    }
    centers = new_centers;
  }
  out_wcss = wcss;
  return assignment;
}

// k-means++ seeding (Arthur & Vassilvitskii, 2007): the first centroid is
// picked uniformly at random, each subsequent one is picked among the
// remaining points with probability proportional to its squared distance to
// the nearest centroid already chosen. Spreads the initial centroids out,
// which is the whole point of testing it here: does it let nstart be much
// smaller while still matching stats::kmeans's quality?
arma::uvec kmeanspp_init(const arma::mat& X, arma::uword k) {
  arma::uword n = X.n_rows;
  arma::uvec chosen(k);
  arma::vec min_dist2(n);
  min_dist2.fill(arma::datum::inf);

  chosen(0) = arma::as_scalar(arma::randi<arma::uvec>(1, arma::distr_param(0, n - 1)));

  for (arma::uword c = 1; c < k; ++c) {
    arma::uword last = chosen(c - 1);
    for (arma::uword i = 0; i < n; ++i) {
      double d2 = arma::accu(arma::square(X.row(i) - X.row(last)));
      if (d2 < min_dist2(i)) min_dist2(i) = d2;
    }
    double total = arma::accu(min_dist2);
    if (total <= 0.0) {
      chosen(c) = arma::as_scalar(arma::randi<arma::uvec>(1, arma::distr_param(0, n - 1)));
      continue;
    }
    double r = arma::randu() * total;
    double cum = 0.0;
    arma::uword sel = n - 1;
    for (arma::uword i = 0; i < n; ++i) {
      cum += min_dist2(i);
      if (cum >= r) { sel = i; break; }
    }
    chosen(c) = sel;
  }
  return chosen;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List kmeans_arma(const arma::mat& X, int k, int nstart = 30, int max_iter = 100,
                       bool kmeanspp = false) {
  arma::uword n = X.n_rows;
  double best_wcss = arma::datum::inf;
  arma::uvec best_assignment;

  for (int s = 0; s < nstart; ++s) {
    arma::uvec init_idx = kmeanspp ? kmeanspp_init(X, static_cast<arma::uword>(k))
                                   : arma::randperm(n, static_cast<arma::uword>(k));
    double wcss;
    arma::uvec assignment = kmeans_one_run(X, static_cast<arma::uword>(k),
                                           static_cast<arma::uword>(max_iter), init_idx, wcss);
    if (wcss < best_wcss) {
      best_wcss = wcss;
      best_assignment = assignment;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("cluster") = best_assignment + 1, // 1-indexed for R
    Rcpp::Named("wcss") = best_wcss
  );
}
