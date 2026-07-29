# Initialization and refinement strategies for clustering inference in `normalblockr`

This note documents, for reuse in a manuscript's Methods/Discussion, the strategies `normalblockr` uses to initialize and refine the unknown clustering of variables when fitting a Normal-Block model with `NormalBlockUnknownClusters` (and its zero-inflated counterpart). Source: `R/NormalBlockBase.R`, `R/NormalBlockCollectionClusters.R`, `R/utils.R`; empirical support: `inst/clustering_initialization_benchmark/clustering_initialization_benchmark.qmd`.

## 1. Initialization heuristics

### 1.1 Why initialization matters

`NormalBlockUnknownClusters` fits a Normal-Block model by variational EM (VEM), jointly estimating the regression coefficients $B$, the residual variances $D$, the cluster-level precision matrix $\Omega$, and a variational posterior over the assignment of each of the $p$ variables to one of $q$ clusters. This is a non-convex, multi-block coordinate-ascent problem, structurally analogous to $k$-means or a Gaussian mixture: different starting clusterings can converge to substantially different fixed points of the objective (the evidence lower bound, ELBO). Before any (V)EM iteration can run, the model therefore needs an initial hard clustering of the $p$ variables into $q$ groups, derived from the ordinary-least-squares (or, for zero-inflated data, ZI-weighted) residuals $R \in \mathbb{R}^{n \times p}$.

The quantity that should drive this initial clustering is not the residuals' raw values but the *covariance structure* they induce: for two variables $j, j'$ assigned to the same cluster $k$, the model implies $\mathrm{Cov}(Y_j, Y_{j'}) = \mathrm{Var}(W_k)$, a single value shared by every pair in the cluster, regardless of each variable's own idiosyncratic noise variance $D_j$. A heuristic that clusters on residual proximity in the ordinary Euclidean sense is targeting something related but distinct from this.

### 1.2 The heuristic menu

Four heuristics are registered (`private$clustering_methods` in `NormalBlockBase.R`), selectable via `NB_control(clustering_init = ...)`:

| method | target | fit to the model's assumption | cost |
|---|---|---|---|
| `kmeans` | residual *values*, Euclidean | weak -- not the model's actual clustering criterion | cheap |
| `ward2` | pairwise *correlation* ($1 - \mathrm{cor}(R)$ distance, Ward's minimum-variance linkage) | partial -- confounded by heterogeneous $D_j$, since correlation divides by $\sqrt{\mathrm{Var}(W_k) + D_j}$ | cheap |
| `sbm` | block-constant *covariance*, fit as a Gaussian stochastic block model on $\mathrm{cov}(R)$ | exact match -- a Gaussian SBM literally assumes constant within/between-block covariance | expensive |
| `spectral` | top-$q$ eigenvectors of $\mathrm{cov}(R)$, row-normalized (Ng-Jordan-Weiss), then $k$-means | same target as `sbm`, via a cheap spectral proxy instead of a combinatorial block-model fit | cheap |

`ward2`'s tree is built once as `hclust(dist(1 - cor(R)), method = "ward.D2")` and cut at $q$; it is also the deterministic fallback used whenever a chosen heuristic collapses to fewer than $q$ distinct clusters (guaranteeing exactly $q$ non-empty groups via `cutree()`). For a collection of models spanning a range of $q$ values, `sbm` is not re-run independently for every $q$: a single `sbm::estimateSimpleSBM` exploration is run over the whole range once and reused, with `ward2` filling in any $q$ the exploration does not reach.

### 1.3 No single heuristic dominates

A benchmark across three real datasets (breast-cancer reverse-phase protein arrays, $p=163$; French stream fish biomass, $p=46$; university webpage term frequencies, $p=100$) fitting every $(q, \text{heuristic})$ combination independently and ranking by BIC and ICL found no uniformly best heuristic: `kmeans` has the best raw BIC rank overall but is also the least reliable by a monotonicity diagnostic specific to this model class (see below); `ward2` has a marginally worse raw rank but is by far the most reliable; `sbm` and `spectral` are competitive on some datasets and weak on others. `ward2` was retained as the package default precisely because it gives the best balance of raw fit quality and reliability, not because it is uniformly best.

The reliability diagnostic exploits a structural property of the model: at a fixed sparsity level, deviance ($-2 \times$ log-likelihood/ELBO) should be non-increasing in $q$, since a finer partition into blocks is always at least as expressive as a coarser one. A heuristic whose deviance path *increases* somewhere reveals that the VEM converged to a worse local optimum at that $q$ than it did at a smaller one -- an initialization failure, not merely a mediocre fit. Counting these violations gives a reliability measure that a raw BIC/ICL ranking, evaluated in isolation at each $q$, misses entirely.

### 1.4 Multi-start initialization: `best_of_inits()`

Because different heuristics can converge to substantially different VEM local optima at the *same* $q$ (on real data, two heuristics' converged clusterings at fixed $q$ can differ with an Adjusted Rand Index as low as $\approx 0.5$, and their final ELBOs by several hundred nats), `NormalBlockBase` exposes `best_of_inits()`: a multi-start procedure that tries several heuristics and keeps the best-ELBO converged fit. To keep this affordable, every candidate heuristic is first screened with a short partial VEM run (a handful of iterations), and only the best-screened candidates are fully retrained to convergence -- the same coarse-to-fine idiom used elsewhere in the package for split/merge candidate scoring (Section 2). Empirically, the ranking induced by a short screening run already reproduces the ranking a full run would give, at a fraction of the cost. `sbm` is excluded from the default candidate set for `best_of_inits()`, since its cost is dominated by the block-model fit itself rather than by the VEM iterations the screening step economizes on.

### 1.5 A rejected alternative: local discrete reassignment

Variational inference for network/block models (e.g. stochastic block models, SBM) is sometimes refined post-convergence by a discrete, Kernighan-Lin-style local search: greedily reassigning single nodes between blocks and accepting moves that improve a criterion, to escape local optima of the coordinate-ascent recursion. This was considered for `normalblockr`'s per-variable cluster posterior and rejected on both theoretical and empirical grounds specific to this model's structure.

In a SBM, a node's variational posterior over block membership is informed only by its (typically few) incident edges, so it can remain genuinely fractional/ambiguous at convergence, leaving room for a local reassignment to find an improving move. In the Normal-Block model, by contrast, a variable's posterior aggregates evidence over *all* $n$ samples (its cluster-membership update at the E-step is a closed-form softmax that is, moreover, already the exact block-coordinate maximizer of the ELBO given the current variational factor means) -- and it concentrates correspondingly much faster. Empirically, this posterior is already essentially one-hot, with log-odds margins between the best and second-best cluster in the tens to hundreds, once $n \gtrsim 30$, across a wide range of signal-to-noise ratios, numbers of clusters, and both simulated and real data. There is therefore essentially never an "on the fence" variable for a local swap to find or exploit; a naive swap evaluated with the variational factor means held fixed is in fact guaranteed, by the concavity of the entropy-regularized assignment objective, to never improve on the current (already block-optimal) soft assignment. Genuine local optima do exist in this model (different heuristics do converge to different fixed points, Section 1.4), but they manifest as different global initialization basins rather than as locally ambiguous individual assignments -- which is why a multi-start strategy over initializations, rather than a post-hoc local discrete search, is the effective remedy here.

## 2. The split/merge refinement strategy (`refine()`)

### 2.1 Motivation

`normal_block()` can return a collection of independently fitted models, one per candidate number of clusters $q$ in a user-specified range. Each model in this collection is cold-started from its own heuristic clustering (Section 1) and optimized on its own -- so, even with a good single-heuristic choice, a given $q$ can still settle into a milder local optimum than an *incremental* search, seeded from an adjacent, already-fitted $q$, would find. `NormalBlockCollectionClusters$refine()` implements this incremental, neighbor-seeded search as a post-hoc pass over an already-fitted collection.

### 2.2 Split and merge moves

Two deterministic, warm-starting structural moves are defined on a converged model:

- **`split(index)`**: splits cluster `index` into two. Member variables are divided by the median of their estimated residual variance $1/\hat D_j$ within that cluster; the higher-variance half keeps the original cluster label, the lower-variance half is reassigned to a new cluster $q+1$ (implemented as a soft split of the posterior membership matrix, offset by a negligible $\epsilon$ to avoid a literal zero probability). The new cluster's variational mean/variance are copied from the parent cluster's, and the precision matrix $\Omega$ is re-derived in closed form from the resulting (already split-consistent) variational moments -- so the resulting $(q+1)$-cluster candidate is warm-started directly from the parent's converged state, not from a fresh heuristic clustering.
- **`merge(indices)`**: the reverse move -- sums the two clusters' posterior membership columns and averages their variational means/variances, producing a warm-started $(q-1)$-cluster candidate; $\Omega$ is likewise re-derived from the merged moments.

`candidates_split()` generates one split candidate per non-singleton cluster of a given model; `candidates_merge()` generates merge candidates for every pair of clusters (capped, for large $q$, to the pairs with the largest $|\Omega_{jk}|$, i.e. the most strongly coupled clusters, since merging near-independent blocks is rarely competitive). Every candidate is given a short, cheap partial re-optimization (a handful of VEM iterations) before being scored, mirroring the screening step of `best_of_inits()` (Section 1.4).

### 2.3 The `refine()` pass

Given an already-fitted collection over a set of $q$ values, `refine()` visits every $q$ that has an adjacent value $q-1$ (respectively $q+1$) already present in the collection, and:

1. builds split candidates (respectively merge candidates) from the *already-fitted* neighboring model;
2. picks the most promising candidate by ICL, fully re-optimizes it;
3. replaces the current model at that $q$ only if the candidate's deviance is *strictly* lower than the current model's -- otherwise the original is kept unchanged.

Because a candidate only ever replaces the incumbent when it is strictly better, `refine()` can only improve, never worsen, any model it touches; run unconditionally over the whole collection, it requires no prior diagnostic of *which* $q$ might need it. Empirically, this recovers meaningful gains over independently cold-started collections and can shift which $q$ a criterion such as ICL selects.

### 2.4 A structural consequence: `refine()` propagates across a wide $q$ range

Because a split candidate at $q$ is seeded from $q-1$'s *already-refined* model within the same collection (not from that $q$'s own original heuristic clustering), a single `refine()` call over a wide, contiguous range of $q$ effectively propagates one basin of attraction forward (and, symmetrically via merges, backward) through the whole chain. A practical consequence, confirmed empirically: where two collections built from different initialization strategies happen to agree (or one heuristic's own exploration falls back to another, as `sbm` does to `ward2` beyond the range it explored) at one end of the range, their `refine()`-d collections can end up identical from that point onward, regardless of what each started from at every other $q$. This means that, over a wide contiguous range, a single `refine()` call absorbs a substantial part of what a more expensive multi-start initialization (Section 1.4) would otherwise have fixed -- to the point that, on one of the three benchmark datasets, the cheapest single heuristic (`kmeans`) combined with `refine()` outright beats the multi-start strategy combined with `refine()` on both BIC and ICL, at a fraction of the cost. Multi-start initialization retains a clear, uncontested advantage only when `refine()` is not subsequently applied, or when the $q$ range considered is narrow or non-contiguous (leaving `refine()` little or nothing to propagate).
