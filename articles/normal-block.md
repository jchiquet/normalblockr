# Analyzing multivariate Gaussian data with the Normal-Block model (variance clustering) - First steps

## Preliminaries

This vignette illustrates the use of the
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
function and the methods accompanying its R6 classes.

From a statistical point of view, the
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
function fits a multivariate Normal-Block model (a Gaussian graphical
model with a latent clustering structure) to a table of observations,
possibly after correcting for the effect of covariates. Depending on the
arguments given, the function uses a clustering supplied by the user or
infers one, includes zero-inflation or not, infers a sparse network or
not, and so on. Parameter inference can either be done with an
integrated variational Expectation-Maximization approach (recommended)
or with a faster heuristic approach.

This vignette covers the default family, where the clustering structures
the *covariance* of the variables (`model = "var"`). A complementary
family clusters variables by their *response to covariates* instead
(`model = "mean"`); see
[`vignette("mean-block-breast-cancer")`](https://jchiquet.github.io/normalblockr/articles/mean-block-breast-cancer.md).

Theoretical explanations for this model and full estimation details
(criteria, E/M updates) can be found in Tous and Chiquet (2026).

#### Requirements

``` r

library(normalblockr)
library(pheatmap)
library(paletteer)
```

## Data simulation

We illustrate the analysis with two simulated datasets.

Fix the seed to make the results reproducible.

``` r

set.seed(1)
```

Data are simulated with
[`generate_normal_block_var_data()`](https://jchiquet.github.io/normalblockr/reference/generate_normal_block_var_data.md).

### Fix simulation parameters

Several parameters need to be defined to simulate the data:

``` r

n = 100 # Number of samples (number of rows in the matrix of observations)
p = 40  # Number of entities observed (number of columns in the matrix of observations)
d = 2   # Number of covariates
q = 6   # Number of clusters
kappa = 0 # Mean zero-inflation probability (can also be a vector to define one ZI-probability for each variable). kappa = 0 means that there is no zero-inflation.
omega_structure = "erdos-renyi" # Network structure.
u_v = c(0.3, 0.1) # Parameters to generate an association matrix from a graph, details given in the bibliography.
SNR = 0.75 # Signal to Noise Ratio, defines the relative weight of the covariates and the variance
alpha = rep(1/q, q) # Vector giving probabilities of belonging to each cluster
range_X = c(0, 10)  # Min and max values for the covariates
range_D = c(0.5, 1.5) # Min and max values for the individual entities variances
```

### Simulation 1 (no zero-inflation)

[`generate_normal_block_var_data()`](https://jchiquet.github.io/normalblockr/reference/generate_normal_block_var_data.md)
generates data under the Normal-Block model. It returns a list
containing the simulated covariates $`X`$ and observations $`Y`$, and
the simulation parameters (including the clustering $`C`$).

``` r

my_nb_data <- generate_normal_block_var_data(n, p, d, q, kappa, omega_structure, u_v,
                                            SNR, alpha, range_X, range_D)
```

``` r

pheatmap::pheatmap(my_nb_data$Y,
                   color =paletteer::paletteer_c("ggthemes::Orange-Gold", n = 100), cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE)
```

![](normal-block_files/figure-html/simulation1-visualization-1.png)

### Simulation 2 (with zero-inflation)

To generate zero-inflated data, we set-up kappa as a vector of
zero-inflation probabilities. One needs to ensure that these values are
between 0 and 1.

``` r

kappa_zi <- rnorm(p, mean = 0.7, sd = 0.05)
kappa_zi <- unlist(lapply(kappa_zi, f <- function(x) return(max(0, min(x, 0.9)))))
my_nb_data_zi <- generate_normal_block_var_data(n = n, p = p, d = d, q = 5, kappa = kappa_zi,
                                                omega_structure = omega_structure, u_v = u_v,
                                                SNR = SNR, alpha = rep(1/5, 5),
                                                range_X = range_X, range_D = range_D)
```

For a better visualization of the zero-inflated data, we set-up the 0 to
be shown in white.

``` r

min_val <- min(my_nb_data_zi$Y) ; max_val <- max(my_nb_data_zi$Y)
orange_gold_pal <- paletteer::paletteer_c("ggthemes::Orange-Gold", n = 100)
n_breaks <- 100
zero_pos <- round((0 - min_val) / (max_val - min_val) * n_breaks) + 1
custom_pal <- c(orange_gold_pal[1:(zero_pos - 1)], "white", orange_gold_pal[zero_pos:n_breaks])
pheatmap::pheatmap(my_nb_data_zi$Y,
                   color = custom_pal,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   show_rownames = FALSE)
```

![](normal-block_files/figure-html/simulation2-visualization-1.png)

## Prepare the data for a Normal-Block analysis

A specific data object of class `NormalBlockData` needs to be created to
analyse the data with `normalblockr`.

One can use $`Y`$ and $`X`$ alone to create the `NormalBlockData`
object.

``` r

my_data    <- NormalBlockData$new(my_nb_data$Y, my_nb_data$X)
my_data_zi <- NormalBlockData$new(my_nb_data_zi$Y, my_nb_data_zi$X)
```

Alternatively, if one only wants to include some of the covariates
contained in $`X`$, a formula can be used to specify which one:

``` r

colnames(my_data$X) <- c("X1", "X2")
my_data_alt    <- NormalBlockData$new(my_data$Y, my_data$X, formula = ~ 0 + X1)
```

## Run a Normal-Block analysis

All Normal-Block analyses are run with the
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
function, called with different arguments depending on whether the
clustering (or the number of clusters) is known, and on the requested
level of sparsity, among other factors. See
[`?normal_block`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
for full details.

By default, the parameters are inferred using a variational
Expectation-Maximization approach, with up to 500 iterations. Finer
control of the optimization is possible through the `control` argument
of
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md),
a list generated by
[`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md).

### With non-zero-inflated data

#### Fixed clustering

When the variables’ clustering is known, it can be given directly as an
input to
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md).

``` r

my_NB <- normal_block(data = my_data,
                      blocks = my_nb_data$parameters$C)
#> Fitting a diagonal normal-block-var model with fixed clusters 
#> 
#> DONE
```

Convergence of the model can be checked with
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

``` r

plot(my_NB)
```

![](normal-block_files/figure-html/simulation1-NB1-plot-1.png)

A summary of the results is accessible via
[`print()`](https://rdrr.io/r/base/print.html).

``` r

print(my_NB)
#> A diagonal normal-block-var model with fixed clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL    EBIC niter
#>       141 6      15        0 -3018.334 6036.668 6685.997 6685.997 6739.75    11
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict()
```

[`summary()`](https://rdrr.io/r/base/summary.html) adds the cluster
sizes and the density of the network, and the usual extractors give the
estimated parameters: [`coef()`](https://rdrr.io/r/stats/coef.html) the
regression coefficients B (one column per variable, on the original
scale of $`Y`$), [`sigma()`](https://rdrr.io/r/stats/sigma.html) the
covariance of the latent cluster signals ($`\Omega^{-1}`$, one row and
column per cluster),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) the fitted and
predicted means.

``` r

summary(my_NB)
#> A diagonal normal-block-var model with fixed clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL    EBIC niter
#>       141 6      15        0 -3018.334 6036.668 6685.997 6685.997 6739.75    11
#> ===========================================================================
#> * Cluster sizes: 8, 4, 2, 11, 8, 7
#> * Network: 15 edge(s), density = 1
coef(my_NB)[, 1:6]
#>           [,1]        [,2]      [,3]      [,4]      [,5]      [,6]
#> [1,] 0.5197270 0.133801167 0.2137110 0.1418543 0.2401034 0.2453367
#> [2,] 0.3284179 0.004461164 0.4911172 0.1063739 0.4263735 0.3967010
sigma(my_NB)
#>              Cluster_1   Cluster_2    Cluster_3  Cluster_4    Cluster_5
#> Cluster_1  0.625145335 -0.45789157 -0.330340796  0.4161182  0.004651783
#> Cluster_2 -0.457891568  0.60628347  0.292195595 -0.4287843 -0.045843878
#> Cluster_3 -0.330340796  0.29219560  0.579505187 -0.3188967 -0.003757869
#> Cluster_4  0.416118245 -0.42878429 -0.318896689  0.6200838 -0.153078180
#> Cluster_5  0.004651783 -0.04584388 -0.003757869 -0.1530782  0.482793595
#> Cluster_6 -0.314160969  0.24460537  0.036044440 -0.3342293  0.105034825
#>             Cluster_6
#> Cluster_1 -0.31416097
#> Cluster_2  0.24460537
#> Cluster_3  0.03604444
#> Cluster_4 -0.33422931
#> Cluster_5  0.10503483
#> Cluster_6  0.54858562
```

Within a cluster, the model makes every pair of variables co-vary
positively. `within_block_correlations()` compares the empirical
correlations of the residuals with those the model implies, block by
block, to check this on a given clustering:

``` r

my_NB$within_block_correlations()
#>   block size  mean_cor   min_cor prop_negative model_cor
#> 1     1    8 0.7728622 0.6740550             0 0.7674209
#> 2     2    4 0.7707950 0.7116610             0 0.7572499
#> 3     3    2 0.7695656 0.7695656             0 0.7676822
#> 4     4   11 0.8170062 0.7463831             0 0.8141878
#> 5     5    8 0.6946831 0.6030269             0 0.6857087
#> 6     6    7 0.7171248 0.6187491             0 0.7111590
```

The inter-cluster association network can be visualized with the
`plot_network()` function, or extracted as a matrix with
`latent_network()`.

``` r

 my_NB$plot_network()
```

![](normal-block_files/figure-html/simulation1-NB1-plot_network-1.png)

#### Fixed number of clusters

When the variables’ clustering is unknown, the number of clusters can
simply be fixed via the `blocks` argument.

``` r

my_NB <- normal_block(data = my_data,
                      blocks = q)
#> Fitting a diagonal normal-block-var model with 6 unknown clusters 
#> 
#> DONE
print(my_NB)
#> A diagonal normal-block-var model with 6 unknown clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL     EBIC
#>       146 6      15        0 -3086.581 6173.162 6840.935 6840.935 6894.688
#>  niter
#>     13
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict()
plot(my_NB)
```

![](normal-block_files/figure-html/simulation1-NB2-1.png) The clustering
inferred by
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
can be compared with the true clustering, when it is known:

``` r

aricode::ARI(my_NB$clustering, apply(my_nb_data$parameters$C, 1, which.max))
#> [1] 1
```

#### Unknown number of clusters

When the number of clusters is unknown,
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
can be given a range of candidate values instead, returning a collection
of Normal-Block models, one per number of clusters.

``` r

my_NB_unknown <- normal_block(data = my_data,
                              blocks = 2:8)
#> Fitting a collection of diagonal normal-block-var models over q 
#>   refine(): improved 5 of 7 models
#> 
#> DONE
```

Each model of the collection is first fitted from its own initial
clustering (a $`k`$-means on the residuals, see
[`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)),
then the whole collection is polished by `refine()`: for every $`q`$, it
tries to split a cluster of the model with $`q-1`$ clusters and to merge
two clusters of the model with $`q+1`$, and keeps a candidate when it
lowers the deviance. Since the models are nested, a deviance that
increases with $`q`$ signals a poor local optimum:
[`print()`](https://rdrr.io/r/base/print.html) reports such cases
(`deviance_violations`).

``` r

my_NB_unknown
#> A collection of diagonal normal-block-var models over q 
#> ===========================================================================
#>   7 model(s) explored
#>     q ranging from 2 to 8 
#> ===========================================================================
#> * Useful fields
#>     $models, $criteria
#> * Useful methods
#>     print(), summary(), plot(), $get_best_model()
#> * Warning: the deviance increases with q at q = 7, 8 
#>     (nested models: a sign of poor local optima that refine() did not repair).
```

Here, a violation remains at q = 7, 8, at the edge of the explored
range, where refine() can only seed the model from its smaller
neighbour; it does not affect the selected model.

A specific model can be selected either by its number of clusters or as
the best model for a given criterion (BIC, ICL or EBIC). The value of
each criterion, for every model in the collection, can be visualized
with [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

``` r

plot(my_NB_unknown)
```

![](normal-block_files/figure-html/simulation1-NB3-plot-1.png)

Model selection is done with
[`get_model()`](https://jchiquet.github.io/normalblockr/reference/get_model.md)
when a specific number of clusters is required, and with
`get_best_model()` to select the best model for a given criterion.

``` r

myNB_6   <- my_NB_unknown$get_model(6)
myNB_BIC <- my_NB_unknown$get_best_model("BIC")
myNB_BIC$q
#> [1] 6
```

The criteria above rest on penalties. Cross-validation offers a
penalty-free alternative:
[`normal_block_cv()`](https://jchiquet.github.io/normalblockr/reference/normal_block_cv.md)
fits the collection on all folds but one and computes the exact
log-likelihood of the left-out observations (for the variance-block
model, the marginal over the latent factors, computed at the selected
partition). It typically plateaus beyond the true number of clusters, so
the one-standard-error rule picks the smallest $`q`$ within one standard
error of the best.

``` r

set.seed(1)
my_NB_cv <- normal_block_cv(my_data, blocks = 2:8, folds = 5)
my_NB_cv
#> Cross-validated log-likelihood, 5 folds (larger is better)
#>  q    loglik       se
#>  2 -7411.048 85.71043
#>  3 -7077.635 76.53842
#>  4 -6775.717 67.59015
#>  5 -6553.907 50.90229
#>  6 -6493.623 49.79544
#>  7 -6498.672 49.55779
#>  8 -6496.936 49.39431
#> Best: q = 6 (one-standard-error rule; maximum at q = 6)
plot(my_NB_cv)
```

![](normal-block_files/figure-html/simulation1-NB3-cv-1.png)

Here BIC, ICL and cross-validation all recover the 6 simulated clusters.
By default the columns of $`Y`$ are rescaled for the variance-block
model, which usually helps on real data (see
[`?NormalBlockData`](https://jchiquet.github.io/normalblockr/reference/NormalBlockData.md)).

#### Playing with sparsity

By default, no penalty is applied to the association network. A penalty
can be added via the `sparsity` argument of
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md),
with any of the parametrizations seen above (fixed clustering, fixed
number of clusters, or unknown number of clusters). The examples below
use the known clustering, so that the inferred network can be compared
with the simulated one, which has 9 edges out of 15 possible ones.

To use a fixed $`\ell_1`$ penalty (Tous and Chiquet 2026) on the
network, pass that value to the `sparsity` argument. The larger the
penalty, the sparser the network:

``` r

my_NB_sparse_low <- normal_block(data = my_data,
                                 blocks = my_nb_data$parameters$C,
                                 sparsity = 0.05,
                                 control = NB_control(verbose = FALSE))
my_NB_sparse_high <- normal_block(data = my_data,
                                  blocks = my_nb_data$parameters$C,
                                  sparsity = 0.2,
                                  control = NB_control(verbose = FALSE))
c(low = my_NB_sparse_low$n_edges, high = my_NB_sparse_high$n_edges)
#>  low high 
#>   12    8
```

It is usually hard to know a priori which penalty is best.
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
can instead explore a range of penalties by setting `sparsity = TRUE`,
which returns a collection of Normal-Block models, one per penalty. The
criteria can then be plotted as a function of the penalty, using
[`plot()`](https://rdrr.io/r/graphics/plot.default.html).

``` r

my_NB_sparse <- normal_block(data = my_data,
                             blocks = my_nb_data$parameters$C,
                             sparsity = TRUE,
                             control = NB_control(verbose = FALSE))
plot(my_NB_sparse)
```

![](normal-block_files/figure-html/simulation1-NB-changing-sparsity-1.png)

A penalty is selected as a number of clusters is, with
[`get_model()`](https://jchiquet.github.io/normalblockr/reference/get_model.md)
and `get_best_model()`. The penalized criteria tend to favour weak
penalties, hence dense networks. Stability selection (StARS) refits the
path on subsamples of the observations, the clustering being held fixed,
and keeps the smallest penalty whose network is stable enough; it is
usually better suited to reading the network.

``` r

myNB_sparse_BIC   <- my_NB_sparse$get_best_model("BIC")
set.seed(1)
myNB_sparse_StARS <- my_NB_sparse$get_best_model("StARS")
c(BIC = myNB_sparse_BIC$n_edges, StARS = myNB_sparse_StARS$n_edges)
#>   BIC StARS 
#>    14    15
```

The edges selected by StARS can be compared with those of the simulated
network:

``` r

true_net  <- my_nb_data$parameters$Omega != 0
stars_net <- as.matrix(myNB_sparse_StARS$latent_network("support")) != 0
table(simulated = true_net[upper.tri(true_net)], StARS = stars_net[upper.tri(stars_net)])
#>          StARS
#> simulated TRUE
#>     FALSE    6
#>     TRUE     9
```

``` r

myNB_sparse_StARS$plot_network(layout = igraph::layout_with_fr)
```

![](normal-block_files/figure-html/simulation1-NB-stars-network-1.png)

The penalty and the number of clusters can also vary jointly, by passing
a range of values to `blocks` together with `sparsity = TRUE`.
`get_model(q)` then returns the path of penalties for $`q`$ clusters,
from which a penalty is selected as above.

``` r

my_NB_sparse_unknown <- normal_block(data = my_data,
                                     blocks = 5:7,
                                     sparsity = TRUE,
                                     control = NB_control(verbose = FALSE))
my_NB_sparse_unknown$get_model(6)
#> A collection of diagonal normal-block-var models over sparsity penalties, with fixed q 
#> ===========================================================================
#>   30 model(s) explored
#>     sparsity ranging from 0.00457 to 0.457 
#> ===========================================================================
#> * Useful fields
#>     $models, $criteria
#> * Useful methods
#>     print(), summary(), plot(), $get_best_model()
```

### With zero-inflated data

The process is similar with zero-inflated data, but
`zero_inflation = TRUE` must be passed to
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md).
The example below uses the simulated number of blocks and no penalty on
the network. See
[`vignette("zero-inflated-normal-block")`](https://jchiquet.github.io/normalblockr/articles/zero-inflated-normal-block.md)
for the model and a real-data example.

``` r

my_NB_zi <- normal_block(data = my_data_zi,
                         blocks = 5,
                         zero_inflation = TRUE)
#> Fitting a zero-inflated diagonal normal-block-var model with 5 unknown clusters 
#> 
#> DONE
```

When the data is zero-inflated, the inference process may take longer.

``` r

plot(my_NB_zi)
```

![](normal-block_files/figure-html/simulation2-NB-plot-1.png)

The clustering inference may also be harder.

``` r

aricode::ARI(my_NB_zi$clustering,
             apply(my_nb_data_zi$parameters$C, 1, which.max))
#> [1] 0.5816134
```

## References

Tous, Jeanne, and Julien Chiquet. 2026. “An Integrated Method for
Clustering and Association Network Inference.” *Computational Statistics
& Data Analysis* 219: 108347.
<https://doi.org/10.1016/j.csda.2026.108347>.
