# Summarize a Normal-Block Model

Summarizes a fitted normal-block model: model type, goodness-of-fit
criteria, cluster sizes, and the density of the inferred network between
blocks.

## Usage

``` r
# S3 method for class 'NormalBlockVarBase'
summary(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockVarBase.

- ...:

  not used, only here for S3 compatibility

## Value

An object of class \`summary.NormalBlockVarBase\` (a list with the
model's \`who_am_I\`, \`criteria\`, \`cluster_sizes\` and network
\`density\`), printed with a dedicated
\[print.summary.NormalBlockVarBase()\] method.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
summary(model)
#> A diagonal normal-block model with 3 unknown blocks .
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        48 3       3        0 -580.406 1160.812 1348.589 1125.142 1355.181     8
#> ===========================================================================
#> * Cluster sizes: 9, 6, 5 
#> * Network: 3 edge(s), density = 1
```
