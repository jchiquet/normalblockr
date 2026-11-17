# Summarize a Normal-Block Model

Summarizes a fitted normal-block model: model type, goodness-of-fit
criteria, cluster sizes, and the density of the inferred network between
blocks.

## Usage

``` r
# S3 method for class 'NormalBlockBase'
summary(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

- ...:

  not used, only here for S3 compatibility

## Value

An object of class \`summary.NormalBlockBase\` (a list with the model's
\`who_am_I\`, \`criteria\`, \`cluster_sizes\` and network \`density\`),
printed with a dedicated \[print.summary.NormalBlockBase()\] method.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
summary(model)
#> A diagonal normal-block-var model with 3 unknown blocks .
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        48 3       3        0 -504.149 1008.299 1194.243 1194.243 1200.835    13
#> ===========================================================================
#> * Cluster sizes: 3, 8, 9 
#> * Network: 3 edge(s), density = 1
```
