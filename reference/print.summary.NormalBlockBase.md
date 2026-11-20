# Print a Normal-Block Model Summary

Print method for objects returned by \[summary.NormalBlockBase()\].

## Usage

``` r
# S3 method for class 'summary.NormalBlockBase'
print(x, ...)
```

## Arguments

- x:

  An object of class \`summary.NormalBlockBase\`.

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
print(summary(model))
#> A diagonal normal-block-var model with 3 unknown clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        48 3       3        0 -670.207 1340.413 1526.358 1526.358 1532.949    11
#> ===========================================================================
#> * Cluster sizes: 10, 4, 6
#> * Network: 3 edge(s), density = 1
```
