# Print a Collection Summary

Print method for objects returned by
\[summary.NormalBlockCollection()\].

## Usage

``` r
# S3 method for class 'summary.NormalBlockCollection'
print(x, ...)
```

## Arguments

- x:

  An object of class \`summary.NormalBlockCollection\`.

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
collection <- normal_block(data, blocks = 2:4, control = NB_control(verbose = FALSE))
print(summary(collection))
#> A collection of diagonal normal-block-var models over q 
#> ===========================================================================
#>     q ranging from 2 to 4 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -710.083 1420.166 1591.379 1591.580 1592.765    14
#>        48 3       3        0 -658.902 1317.804 1503.748 1503.748 1510.340     4
#>        53 4       6        0 -623.475 1246.950 1451.539 1451.539 1468.174     5
```
