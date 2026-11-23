# Print a Collection Summary

Print method for objects returned by
[`summary.NormalBlockCollection()`](https://jchiquet.github.io/normalblockr/reference/summary.NormalBlockCollection.md).

## Usage

``` r
# S3 method for class 'summary.NormalBlockCollection'
print(x, ...)
```

## Arguments

- x:

  An object of class `summary.NormalBlockCollection`.

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns `x`.

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
#>        44 2       1        0 -641.534 1283.067 1454.280 1454.280 1455.666     7
#>        48 3       3        0 -538.911 1077.823 1263.767 1263.767 1270.359    10
#>        53 4       6        0 -515.052 1030.104 1234.693 1235.180 1251.328    52
```
