# Print a Sequential Mean-then-Variance Fit

Reports both stages and, when aricode is available, how related the two
partitions are.

## Usage

``` r
# S3 method for class 'normal_block_sequential'
print(x, ...)
```

## Arguments

- x:

  an object of class \`normal_block_sequential\`

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex   <- generate_normal_block_mean_data(n = 80, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
fit  <- normal_block_sequential(data, blocks_mean = 3, blocks_var = 2)
print(fit)
#> A sequential mean-then-variance normal-block fit
#> ===========================================================================
#>   mean-block stage    : 3 clusters -- diagonal normal-block-mean model with 3
#>                         unknown clusters
#>   variance-block stage: 2 clusters -- diagonal normal-block-var model with 2
#>                         unknown clusters
#>   ARI between the two partitions: -0.015 
#>   (near 0 means the two structures are unrelated, as is usual)
#> ===========================================================================
#> * Useful fields
#>     $mean, $var (both ordinary fitted models), $residuals 
```
