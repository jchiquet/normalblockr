# Print a Cross-Validated Log-Likelihood

Print a Cross-Validated Log-Likelihood

## Usage

``` r
# S3 method for class 'normal_block_cv'
print(x, ...)
```

## Arguments

- x:

  an object of class \`normal_block_cv\`, from \[normal_block_cv()\]

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
set.seed(1)
cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
print(cv)
#> Cross-validated log-likelihood, 3 folds (larger is better)
#>  q    loglik       se
#>  2 -1682.441 27.55630
#>  3 -1614.087 24.04311
#>  4 -1588.535 23.02608
#> Best: q = 4 (one-standard-error rule; maximum at q = 4)
```
