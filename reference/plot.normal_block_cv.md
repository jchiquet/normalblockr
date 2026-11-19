# Plot a Cross-Validated Log-Likelihood

The held-out log-likelihood against the number of clusters, with one
standard error on either side.

## Usage

``` r
# S3 method for class 'normal_block_cv'
plot(x, ...)
```

## Arguments

- x:

  an object of class \`normal_block_cv\`, from \[normal_block_cv()\]

- ...:

  not used, only here for S3 compatibility

## Value

a \[\`ggplot2::ggplot\`\] graph

## Examples

``` r
ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
set.seed(1)
cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
plot(cv)
```
