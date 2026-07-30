# Bayesian Information Criterion for a Collection of Normal-Block Models

Returns the (variational) BIC of every model in a collection of
normal-block models.

## Usage

``` r
# S3 method for class 'NormalBlockVarCollection'
BIC(object, ...)
```

## Arguments

- object:

  An object inheriting from NormalBlockVarCollection.

- ...:

  not used, only here for S3 compatibility

## Value

A numeric vector of BIC values, one per model in the collection.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
BIC(models)
#> [1] 1851.834 1581.300 1580.864 1611.773
```
