# Bayesian Information Criterion for a Normal-Block Model

Extracts the (variational) BIC of a fitted normal-block model, computed
as \`deviance + log(n) \* nb_param\` (lower is better).

## Usage

``` r
# S3 method for class 'NormalBlockVarBase'
BIC(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockVarBase.

- ...:

  not used, only here for S3 compatibility

## Value

A scalar: the (variational) BIC.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
BIC(model)
#> [1] 1373.274
```
