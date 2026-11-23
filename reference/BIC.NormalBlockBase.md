# Bayesian Information Criterion for a Normal-Block Model

Extracts the (variational) BIC of a fitted normal-block model, computed
as `deviance + penalty` (lower is better). The penalty charges each
group of parameters the log of the number of observations informing it
(see the `penalty_terms` field of
[NormalBlockBase](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)),
so it differs from [`stats::BIC()`](https://rdrr.io/r/stats/AIC.html)'s
default `log(nobs) * df` computed from
[`logLik()`](https://rdrr.io/r/stats/logLik.html).

## Usage

``` r
# S3 method for class 'NormalBlockBase'
BIC(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

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
#> [1] 1351.51
```
