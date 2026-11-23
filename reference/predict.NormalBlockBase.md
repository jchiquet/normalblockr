# Predict Method for Normal-Block Models

Predicts the mean of the observations Y for new covariates X, on the
original scale of Y: `B' X` for the variance-block models (the latent
block signals, centered, average out; unlike
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), which adds
their posterior mean for the training observations), `C B' X` for the
mean-block models. For zero-inflated models, the mean of the Gaussian
component.

## Usage

``` r
# S3 method for class 'NormalBlockBase'
predict(object, new_X, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

- new_X:

  New set of covariates.

- ...:

  not used, only here for S3 compatibility

## Value

A n\*p prediction matrix for new observations

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
predict(model, ex_data$X[1:3, , drop = FALSE])[, 1:5]
#>          [,1]     [,2]      [,3]     [,4]     [,5]
#> [1,] 5.688566 6.215028 11.805181 5.413255 9.079170
#> [2,] 4.111972 4.492524  8.533358 3.912964 6.562865
#> [3,] 5.453197 5.957875 11.316731 5.189277 8.703512
```
