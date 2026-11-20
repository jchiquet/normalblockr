# Predict Method for Normal-Block Models

Predicts the mean of the observations Y for new covariates X, on the
original scale of Y: \`B' X\` for the variance-block models (the latent
block signals, centered, average out; unlike \[fitted()\], which adds
their posterior mean for the training observations), \`C B' X\` for the
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
#>          [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 5.099641 8.712870 8.849630 2.240129 9.761127
#> [2,] 2.297330 3.925049 3.986658 1.009153 4.397277
#> [3,] 2.348596 4.012638 4.075622 1.031672 4.495404
```
