# Predict Method for Normal-Block Models

Predicts the mean of the observations Y for new covariates X: B' X for
the variance-block models, C B' X for the mean-block models.

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
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 11.567455 5.5188298 10.020763 3.5151962 8.9449171
#> [2,]  1.186254 0.5659616  1.027639 0.3604869 0.9173103
#> [3,]  6.513361 3.1075228  5.642455 1.9793240 5.0366717
```
