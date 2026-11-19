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
