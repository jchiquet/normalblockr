# Predict Method for Variance-Block Models

Predicts observations Y for new covariates X. Specific to the
variance-block family: the mean-block models have their own formula (mu
= C B' X).

## Usage

``` r
# S3 method for class 'NormalBlockBase'
predict(object, new_X, ...)
```

## Arguments

- object:

  An object of class NormalBlockVarBase.

- new_X:

  New set of covariates.

- ...:

  not used, only here for S3 compatibility

## Value

A n\*p prediction matrix for new observations
