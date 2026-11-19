# Extract Model Coefficients

Extract the regression coefficients B of a normal-block model. For the
variance-block models, one column per variable, on the original scale of
Y even when the variables were rescaled for the fit. For the mean-block
models, one column per cluster, on the scale of the fit (the original
one by default), since a profile is shared by variables of different
scales.

## Usage

``` r
# S3 method for class 'NormalBlockBase'
coef(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

- ...:

  not used, only here for S3 compatibility

## Value

A d x p (variance-block) or d x q (mean-block) matrix.
