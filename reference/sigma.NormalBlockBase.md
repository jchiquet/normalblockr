# Extracts the latent-block covariance matrix from objects returned by \[NormalBlockBase()\] and its variants

Extract the covariance matrix between latent blocks (the inverse of
\`Omegaq\`) from a NormalBlockBase object.

## Usage

``` r
# S3 method for class 'NormalBlockBase'
sigma(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

- ...:

  not used, only here for S3 compatibility

## Value

The q x q covariance matrix between latent blocks.
