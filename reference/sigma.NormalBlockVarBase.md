# Extract the Latent-Block Covariance Matrix

Extract the covariance matrix between latent blocks (the inverse of
\`Omegaq\`) from a \[NormalBlockVarBase()\] object.

## Usage

``` r
# S3 method for class 'NormalBlockVarBase'
sigma(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockVarBase.

- ...:

  not used, only here for S3 compatibility

## Value

The q x q covariance matrix between latent blocks.
