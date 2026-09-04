# Extract the Covariance Matrix

Extract the covariance matrix \`Omega^-1\`: between latent blocks (q x
q) for the variance-block models, between variables (p x p) for the
mean-block models.

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

The covariance matrix, of size q x q or p x p depending on the model.
