# Extract the Covariance Matrix

Extract the covariance matrix \`Omega^-1\`. For the variance-block
models, the q x q covariance of the latent block signals, on the scale
of the fit (it has no original-scale counterpart when the variables of a
block were rescaled by different factors). For the mean-block models,
the p x p residual covariance between variables, on the original scale
of Y.

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

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
sigma(model)
#>            Cluster_1  Cluster_2  Cluster_3
#> Cluster_1  0.3440074 -0.2667135 -0.2148215
#> Cluster_2 -0.2667135  0.4472436  0.1823370
#> Cluster_3 -0.2148215  0.1823370  0.2554245
```
