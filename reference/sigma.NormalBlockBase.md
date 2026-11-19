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

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
sigma(model)
#>               Cluster_1   Cluster_2     Cluster_3
#> Cluster_1  6.229576e-01 -0.34215315 -9.907212e-05
#> Cluster_2 -3.421532e-01  0.35311703 -4.111713e-02
#> Cluster_3 -9.907212e-05 -0.04111713  2.547881e-01
```
