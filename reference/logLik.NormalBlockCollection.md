# Extract Log-Likelihood of a Collection of Normal-Block Models

Returns the log-likelihood of every model in a collection of
normal-block models (see \[NormalBlockVarCollectionClusters\],
\[NormalBlockVarCollectionSparsity\],
\[NormalBlockVarCollectionClustersSparsity\]).

## Usage

``` r
# S3 method for class 'NormalBlockCollection'
logLik(object, ...)
```

## Arguments

- object:

  An object inheriting from NormalBlockCollection.

- ...:

  not used, only here for S3 compatibility

## Value

A numeric vector of log-likelihood values, one per model in the
collection.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
logLik(models)
#> [1] -566.5896 -423.2197 -403.1926 -403.1973
```
