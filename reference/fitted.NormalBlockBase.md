# Extract Fitted Values

Extract the fitted values of a normal-block model, on the original scale
of Y. For the variance-block models, they include the posterior mean of
the latent block signals, \`X B + E\[W\] C'\`; for the mean-block
models, they are \`X B C'\` (with the posterior cluster probabilities
when the clustering is unknown). Zero-inflated models set the observed
zeros to 0.

## Usage

``` r
# S3 method for class 'NormalBlockBase'
fitted(object, ...)
```

## Arguments

- object:

  An object of class NormalBlockBase.

- ...:

  not used, only here for S3 compatibility

## Value

A matrix of fitted values extracted from the object.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
head(fitted(model)[, 1:5])
#>            [,1]      [,2]       [,3]      [,4]         [,5]
#> [1,]  3.3088902  7.490816  2.7400874 5.6863548  2.553358452
#> [2,] -0.7366125  1.346913 -1.7122577 0.5197714 -1.303165497
#> [3,]  1.0538357  1.874721  1.0595307 1.5083323  0.937759551
#> [4,]  0.9230792  1.639516  0.9290168 1.3196457  0.822038138
#> [5,]  0.2044879  1.092832 -0.0609923 0.7245398  0.004264997
#> [6,]  5.4922124 11.554270  4.8695978 8.9175841  4.452461291
```
