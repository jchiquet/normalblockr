# Summarize a Collection of Normal-Block Models

Summarizes a collection of normal-block models: model type, the full
criteria table, and the range of q/sparsity explored.

## Usage

``` r
# S3 method for class 'NormalBlockCollection'
summary(object, ...)
```

## Arguments

- object:

  An object inheriting from NormalBlockCollection.

- ...:

  not used, only here for S3 compatibility

## Value

An object of class \`summary.NormalBlockCollection\` (a list with the
collection's \`who_am_I\`, \`criteria\`, \`q_range\` and
\`sparsity_range\`), printed with a dedicated
\[print.summary.NormalBlockCollection()\] method.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
summary(models)
#> A diagonal normal-block-var model with unknown q 
#> ===========================================================================
#>     q ranging from 2 to 5 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -770.996 1541.993 1713.205 1713.208 1714.592     9
#>        48 3       3        0 -562.753 1125.506 1311.451 1311.451 1318.042     6
#>        53 4       6        0 -533.046 1066.092 1270.680 1270.680 1287.315    17
#>        59 5      10        0 -533.047 1066.094 1293.238 1293.238 1325.427    17
```
