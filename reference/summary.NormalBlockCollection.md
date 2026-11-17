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
#>        44 2       1        0 -607.284 1214.569 1385.782 1385.782 1387.168     8
#>        48 3       3        0 -521.591 1043.181 1229.126 1229.126 1235.717    10
#>        53 4       6        0 -489.346  978.692 1183.280 1183.281 1199.916    14
#>        59 5      10        0 -489.349  978.698 1205.842 1205.843 1238.031    36
```
