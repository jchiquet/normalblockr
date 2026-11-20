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
#> A collection of diagonal normal-block-var models over q 
#> ===========================================================================
#>     q ranging from 2 to 5 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -733.207 1466.415 1637.627 1637.627 1639.014    10
#>        48 3       3        0 -628.158 1256.317 1442.261 1442.261 1448.853     7
#>        53 4       6        0 -596.138 1192.277 1396.865 1396.865 1413.501     6
#>        59 5      10        0 -596.241 1192.481 1419.626 1419.626 1451.814    10
```
