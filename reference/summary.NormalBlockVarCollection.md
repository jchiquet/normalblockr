# Summarize a Collection of Normal-Block Models

Summarizes a collection of normal-block models: model type, the full
criteria table, and the range of q/sparsity explored.

## Usage

``` r
# S3 method for class 'NormalBlockVarCollection'
summary(object, ...)
```

## Arguments

- object:

  An object inheriting from NormalBlockVarCollection.

- ...:

  not used, only here for S3 compatibility

## Value

An object of class \`summary.NormalBlockVarCollection\` (a list with the
collection's \`who_am_I\`, \`criteria\`, \`q_range\` and
\`sparsity_range\`), printed with a dedicated
\[print.summary.NormalBlockVarCollection()\] method.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
models <- normal_block(data, blocks = 2:5, control = NB_control(verbose = FALSE))
summary(models)
#> A diagonal normal-block model with unknown q 
#> ===========================================================================
#>     q ranging from 2 to 5 
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        44 2       1        0 -664.666 1329.333 1501.462 1335.045 1502.848     9
#>        48 3       3        0 -570.040 1140.081 1327.858 1097.444 1334.450    13
#>        53 4       6        0 -570.042 1140.085 1347.422 1201.283 1364.057    22
#>        59 5      10        0 -570.041 1140.081 1370.891 1316.179 1403.079    34
```
