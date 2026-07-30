# Print a Normal-Block Model

Print a short summary of a fitted normal-block model: model type,
goodness-of-fit criteria, and the useful fields/methods to explore it
further.

## Usage

``` r
# S3 method for class 'NormalBlockVarBase'
print(x, ...)
```

## Arguments

- x:

  An object of class NormalBlockVarBase.

- ...:

  not used, only here for S3 compatibility

## Value

Invisibly returns \`x\`.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
print(model)
#> A diagonal normal-block model with 3 unknown blocks .
#> ===========================================================================
#>  nb_param q n_edges sparsity   loglik deviance      BIC      ICL     EBIC niter
#>        48 3       3        0 -583.199 1166.399 1354.176 1118.873 1360.767     8
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
```
