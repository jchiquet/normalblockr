# Check if an Object is a Normal-Block Model

Checks if a model is of class \[NormalBlockBase()\], i.e. any
normal-block model (variance-block or mean-block family).

## Usage

``` r
isNB(object)
```

## Arguments

- object:

  An R object.

## Value

A boolean telling whether object inherits from the NormalBlockBase
class.

## Examples

``` r
ex_data <- generate_normal_block_var_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
model <- normal_block(data, blocks = 3, control = NB_control(verbose = FALSE))
isNB(model)
#> [1] TRUE
```
