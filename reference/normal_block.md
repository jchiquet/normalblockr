# Normal-block model

Fit a normal-block model with a variational or heuristic algorithm

## Usage

``` r
normal_block(
  data,
  blocks,
  sparsity = 0,
  zero_inflation = FALSE,
  control = NB_control()
)
```

## Arguments

- data:

  NormalBlockData object, contains the matrix of responses (Y, n x p)
  and the design matrix (X, n x d), must be created with
  NormalBlockData\$new.

- blocks:

  either a integer (number of blocks), a vector of integer (list of
  possible number of block) or a p \* q matrix (for indicating block
  membership when its known)

- sparsity:

  either TRUE to run the optimization for different sparsity penalty
  values OR float to run model with a single sparsity penalty value

- zero_inflation:

  boolean to indicate if Y is zero-inflated and adjust fitted model as a
  consequence

- control:

  a list-like structure for detailed control on parameters should be
  generated with NB_control().

## Value

an R6 object with one of the model classes (or a collection of model
objects).

## Examples

``` r
## Normal Data
ex_data <- generate_normal_block_data(n=50, p=50, d=1, q=3)
#> Error in generate_normal_block_data(n = 50, p = 50, d = 1, q = 3): could not find function "generate_normal_block_data"
data <- NormalBlockData$new(ex_data$Y, ex_data$X)
#> Error: object 'ex_data' not found
my_normal_block <- normal_block(data, blocks = 1:6)
#> Error in mydata$Y: object of type 'closure' is not subsettable
if (FALSE) { # \dontrun{
my_normal_block$plot(c("deviance", "BIC", "ICL"))
Y_hat <- my_normal_block$get_best_model()$fitted
plot(data$Y, Y_hat, log = "xy"); abline(0,1)
} # }
## Normal Data with Zero Inflation
ex_data_zi <- generate_normal_block_data(n=50, p=50, d=1, q=3, kappa = rep(0.5,50))
#> Error in generate_normal_block_data(n = 50, p = 50, d = 1, q = 3, kappa = rep(0.5,     50)): could not find function "generate_normal_block_data"
zidata <- NormalBlockData$new(ex_data_zi$Y, ex_data_zi$X)
#> Error: object 'ex_data_zi' not found
my_normal_block <- normal_block(zidata, blocks = 1:6, zero_inflation = TRUE)
#> Error: object 'zidata' not found
```
