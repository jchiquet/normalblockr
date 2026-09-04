# Create a Normal-Block Model Object

Creates the appropriate normal-block model (or collection of models)
depending on the parametrization.

## Usage

``` r
get_model(
  data,
  blocks,
  sparsity = 0,
  zero_inflation = FALSE,
  control = NB_control(),
  model = c("var", "mean")
)
```

## Arguments

- data:

  contains the matrix of responses (Y) and the design matrix (X).

- blocks:

  either an integer (number of blocks), a vector of integer (list of
  possible number of block) or a p \* q matrix (for indicating block
  membership when its known)

- sparsity:

  boolean to say whether the model should have a changing penalty OR
  float to run model with a single penalty value

- zero_inflation:

  boolean to indicate if Y is zero-inflated and adjust fitted model as a
  consequence

- control:

  a list-like structure for detailed control on parameters should be
  generated with NB_control() for collections of sparse models

- model:

  which model family to fit, "var" (the default) or "mean" – see
  \[normal_block()\]
