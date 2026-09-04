# Generate Normal Block Mean Data

A function to draw data from the normal block model (see details). The
function returns both the generated data and the corresponding model
parameters, in a list.

## Usage

``` r
generate_normal_block_mean_data(
  n = 100,
  p = 40,
  d = 1,
  q = 3,
  kappa = 0,
  omega_structure = "erdos-renyi",
  u_v = c(0.3, 0.1),
  SNR = 5,
  alpha = rep(1/q, q),
  range_X = c(0, 10)
)
```

## Arguments

- n:

  number of individuals. Default to 100.

- p:

  number of variables. Default to 40.

- d:

  number of covariates. Default to 1.

- q:

  number of groups. Default to 3.

- kappa:

  vector (or scalar) of variable-wise probability of zero inflation.
  Default to 0.

- omega_structure:

  the structure of the graph on which the precision matrix between
  variables is built. Can be a symmetric matrix with p rows/columns or a
  character picked in "erdos-renyi", "preferential_attachment",
  "community" in which case a graph is drawn with sensible generation
  parameters. See generate_precision_matrix for details.

- u_v:

  two-size vector of positive numbers v and u controlling the generation
  of the precision matrix Omega: v scales the off-diagonal elements of
  the precision matrix (magnitude of partial correlations), and u is a
  positive number added to the diagonal elements to ensure
  positive-definiteness. The default value is c(0.3, 0.1).

- SNR:

  Signal to noise ratio: magnitude of the regression parameters B will
  be adjusted so that tr(var(XB)) and tr(Sigma) match the desired SNR.

- alpha:

  the q-size vector of group proportion. Default to rep(1/q, q)

- range_X:

  A 2-size vector defining the range of the uniform distribution used to
  draw values in X, the regressor matrix. Default is c(0, 10)

## Value

A named list with the following element - Y a matrix of responses - X a
regressor/design matrix - a list of model parameters, encompassing - B:
matrix of regression coefficients - C: matrix of group membership -
Omega: precision matrix of the variables - Sigma: covariance matrix of
the variables - kappa: vector of ZI inflation probabilities (one per
variable)
