# R6 class for a generic normal model

R6 class for a generic normal model

R6 class for a generic normal model

## Public fields

- `Y`:

  the matrix of responses (rescaled column-wise if \`scale = TRUE\`)

- `Y_scale`:

  the per-column standard deviation Y was divided by (all 1's if \`scale
  = FALSE\`)

- `X`:

  the matrix of covariates

- `X0`:

  the matrix of zero-inflation covariates, if applicable

- `formula`:

  describes the relationship between Y and X, and X0 if applicable,
  useful if not all of X's or X0's covariates should be used, should be
  formatted ~ X1 + X2... \| Z1 + Z2... with the Normal formula before
  the \| and the ZI formula after the \|

- `n`:

  sample size

- `d`:

  number of covariates

- `d0`:

  number of zero-inflation covariates, if applicable

- `p`:

  number of variables

- `XtXm1`:

  inverse of XtX, useful for inference

- `XtY`:

  useful for inference

- `npY`:

  total number of non zeros in Y

- `nY`:

  total number of non zeros for each column/variable in Y

- `zeros`:

  where are the zero in Y

- `zeros_bar`:

  where are the non-zeros in Y

## Methods

### Public methods

- [`NormalBlockData$new()`](#method-NormalBlockData-new)

- [`NormalBlockData$clone()`](#method-NormalBlockData-clone)

------------------------------------------------------------------------

### Method `new()`

Create a new \[\`NormalBlockData\`\] object.

#### Usage

    NormalBlockData$new(Y, X, X0 = NULL, formula = NULL, scale = TRUE)

#### Arguments

- `Y`:

  the matrix of responses (called Y in the model).

- `X`:

  design matrix (called X in the model).

- `X0`:

  zero-inflation design matrix, if applicable.

- `formula`:

  describes the relationship between Y and X, useful if not all of X's
  covariates should be used.

- `scale`:

  whether to rescale each column of Y by its own standard deviation (no
  centering). Default TRUE – see the class-level documentation for the
  rationale and its limits.

------------------------------------------------------------------------

### Method `clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockData$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
