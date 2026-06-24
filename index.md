# normalblockr

Normal-Block is a graphical model designed for the multivariate analysis
of continuous data. It clusters variables and builds on Graphical-Lasso
to infer a network of statistical dependencies between clusters.

## Installation

You can install the development version of normalblockr from
[GitHub](https://github.com/) with:

``` r

# install.packages("pak")
pak::pak("jeannetous/normalblockr")
```

## Usage

All fitting is done using the function normal_block. An object of the
class normal_data must be created (using observations Y and covariates
X) and used as an input to normal_block. normal_block parameters allow
to choose between sparse / not sparse models, fixed or unknown blocks…

The package comes with a small artificial data set simulated under the
model.

``` r

testdata <- readRDS("tests/testthat/testdata/testdata_normal.RDS")
Y        <- testdata$Y ; X <- testdata$X
C        <- testdata$parameters$C ; Q <- ncol(C)
data     <- NormalBlockData$new(Y, X) 
model    <- normal_block(data, blocks = C)
```

### Sparse model

``` r

model    <- normal_block(data, blocks = C, sparsity = TRUE)
```

### With unobserved clustering

``` r

model    <- normal_block(data, blocks = 2:5, sparsity = TRUE)
model$plot()
model_BIC <- model$get_best_model("BIC")
```

### With zero-inflated data

``` r

testdata_zi <- readRDS("tests/testthat/testdata/testdata_normal_zi.RDS")
Y        <- testdata_zi$Y ; X <- testdata_zi$X
C        <- testdata_zi$parameters$C ; Q <- ncol(C)
data_zi  <- NormalBlockData$new(Y, X) 
model    <- normal_block(data_zi, blocks = C, zero_inflation = TRUE)
```

### References

Please cite our work using the following reference: Tous, J., & Chiquet,
J. (2025). An integrated method for clustering and association network
inference. arXiv preprint arXiv:2503.22467.
