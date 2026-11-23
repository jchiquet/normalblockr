# Cross-Validated Log-Likelihood of Normal-Block Models

Splits the observations into folds, fits
[`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
on all but one fold and evaluates the held-out log-likelihood of the
left-out fold (see `heldout_loglik()` in
[NormalBlockBase](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)),
for every model of the collection. It selects the number of clusters
without a penalty, and compares models of different families, shapes of
Sigma or scalings.

## Usage

``` r
normal_block_cv(
  data,
  blocks,
  model = c("var", "mean"),
  zero_inflation = FALSE,
  sparsity = 0,
  folds = 5,
  control = NB_control(verbose = FALSE)
)
```

## Arguments

- data:

  a
  [NormalBlockData](https://jchiquet.github.io/normalblockr/reference/NormalBlockData.md)
  object

- blocks:

  as in
  [`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md):
  a range of numbers of clusters, a single one, or a known clustering (p
  x q indicator matrix)

- model, zero_inflation, sparsity:

  as in
  [`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md).
  A sparsity path (`sparsity = TRUE`) is not supported: its penalty grid
  depends on the data, hence on the fold; pick the penalty with StARS
  instead (see
  [NormalBlockVarCollectionSparsity](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionSparsity.md)).

- folds:

  the number of folds, or a vector of `data$n` fold labels. With a
  number, the folds are drawn at random, stratified on the distinct rows
  of the design X when it is categorical (e.g. a factor coded by
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)), so
  that every training set sees every category.

- control:

  a list of controls, see
  [`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)

## Value

an object of class `normal_block_cv`: a data frame with one row per
model (`q`), the total held-out log-likelihood over the `n` observations
(`loglik`, larger is better) and its standard error (`se`, from the
spread of the per-observation values); the per-observation values are
kept in the `per_observation` attribute (n x number of models) and the
folds in `folds`. Its `best` attribute is the q chosen by the
one-standard-error rule, the smallest q within one standard error of the
maximum; `best_max` is the maximum itself, which tends to over-select.

## Examples

``` r
ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
set.seed(1)
cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
cv
#> Cross-validated log-likelihood, 3 folds (larger is better)
#>  q    loglik       se
#>  2 -1597.146 31.19902
#>  3 -1582.757 30.23785
#>  4 -1595.184 32.54808
#> Best: q = 2 (one-standard-error rule; maximum at q = 3)
```
