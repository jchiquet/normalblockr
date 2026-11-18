# Cross-validated log-likelihood of normal-block models

Splits the observations into folds, fits \[normal_block()\] on all but
one fold and evaluates the held-out log-likelihood of the left-out fold
(\`heldout_loglik()\`, see \[NormalBlockBase\]), for every model of the
collection. It gives a criterion to choose the number of clusters that
does not rest on a penalty: the held-out log-likelihood is exact for
both families (for variance-block models, the marginal over the latent W
is computed with the Woodbury identity at the MAP partition), so it
compares models of different families, shapes of Sigma or scalings on
the same footing. Observations, not variables, are split: every fold
sees all the variables, whose clustering is re-estimated on each
training set.

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

  a \[NormalBlockData\] object

- blocks:

  as in \[normal_block()\]: a range of numbers of clusters, a single
  one, or a known clustering (p x q indicator matrix)

- model, zero_inflation, sparsity:

  as in \[normal_block()\]. A sparsity path (\`sparsity = TRUE\`) is not
  supported: its penalty grid depends on the data, hence on the fold;
  pick the penalty with StARS instead (see
  \[NormalBlockVarCollectionSparsity\]).

- folds:

  the number of folds, or a vector of \`data\$n\` fold labels. With a
  number, the folds are drawn at random, stratified on the distinct rows
  of the design X when it is categorical (e.g. a factor coded by
  \`model.matrix()\`), so that every training set sees every category.

- control:

  a list of controls, see \[NB_control()\]

## Value

an object of class \`normal_block_cv\`: a data frame with one row per
model (\`q\`), the total held-out log-likelihood over the \`n\`
observations (\`loglik\`, larger is better) and its standard error
(\`se\`, from the spread of the per-observation values); the
per-observation values are kept in the \`per_observation\` attribute (n
x number of models) and the folds in \`folds\`. Its \`best\` attribute
is the q chosen by the one-standard-error rule (Breiman et al., 1984):
the smallest q whose log-likelihood is within one standard error of the
largest. The held-out log-likelihood typically plateaus beyond the true
number of clusters, since superfluous clusters barely hurt prediction,
so its plain maximum (attribute \`best_max\`) tends to over-select.

## Examples

``` r
ex   <- generate_normal_block_var_data(n = 60, p = 15, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
set.seed(1)
cv <- normal_block_cv(data, blocks = 2:4, folds = 3)
cv
#> Cross-validated log-likelihood, 3 folds (larger is better)
#>  q    loglik       se
#>  2 -1640.221 32.06930
#>  3 -1598.375 29.03992
#>  4 -1582.373 28.29224
#> Best: q = 3 (one-standard-error rule; maximum at q = 4)
```
