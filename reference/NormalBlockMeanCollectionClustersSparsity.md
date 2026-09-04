# Collection of Mean-Block Models over Cluster Counts and Sparsity Levels

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
over both a range of cluster counts (q) and a sparsity path, i.e. one
\[NormalBlockMeanCollectionSparsity\] per q. Mirrors
\[NormalBlockVarCollectionClustersSparsity\], minus its SBM-path
shortcut for the initial clustering (ill-suited to this family, see
\[NormalBlockMeanCollectionClusters\]).

## Super classes

[`NormalBlockCollection`](NormalBlockCollection.md) -\>
[`NormalBlockCollectionClustersSparsity`](NormalBlockCollectionClustersSparsity.md)
-\> `NormalBlockMeanCollectionClustersSparsity`

## Active bindings

- `who_am_I`:

  a method to print what model is being fitted

## Methods

### Public methods

- [`NormalBlockMeanCollectionClustersSparsity$new()`](#method-NormalBlockMeanCollectionClustersSparsity-initialize)

- [`NormalBlockMeanCollectionClustersSparsity$clone()`](#method-NormalBlockMeanCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClustersSparsity$get_best_model()`](NormalBlockCollectionClustersSparsity.html#method-get_best_model)
- [`NormalBlockCollectionClustersSparsity$get_model()`](NormalBlockCollectionClustersSparsity.html#method-get_model)
- [`NormalBlockCollectionClustersSparsity$plot()`](NormalBlockCollectionClustersSparsity.html#method-plot)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionClustersSparsity$new()`

Create a new \[\`NormalBlockMeanCollectionClustersSparsity\`\] object.

#### Usage

    NormalBlockMeanCollectionClustersSparsity$new(
      mydata,
      q_list,
      control = NB_control()
    )

#### Arguments

- `mydata`:

  object of NormalBlockData class, with responses and design matrix

- `q_list`:

  list of q values (number of groups) in the collection

- `control`:

  structured list of parameters to handle sparsity control

#### Returns

A new \[\`NormalBlockMeanCollectionClustersSparsity\`\] object

------------------------------------------------------------------------

### `NormalBlockMeanCollectionClustersSparsity$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanCollectionClustersSparsity$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 60, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:4, sparsity = TRUE, model = "mean",
                       control = NB_control(n_sparsity_penalties = 4))
#> Fitting a collection of normal-block-mean models with different values of q and different penalties 
#>   number of blocks = 2                penalty = 0.553603              penalty = 0.1192702             penalty = 0.02569597                penalty = 0.00553603                number of blocks = 3                penalty = 0.553603              penalty = 0.1192702             penalty = 0.02569597                penalty = 0.00553603                number of blocks = 4                penalty = 0.553603              penalty = 0.1192702             penalty = 0.02569597                penalty = 0.00553603           
#> DONE
models$plot("BIC")

models$get_best_model("BIC")
#> A full normal-block-mean model with 3 unknown blocks .
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL     EBIC
#>        90 3      65    0.119 -1186.023 2372.045 2740.536 2741.148 2883.356
#>  niter
#>      8
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
```
