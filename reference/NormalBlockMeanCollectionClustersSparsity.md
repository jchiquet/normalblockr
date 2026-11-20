# Collection of Normal-Block-Mean Models over Cluster Counts and Sparsity Levels

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
over both a range of cluster counts (q) and a sparsity path, i.e. one
\[NormalBlockMeanCollectionSparsity\] per q. Mirrors
\[NormalBlockVarCollectionClustersSparsity\].

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionClustersSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.md)
-\> `NormalBlockMeanCollectionClustersSparsity`

## Methods

### Public methods

- [`NormalBlockMeanCollectionClustersSparsity$clone()`](#method-NormalBlockMeanCollectionClustersSparsity-clone)

Inherited methods

- [`NormalBlockCollection$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-optimize)
- [`NormalBlockCollection$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-print)
- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClustersSparsity$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-get_best_model)
- [`NormalBlockCollectionClustersSparsity$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-get_model)
- [`NormalBlockCollectionClustersSparsity$initialize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-initialize)
- [`NormalBlockCollectionClustersSparsity$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.html#method-plot)

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
#> Fitting a collection of full normal-block-mean models over q and sparsity penalties 
#>   number of blocks = 2           
     penalty = 3.916143           
     penalty = 0.8437074           
     penalty = 0.1817713           
     penalty = 0.03916143           
     number of blocks = 3           
     penalty = 3.916143           
     penalty = 0.8437074           
     penalty = 0.1817713           
     penalty = 0.03916143           
     number of blocks = 4           
     penalty = 3.916143           
     penalty = 0.8437074           
     penalty = 0.1817713           
     penalty = 0.03916143           

#> DONE
models$plot("BIC")

models$get_best_model("BIC")
#> A full normal-block-mean model with 3 unknown blocks.
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL     EBIC
#>        86 3      61    0.844 -2236.571 4473.143 4832.046 4833.154 5197.526
#>  niter
#>      7
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
```
