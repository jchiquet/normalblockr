# Collection of Normal-Block-Mean Models over a Range of Cluster Counts

R6 class for a collection of mean-block models (\[NormalBlockMeanBase\])
with different numbers of clusters (q). As for
\[NormalBlockVarCollectionClusters\], the q-independent part of the
initial clustering heuristic is computed once for the whole collection.

## Super classes

[`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
-\>
[`NormalBlockCollectionClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.md)
-\> `NormalBlockMeanCollectionClusters`

## Methods

### Public methods

- [`NormalBlockMeanCollectionClusters$clone()`](#method-NormalBlockMeanCollectionClusters-clone)

Inherited methods

- [`NormalBlockCollection$summary()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.html#method-summary)
- [`NormalBlockCollectionClusters$get_best_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-get_best_model)
- [`NormalBlockCollectionClusters$get_model()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-get_model)
- [`NormalBlockCollectionClusters$initialize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-initialize)
- [`NormalBlockCollectionClusters$optimize()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-optimize)
- [`NormalBlockCollectionClusters$plot()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-plot)
- [`NormalBlockCollectionClusters$print()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-print)
- [`NormalBlockCollectionClusters$refine()`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.html#method-refine)

------------------------------------------------------------------------

### `NormalBlockMeanCollectionClusters$clone()`

The objects of this class are cloneable with this method.

#### Usage

    NormalBlockMeanCollectionClusters$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
ex <- generate_normal_block_mean_data(n = 50, p = 20, d = 1, q = 3)
data <- NormalBlockData$new(ex$Y, ex$X)
models <- normal_block(data, blocks = 2:5, model = "mean")
#> Fitting a collection of diagonal normal-block-mean models over q 
#>   refine(): improved 3 of 4 models
#> 
#> DONE
models$plot(c("BIC", "ICL"))

models$get_best_model()
#> A diagonal normal-block-mean model with 3 unknown clusters.
#> ===========================================================================
#>  nb_param q n_edges sparsity    loglik deviance      BIC      ICL     EBIC
#>        25 3       0        0 -1984.862 3969.724 4074.679 4077.595 4074.679
#>  niter
#>      1
#> ===========================================================================
#> * Useful fields
#>     $model_par, $posterior_par / $var_par, $clustering 
#>     $loglik, $BIC, $ICL, $objective, $nb_param, $criteria
#> * Useful S3 methods
#>     print(), summary(), plot(), coef(), sigma(), fitted(), predict() 
```
