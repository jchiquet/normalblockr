# Package index

## User-facing functions

Main entry points for fitting a Normal-Block model and simulating data

- [`normal_block()`](https://jchiquet.github.io/normalblockr/reference/normal_block.md)
  : Fit Normal-Block Models
- [`normal_block_sequential()`](https://jchiquet.github.io/normalblockr/reference/normal_block_sequential.md)
  : Cluster Variables in the Mean, then in the Residual Covariance
- [`normal_block_cv()`](https://jchiquet.github.io/normalblockr/reference/normal_block_cv.md)
  : Cross-Validated Log-Likelihood of Normal-Block Models
- [`NB_control()`](https://jchiquet.github.io/normalblockr/reference/NB_control.md)
  : Control Parameters for Normal-Block Models
- [`generate_normal_block_var_data()`](https://jchiquet.github.io/normalblockr/reference/generate_normal_block_var_data.md)
  : Simulate Data from a Normal-Block-Var Model
- [`generate_normal_block_mean_data()`](https://jchiquet.github.io/normalblockr/reference/generate_normal_block_mean_data.md)
  : Simulate Data from a Normal-Block-Mean Model

## Model and data classes

R6 classes returned by normal_block(). NormalBlockData wraps the
responses and design matrix; NormalBlock{Var,Mean}*Clusters (optionally
ZI-prefixed) are single fitted models, with a known or unknown
clustering; NormalBlock{Var,Mean}Collection* are collections of models
explored over a range of q and/or sparsity penalties.

- [`NormalBlockData`](https://jchiquet.github.io/normalblockr/reference/NormalBlockData.md)
  : Data Container for Normal-Block Models
- [`NormalBlockVarKnownClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarKnownClusters.md)
  : Normal-Block-Var Model with Known Clustering
- [`NormalBlockVarUnknownClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarUnknownClusters.md)
  : Normal-Block-Var Model with Unknown Clustering
- [`NormalBlockMeanKnownClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanKnownClusters.md)
  : Normal-Block-Mean Model with Known Clustering
- [`NormalBlockMeanUnknownClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanUnknownClusters.md)
  : Normal-Block-Mean Model with Unknown Clustering
- [`ZINormalBlockVarKnownClusters`](https://jchiquet.github.io/normalblockr/reference/ZINormalBlockVarKnownClusters.md)
  : Zero-Inflated Normal-Block-Var Model with Known Clustering
- [`ZINormalBlockVarUnknownClusters`](https://jchiquet.github.io/normalblockr/reference/ZINormalBlockVarUnknownClusters.md)
  : Zero-Inflated Normal-Block-Var Model with Unknown Clustering
- [`ZINormalBlockMeanKnownClusters`](https://jchiquet.github.io/normalblockr/reference/ZINormalBlockMeanKnownClusters.md)
  : Zero-Inflated Normal-Block-Mean Model with Known Clustering
- [`ZINormalBlockMeanUnknownClusters`](https://jchiquet.github.io/normalblockr/reference/ZINormalBlockMeanUnknownClusters.md)
  : Zero-Inflated Normal-Block-Mean Model with Unknown Clustering
- [`NormalBlockVarCollectionClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionClusters.md)
  : Collection of Normal-Block-Var Models over a Range of Cluster Counts
- [`NormalBlockVarCollectionSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionSparsity.md)
  : Collection of Normal-Block-Var Models over a Sparsity Path
- [`NormalBlockVarCollectionClustersSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarCollectionClustersSparsity.md)
  : Collection of Normal-Block-Var Models over Cluster Counts and
  Sparsity Levels
- [`NormalBlockMeanCollectionClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionClusters.md)
  : Collection of Normal-Block-Mean Models over a Range of Cluster
  Counts
- [`NormalBlockMeanCollectionSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionSparsity.md)
  : Collection of Normal-Block-Mean Models over a Sparsity Path
- [`NormalBlockMeanCollectionClustersSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanCollectionClustersSparsity.md)
  : Collection of Normal-Block-Mean Models over Cluster Counts and
  Sparsity Levels

## Data sets

- [`brca_rppa`](https://jchiquet.github.io/normalblockr/reference/brca_rppa.md)
  : Breast Cancer Proteomics Data (TCGA, RPPA)
- [`onema`](https://jchiquet.github.io/normalblockr/reference/onema.md)
  : French Stream Fish Community Data (ONEMA / OFB Electrofishing
  Surveys)
- [`university`](https://jchiquet.github.io/normalblockr/reference/university.md)
  : University Webpages Text Data (CMU "4 Universities" / WebKB)

## S3 methods

Standard extractors and methods for a fitted model (any NormalBlockBase
subclass) and for a collection of models (any NormalBlockCollection
subclass)

- [`coef(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/coef.NormalBlockBase.md)
  : Extract Model Coefficients
- [`fitted(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/fitted.NormalBlockBase.md)
  : Extract Fitted Values
- [`predict(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/predict.NormalBlockBase.md)
  : Predict Method for Normal-Block Models
- [`sigma(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/sigma.NormalBlockBase.md)
  : Extract the Covariance Matrix
- [`print(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.NormalBlockBase.md)
  : Print a Normal-Block Model
- [`summary(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/summary.NormalBlockBase.md)
  : Summarize a Normal-Block Model
- [`print(`*`<summary.NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.summary.NormalBlockBase.md)
  : Print a Normal-Block Model Summary
- [`plot(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/plot.NormalBlockBase.md)
  : Plot a Normal-Block Model
- [`logLik(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/logLik.NormalBlockBase.md)
  : Extract Log-Likelihood of a Normal-Block Model
- [`BIC(`*`<NormalBlockBase>`*`)`](https://jchiquet.github.io/normalblockr/reference/BIC.NormalBlockBase.md)
  : Bayesian Information Criterion for a Normal-Block Model
- [`plot(`*`<NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/plot.NormalBlockCollection.md)
  : Plot a Collection of Normal-Block Models
- [`print(`*`<NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.NormalBlockCollection.md)
  : Print a Collection of Normal-Block Models
- [`summary(`*`<NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/summary.NormalBlockCollection.md)
  : Summarize a Collection of Normal-Block Models
- [`print(`*`<summary.NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.summary.NormalBlockCollection.md)
  : Print a Collection Summary
- [`logLik(`*`<NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/logLik.NormalBlockCollection.md)
  : Extract Log-Likelihood of a Collection of Normal-Block Models
- [`BIC(`*`<NormalBlockCollection>`*`)`](https://jchiquet.github.io/normalblockr/reference/BIC.NormalBlockCollection.md)
  : Bayesian Information Criterion for a Collection of Normal-Block
  Models
- [`print(`*`<normal_block_sequential>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.normal_block_sequential.md)
  : Print a Sequential Mean-then-Variance Fit
- [`print(`*`<normal_block_cv>`*`)`](https://jchiquet.github.io/normalblockr/reference/print.normal_block_cv.md)
  : Print a Cross-Validated Log-Likelihood
- [`plot(`*`<normal_block_cv>`*`)`](https://jchiquet.github.io/normalblockr/reference/plot.normal_block_cv.md)
  : Plot a Cross-Validated Log-Likelihood

## Internal

Abstract base classes and low-level helpers, documented for developers
and for the methods their subclasses inherit

- [`NormalBlockBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockBase.md)
  : Root Base Class for Normal-Block Models
- [`NormalBlockVarBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockVarBase.md)
  : Base Class for Normal-Block-Var Models
- [`NormalBlockMeanBase`](https://jchiquet.github.io/normalblockr/reference/NormalBlockMeanBase.md)
  : Base Class for Normal-Block-Mean Models
- [`NormalBlockCollection`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollection.md)
  : Base Class for a Collection of Normal-Block Models
- [`NormalBlockCollectionClusters`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClusters.md)
  : Base Class for a Collection of Models over a Range of Cluster Counts
- [`NormalBlockCollectionSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionSparsity.md)
  : Base Class for a Collection of Models over a Sparsity Path
- [`NormalBlockCollectionClustersSparsity`](https://jchiquet.github.io/normalblockr/reference/NormalBlockCollectionClustersSparsity.md)
  : Base Class for a Collection over Cluster Counts and Sparsity Levels
- [`get_model()`](https://jchiquet.github.io/normalblockr/reference/get_model.md)
  : Create a Normal-Block Model Object
- [`isNB()`](https://jchiquet.github.io/normalblockr/reference/isNB.md)
  : Check if an Object is a Normal-Block Model
