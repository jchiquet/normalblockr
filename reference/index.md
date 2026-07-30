# Package index

## User-facing functions

Main entry points for fitting a Normal-Block model and simulating data

- [`normal_block()`](normal_block.md) : Normal-block model
- [`NB_control()`](NB_control.md) : NB_control
- [`generate_normal_block_var_data()`](generate_normal_block_var_data.md)
  : Generate Normal Block Data

## Model and data classes

R6 classes returned by normal_block(). NormalBlockData wraps the
responses and design matrix; NormalBlockVar*Clusters (optionally
ZI-prefixed) are single fitted models, known- or unknown-clustering;
NormalBlockVarCollection* are collections of models explored over a
range of q and/or sparsity penalties.

- [`NormalBlockData`](NormalBlockData.md) : Data Container for
  Normal-Block Models
- [`NormalBlockVarKnownClusters`](NormalBlockVarKnownClusters.md) :
  Normal-Block Model with Known Clustering
- [`NormalBlockVarUnknownClusters`](NormalBlockVarUnknownClusters.md) :
  Normal-Block Model with Unknown Clustering
- [`ZINormalBlockVarKnownClusters`](ZINormalBlockVarKnownClusters.md) :
  Zero-Inflated Normal-Block Model with Known Clustering
- [`ZINormalBlockVarUnknownClusters`](ZINormalBlockVarUnknownClusters.md)
  : Zero-Inflated Normal-Block Model with Unknown Clustering
- [`NormalBlockVarCollectionClusters`](NormalBlockVarCollectionClusters.md)
  : Collection of Normal-Block Models over a Range of Cluster Counts
- [`NormalBlockVarCollectionSparsity`](NormalBlockVarCollectionSparsity.md)
  : Collection of Normal-Block Models over a Sparsity Path
- [`NormalBlockVarCollectionClustersSparsity`](NormalBlockVarCollectionClustersSparsity.md)
  : Collection of Normal-Block Models over Cluster Counts and Sparsity
  Levels

## Data sets

- [`brca_rppa`](brca_rppa.md) : Breast cancer proteomics data (TCGA,
  RPPA)
- [`onema`](onema.md) : French stream fish community data (ONEMA / OFB
  electrofishing surveys)
- [`university`](university.md) : University webpages text data (CMU "4
  Universities" / WebKB)

## S3 methods

Standard extractors and methods for a fitted model (any
NormalBlockVarBase subclass) and for a collection of models (any
NormalBlockVarCollection subclass)

- [`coef(`*`<NormalBlockVarBase>`*`)`](coef.NormalBlockVarBase.md) :
  Extract Model Coefficients
- [`fitted(`*`<NormalBlockVarBase>`*`)`](fitted.NormalBlockVarBase.md) :
  Extract Fitted Values
- [`predict(`*`<NormalBlockVarBase>`*`)`](predict.NormalBlockVarBase.md)
  : Predict Method for Normal-Block Models
- [`sigma(`*`<NormalBlockVarBase>`*`)`](sigma.NormalBlockVarBase.md) :
  Extract the Latent-Block Covariance Matrix
- [`print(`*`<NormalBlockVarBase>`*`)`](print.NormalBlockVarBase.md) :
  Print a Normal-Block Model
- [`summary(`*`<NormalBlockVarBase>`*`)`](summary.NormalBlockVarBase.md)
  : Summarize a Normal-Block Model
- [`print(`*`<summary.NormalBlockVarBase>`*`)`](print.summary.NormalBlockVarBase.md)
  : Print a Normal-Block Model Summary
- [`plot(`*`<NormalBlockVarBase>`*`)`](plot.NormalBlockVarBase.md) :
  Plot a Normal-Block Model
- [`logLik(`*`<NormalBlockVarBase>`*`)`](logLik.NormalBlockVarBase.md) :
  Extract Log-Likelihood of a Normal-Block Model
- [`BIC(`*`<NormalBlockVarBase>`*`)`](BIC.NormalBlockVarBase.md) :
  Bayesian Information Criterion for a Normal-Block Model
- [`print(`*`<NormalBlockVarCollection>`*`)`](print.NormalBlockVarCollection.md)
  : Print a Collection of Normal-Block Models
- [`summary(`*`<NormalBlockVarCollection>`*`)`](summary.NormalBlockVarCollection.md)
  : Summarize a Collection of Normal-Block Models
- [`print(`*`<summary.NormalBlockVarCollection>`*`)`](print.summary.NormalBlockVarCollection.md)
  : Print a Collection Summary
- [`logLik(`*`<NormalBlockVarCollection>`*`)`](logLik.NormalBlockVarCollection.md)
  : Extract Log-Likelihood of a Collection of Normal-Block Models
- [`BIC(`*`<NormalBlockVarCollection>`*`)`](BIC.NormalBlockVarCollection.md)
  : Bayesian Information Criterion for a Collection of Normal-Block
  Models

## Internal

Abstract base classes, low-level helpers and exploratory code kept for
reference; not part of the user-facing API

- [`NormalBlockVarBase`](NormalBlockVarBase.md) : Base Class for
  Normal-Block Models
- [`NormalBlockVarCollection`](NormalBlockVarCollection.md) : Base Class
  for a Collection of Normal-Block Models
- [`SelectionNClusters`](SelectionNClusters.md) : Select the Number of
  Clusters by Split/Merge Search
- [`get_model()`](get_model.md) : Create a Normal-Block Model Object
- [`isNB()`](isNB.md) : Check if an Object is a Normal-Block Model
