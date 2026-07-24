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

- [`NormalBlockData`](NormalBlockData.md) : R6 class for a generic
  normal model
- [`NormalBlockVarKnownClusters`](NormalBlockVarKnownClusters.md) : R6
  class for a normal-block model with known clustering.
- [`NormalBlockVarUnknownClusters`](NormalBlockVarUnknownClusters.md) :
  R6 class for a normal-block model with fixed number of clusters (but
  unknown clustering).
- [`ZINormalBlockVarKnownClusters`](ZINormalBlockVarKnownClusters.md) :
  R6 class for a Zero-Inflated normal-block model with a known
  clustering.
- [`ZINormalBlockVarUnknownClusters`](ZINormalBlockVarUnknownClusters.md)
  : R6 class for zero-inflated normal-block model with a fixed number of
  clusters (but unknown clustering).
- [`NormalBlockVarCollectionClusters`](NormalBlockVarCollectionClusters.md)
  : R6 class for a collection of normal-block models with different
  number of clusters (q) and a fixed sparsity level.
- [`NormalBlockVarCollectionSparsity`](NormalBlockVarCollectionSparsity.md)
  : R6 class for a collection of normal-block models with a fixed
  clustering (blocks) and different sparsity levels.
- [`NormalBlockVarCollectionClustersSparsity`](NormalBlockVarCollectionClustersSparsity.md)
  : R6 class for a collection of normal-block models with different
  number of clusters (q) and different sparsity levels.

## Data sets

- [`brca_rppa`](brca_rppa.md) : Breast cancer proteomics data (TCGA,
  RPPA)
- [`onema`](onema.md) : French stream fish community data (ONEMA / OFB
  electrofishing surveys)
- [`university`](university.md) : University webpages text data (CMU "4
  Universities" / WebKB)

## S3 methods

Standard extractors, common to every fitted model class

- [`coef(`*`<NormalBlockVarBase>`*`)`](coef.NormalBlockVarBase.md) :
  Extracts model coefficients from objects returned by
  \[NormalBlockVarBase()\] and its variants
- [`fitted(`*`<NormalBlockVarBase>`*`)`](fitted.NormalBlockVarBase.md) :
  Extracts fitted values from objects returned by
  \[NormalBlockVarBase()\] and its variants
- [`predict(`*`<NormalBlockVarBase>`*`)`](predict.NormalBlockVarBase.md)
  : Predicts observations Y for new covariates X.
- [`sigma(`*`<NormalBlockVarBase>`*`)`](sigma.NormalBlockVarBase.md) :
  Extracts the latent-block covariance matrix from objects returned by
  \[NormalBlockBase()\] and its variants

## Internal

Abstract base classes, low-level helpers and exploratory code kept for
reference; not part of the user-facing API

- [`NormalBlockVarBase`](NormalBlockVarBase.md) : R6 abstract class for
  a generic sparse Normal Block model
- [`NormalBlockVarCollection`](NormalBlockVarCollection.md) : R6
  abstract class for a collection of normal-block models
- [`SelectionNClusters`](SelectionNClusters.md) : R6 class for selecting
  the number of clusters (q) by forward/backward split-and-merge
  exploration, comparing models with the ICL
- [`get_model()`](get_model.md) : Creates appropriate new normal block
  model depending on the parametrization
- [`isNB()`](isNB.md) : Checks if a model is of class
  \[NormalBlockVarBase()\]
