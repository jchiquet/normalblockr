#ifndef NORMALBLOCKR_NORMAL_BLOCK_MEAN_TYPES_H
#define NORMALBLOCKR_NORMAL_BLOCK_MEAN_TYPES_H

#include "normal_block_mean_known_clusters.h"
#include "normal_block_mean_unknown_clusters.h"

// The two concrete mean-block models. Mirrors normal_block_var_types.h, minus
// the residual-noise axis: the mean-block family carries a full p x p Sigma
// rather than a diagonal/spherical D.
using norm_block_mean_known_clusters   = NormalBlockMeanKnownClusters;
using norm_block_mean_unknown_clusters = NormalBlockMeanUnknownClusters;

#endif // NORMALBLOCKR_NORMAL_BLOCK_MEAN_TYPES_H
