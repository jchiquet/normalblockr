#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_TYPES_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_TYPES_H

#include "zi_normal_block_mean_known_clusters.h"
#include "zi_normal_block_mean_unknown_clusters.h"

// The two zero-inflated mean-block models. Mirrors normal_block_mean_types.h;
// the diagonal/spherical shape of Sigma is a runtime string here, as in the
// non-ZI family, rather than the template policy the variance-block side uses.
using zi_norm_block_mean_known_clusters   = ZINormalBlockMeanKnownClusters;
using zi_norm_block_mean_unknown_clusters = ZINormalBlockMeanUnknownClusters;

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_MEAN_TYPES_H
