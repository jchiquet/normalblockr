#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_TYPES_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_TYPES_H

#include "zi_normal_block_known_clusters.h"
#include "zi_normal_block_unknown_clusters.h"
#include "zi_noise_models.h"

// The zero-inflated counterparts of the four types in normal_block_types.h,
// crossing the same two axes (known / unknown clustering, diagonal /
// spherical residual noise).
using zi_norm_block_cov_diag_noise_known_clusters        = ZINormalBlockKnownClusters<ZIDiagonalNoise>;
using zi_norm_block_cov_spherical_noise_known_clusters   = ZINormalBlockKnownClusters<ZISphericalNoise>;
using zi_norm_block_cov_diag_noise_unknown_clusters      = ZINormalBlockUnknownClusters<ZIDiagonalNoise>;
using zi_norm_block_cov_spherical_noise_unknown_clusters = ZINormalBlockUnknownClusters<ZISphericalNoise>;

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_TYPES_H
