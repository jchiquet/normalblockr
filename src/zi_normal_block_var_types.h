#ifndef NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_TYPES_H
#define NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_TYPES_H

#include "zi_normal_block_var_known_clusters.h"
#include "zi_normal_block_var_unknown_clusters.h"
#include "zi_noise_models.h"

// The zero-inflated counterparts of the four types in normal_block_types.h,
// crossing the same two axes (known / unknown clustering, diagonal /
// spherical residual noise).
using zi_norm_block_var_cov_diag_noise_known_clusters        = ZINormalBlockVarKnownClusters<ZIDiagonalNoise>;
using zi_norm_block_var_cov_spherical_noise_known_clusters   = ZINormalBlockVarKnownClusters<ZISphericalNoise>;
using zi_norm_block_var_cov_diag_noise_unknown_clusters      = ZINormalBlockVarUnknownClusters<ZIDiagonalNoise>;
using zi_norm_block_var_cov_spherical_noise_unknown_clusters = ZINormalBlockVarUnknownClusters<ZISphericalNoise>;

#endif // NORMALBLOCKR_ZI_NORMAL_BLOCK_VAR_TYPES_H
