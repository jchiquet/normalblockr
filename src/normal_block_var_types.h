#ifndef NORMALBLOCKR_NORMAL_BLOCK_VAR_TYPES_H
#define NORMALBLOCKR_NORMAL_BLOCK_VAR_TYPES_H

#include "noise_models.h"
#include "normal_block_var_known_clusters.h"
#include "normal_block_var_unknown_clusters.h"

// The four concrete normal-block models, obtained by crossing the clustering
// axis (known / unknown) with the residual-noise axis (diagonal / spherical).
using norm_block_var_cov_diag_noise_known_clusters        = NormalBlockVarKnownClusters<DiagonalNoise>;
using norm_block_var_cov_spherical_noise_known_clusters   = NormalBlockVarKnownClusters<SphericalNoise>;
using norm_block_var_cov_diag_noise_unknown_clusters      = NormalBlockVarUnknownClusters<DiagonalNoise>;
using norm_block_var_cov_spherical_noise_unknown_clusters = NormalBlockVarUnknownClusters<SphericalNoise>;

#endif // NORMALBLOCKR_NORMAL_BLOCK_VAR_TYPES_H
