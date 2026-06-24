#ifndef NORMALBLOCKR_NORMAL_BLOCK_TYPES_H
#define NORMALBLOCKR_NORMAL_BLOCK_TYPES_H

#include "noise_models.h"
#include "normal_block_known_clusters.h"
#include "normal_block_unknown_clusters.h"

// The four concrete normal-block models, obtained by crossing the clustering
// axis (known / unknown) with the residual-noise axis (diagonal / spherical).
using norm_block_cov_diag_noise_known_clusters        = NormalBlockKnownClusters<DiagonalNoise>;
using norm_block_cov_spherical_noise_known_clusters   = NormalBlockKnownClusters<SphericalNoise>;
using norm_block_cov_diag_noise_unknown_clusters      = NormalBlockUnknownClusters<DiagonalNoise>;
using norm_block_cov_spherical_noise_unknown_clusters = NormalBlockUnknownClusters<SphericalNoise>;

#endif // NORMALBLOCKR_NORMAL_BLOCK_TYPES_H
