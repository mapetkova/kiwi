// Shared helper for the three places (tessellate_voronoi, lloyd_step,
// run_mapping_stage) that build a voro++ container and need a block count.
// See GridBlocks::block_scale_factor (config.hpp) for the full story: voro++'s
// own pre_container_base::guess_optimal() sizes its spatial hash for a
// domain-uniform particle density and completely ignores whatever nx/ny/nz
// are passed in -- this helper is the one place block_scale_factor actually
// takes effect, applied uniformly across all three call sites so they can't
// drift out of sync.
#pragma once

#include <algorithm>
#include <cmath>

#include "voro++.hh"

namespace kiwi {

inline void guess_scaled_blocks(voro::pre_container_base& pcon, double scale_factor,
                                 int& nx, int& ny, int& nz) {
    pcon.guess_optimal(nx, ny, nz);
    if (scale_factor != 1.0) {
        nx = std::max(1, static_cast<int>(std::round(nx * scale_factor)));
        ny = std::max(1, static_cast<int>(std::round(ny * scale_factor)));
        nz = std::max(1, static_cast<int>(std::round(nz * scale_factor)));
    }
}

} // namespace kiwi
