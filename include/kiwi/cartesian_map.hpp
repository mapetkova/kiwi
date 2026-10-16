// Cartesian resampling of the finished Voronoi grid, for quick/simple
// visualization (e.g. a plain numpy slice plot -- see
// tools/plot_cartesian_map.py). NOT part of the pipeline's data contract,
// purely a convenience output.
//
// Direct port of the algorithm in the original grid-3d-cartesian-map.cc: for
// each Cartesian grid-cell center inside cfg.domain, find the enclosing
// Voronoi cell via voro++'s own nearest-site search
// (container::find_voronoi_cell(), which is exactly "the particle whose
// Voronoi cell contains this point" for voro++'s plain, unweighted
// container -- see its own docstring) and copy that cell's density/velocity
// across.
//
// This stage resamples output the Voronoi pipeline ALREADY produced, for
// visualization. It is distinct from sph-to-grid (cartesian_grid.hpp), an
// independent Cartesian mapping that runs SPH particles directly onto
// Cartesian cells, without buffer, Lloyd or tessellation.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"

namespace kiwi {

struct CartesianMap {
    int resolution = 0; // grid cells per axis
    Domain domain;       // the domain this map spans (== cfg.domain at the time it was built)

    // Flattened index for (ix, iy, iz), each in [0, resolution), x-major:
    //   idx = (ix * resolution + iy) * resolution + iz
    // Cell (ix,iy,iz)'s center is at
    //   domain.x_min + (ix+0.5)*dpx,  dpx = (domain.x_max-domain.x_min)/resolution
    // (and analogously for y/z, using the domain's own per-axis extent --
    // only exactly a cube of cells if cfg.domain itself is a cube).
    std::vector<double> density;
    std::vector<double> vx, vy, vz;
    // The Voronoi cell id that "owns" each Cartesian cell (voro++'s particle
    // id from find_voronoi_cell()). 0 (never a valid voro++ id) marks a
    // Cartesian cell where the search found no enclosing Voronoi cell at
    // all -- density/velocity are left at 0 there too.
    std::vector<int> voronoi_cell_id;
};

// Runs the point-set stage (cached) to rebuild the site positions and the
// mapping stage (cached) to get density/velocity per cell id, then resamples
// onto a uniform cfg.cartesian_map.resolution^3 Cartesian grid spanning
// cfg.domain. Parallelized with OpenMP -- see cartesian_map.cpp for why this
// needs one full voro++ container PER THREAD rather than one shared
// container (container::find_voronoi_cell() is not safe to call
// concurrently on a single container).
//
// Cached like every other stage: if cfg.stage_path("cartesian_map.dat")
// already exists, it's loaded directly instead of recomputed.
CartesianMap run_cartesian_map_stage(const Config& cfg);

// Binary format: 8-byte magic "VGTCMAP1", int32 resolution, 6 doubles
// (domain x_min,x_max,y_min,y_max,z_min,z_max), then density/vx/vy/vz
// (resolution^3 doubles each, in that order), then voronoi_cell_id
// (resolution^3 ints). Self-describing -- resolution and domain bounds are
// enough to reconstruct every cell center without any side information.
void save_cartesian_map(const CartesianMap& map, const std::string& path);
CartesianMap load_cartesian_map(const std::string& path);

} // namespace kiwi
