// sph-to-grid: an independent SPH -> Cartesian grid mapping path --
// computed directly, without going through voro++, point_set, tessellation,
// or the Voronoi mapping stage at all. Unlike cartesian-map
// (cartesian_map.hpp), which is a *resample* of the Voronoi pipeline's own
// already-computed output, for visualization convenience, this computes
// density and mass-weighted velocity straight onto a regular Cartesian
// grid -- no buffer, no Lloyd relaxation, no Voronoi tessellation anywhere
// in the path.
//
// Modeled directly on SPLASH's interpolate3D.f90 (github.com/danieljprice/
// splash, "exact_rendering" branch -- Petkova, Laibe & Bonnell 2018, J.
// Comp. Phys, "exact overlap of kernel with wall boundaries"):
//
//  - For each SPH particle, the range of grid-cell indices its kernel could
//    possibly reach is computed by plain arithmetic from the particle's
//    position and smoothing length (kernel support radius 2h), the same way
//    interpolate3D.f90 computes ipixmin/ipixmax/etc. No neighbor search or
//    BFS discovery is needed here at all, unlike the Voronoi mapping stage
//    (mapping.cpp), which has to discover reachable Voronoi cells one face
//    at a time because they aren't arranged on a regular lattice -- a
//    Cartesian grid's neighbor structure is already known in closed form.
//    The index range is a bounding BOX around the kernel's support sphere,
//    not trimmed to the sphere itself (same as SPLASH): cells in the box's
//    corners that the sphere doesn't actually reach just evaluate to ~0
//    contribution rather than being explicitly excluded.
//
//  - Each reachable cell's exact kernel-overlap integral is evaluated with
//    KIWI's kernel-integral machinery
//    (kernel_integral.hpp's compute_face_integral(), compute_wall_constants(),
//    compute_edge_constants(), vertex_integral()) -- unchanged math, just
//    applied to a cell's 6 fixed 4-corner axis-aligned rectangular faces
//    instead of walking an arbitrary voro++ polygon face. This is the same
//    underlying closed-form kernel line integral family as SPLASH's own
//    wallint()/pint3D(), just specialized to a known rectangle rather than
//    generalized to an arbitrary polygon -- see compute_face_integral()'s
//    own doc comment for the winding convention this relies on.
//
// The output format is a placeholder: it is modelled on cartesian_map.hpp's
// self-describing binary format (without the voronoi_cell_id field, which
// doesn't apply here -- there is no Voronoi cell in this path) so that there
// is a runnable, inspectable artifact end-to-end. A POLARIS-native Cartesian
// grid format and a plain self-describing ASCII writer have not been
// designed yet.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"

namespace kiwi {

struct CartesianGridResult {
    // Cells per axis. Equal to cfg.cartesian_grid.resolution on all 3 axes
    // unless resolution_x/resolution_y/resolution_z overrode one -- see
    // CartesianGridConfig's doc comment (config.hpp). Kept as 3 explicit
    // fields (not a single `resolution`) since that override is the whole
    // point: this result can legitimately describe a non-cube-count grid.
    int nx = 0, ny = 0, nz = 0;
    Domain domain; // == cfg.domain at the time this was built

    // Flattened index for (ix, iy, iz), ix in [0,nx), iy in [0,ny), iz in
    // [0,nz), x-major: idx = (ix * ny + iy) * nz + iz -- same layout
    // convention as CartesianMap (cartesian_map.hpp) generalized to
    // independent per-axis counts (for nx==ny==nz==res it is
    // `(ix*res+iy)*res+iz`). Cell (ix,iy,iz)'s center is at
    // domain.x_min + (ix+0.5)*dpx, dpx = (domain.x_max-domain.x_min)/nx (and
    // analogously for y/z). NOTE: this is NOT the same element order as
    // SPLASH's own `splash to grid` output files (those are written
    // z-outer/y-middle/x-inner, i.e. x fastest) -- a comparison script needs
    // to index each format by its own convention, not assume they match.
    std::vector<double> density;
    std::vector<double> vx, vy, vz;
};

// Loads cfg.input_snapshot directly (no point-set/tessellation/mapping
// stage involved at all -- this path doesn't touch voro++) and computes
// density + mass-weighted velocity on a uniform nx x ny x nz grid spanning
// cfg.domain (nx/ny/nz per CartesianGridConfig -- see config.hpp). Parallelized
// with OpenMP over particles, same per-thread-accumulator-then-reduce
// pattern as run_mapping_stage() (see mapping.cpp) -- each thread gets a
// full-size accumulator grid rather than locking a shared one, since nearby
// particles routinely hit the same handful of cells and contention would
// erase most of the parallel speedup.
//
// Cached like every other stage: if cfg.stage_path("cartesian_grid.dat")
// already exists, it's loaded directly instead of recomputed.
CartesianGridResult run_sph_to_grid_stage(const Config& cfg);

// Binary format v2: 8-byte magic "VGTCGRD2", int32 nx, int32 ny, int32 nz,
// 6 doubles (domain bounds: x_min,x_max,y_min,y_max,z_min,z_max), then
// density/vx/vy/vz (nx*ny*nz doubles each, in that order, x-major -- see
// CartesianGridResult's own doc comment for the exact index formula).
// Placeholder -- see this file's header comment; not the real POLARIS
// Cartesian output nor the plain self-describing ASCII format, both still
// to be designed. Older "VGTCGRD1" files (a single int32 resolution) are
// not readable by this loader.
void save_cartesian_grid(const CartesianGridResult& result, const std::string& path);
CartesianGridResult load_cartesian_grid(const std::string& path);

} // namespace kiwi
