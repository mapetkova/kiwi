// Tessellation stage: turn a point set into per-cell geometry (position,
// volume, neighbor list). Deliberately knows nothing about density,
// velocity, or POLARIS -- this is the seam where a Cartesian (or other)
// grid backend can be swapped in later without touching the mapping or
// writer stages.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"
#include "kiwi/point_set.hpp"

namespace kiwi {

struct CellGeometry {
    int id = 0;
    double x = 0, y = 0, z = 0; // cell center (site position after any relaxation)
    double volume = 0;
    std::vector<int> neighbor_ids; // voro++ convention: negative ids denote domain walls
};

struct Tessellation {
    std::vector<CellGeometry> cells;
};

// Builds a Voronoi tessellation of `points` inside `domain`. This is a
// direct port of the con.compute_cell/c_loop_all sequence in grid-file.cc,
// factored out so it can be called independently of reading/writing any
// particular file format.
Tessellation tessellate_voronoi(const PointSet& points, const Domain& domain, const GridBlocks& blocks);

// Parses the format save_tessellation_summary() writes ("id x y z volume
// n_neighbors neighbor_ids..." per line) back into a Tessellation.
Tessellation load_tessellation(const std::string& path);

// Runs the point-set stage (cached) to get the final grid, then tessellates
// it -- unless cfg.stage_path("tessellation.dat") already exists, in which
// case it's parsed and returned directly instead of recomputed. Same
// caching pattern as run_point_set_stage()/run_mapping_stage(): every
// caller that only needs position/volume/neighbor-topology (tessellate,
// write-polaris, write-voronoi) should go through this rather than calling
// tessellate_voronoi() directly, so a full voro++ container build+cell walk
// happens at most once per run directory instead of once per command.
//
// Not used by the mapping stage: run_mapping_stage() needs full per-cell
// face-vertex geometry (CellFaceGeometry) for the kernel integral, a
// strictly richer dataset than this stage's CellGeometry summary captures,
// so it still does its own separate voro++ build.
Tessellation run_tessellation_stage(const Config& cfg);

// Diagnostics mirroring grid-file.cc's own printf checks: cells with
// non-positive volume, and the count of input points that didn't produce a
// computed cell at all (voro++'s compute_cell() can return false for a
// point in a degenerate configuration -- e.g. at a domain wall during Lloyd
// relaxation).
struct TessellationDiagnostics {
    std::vector<int> zero_or_negative_volume_ids;
    std::vector<int> missing_cell_ids; // present in input, absent from output
};
TessellationDiagnostics diagnose(const PointSet& input, const Tessellation& tess);

void save_tessellation_summary(const Tessellation& tess, const std::string& path);

} // namespace kiwi
