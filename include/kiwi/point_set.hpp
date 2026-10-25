// Point-set stage: load particle sites, optionally add a boundary buffer
// shell, optionally relax with N Lloyd iterations, as one config-driven
// stage.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"

namespace kiwi {

// "id x y z" -- id is whatever voro++'s pre_container/container expect as a
// particle identifier (not necessarily contiguous or 0-based; whatever the
// input file used is preserved).
struct PointSet {
    std::vector<int> id;
    std::vector<double> x, y, z;

    size_t size() const { return id.size(); }
    void append(int pid, double px, double py, double pz);
};

PointSet load_point_set(const std::string& path);
void save_point_set(const PointSet& points, const std::string& path);

// Appends the graded-shells buffer points (see GradedShellsBufferParams) to `points`
// and returns the new set (does not mutate the input). Only valid for the
// "graded_shells" strategy -- see BufferConfig.
PointSet add_graded_shells_buffer(const PointSet& points, const Domain& domain,
                                  const GradedShellsBufferParams& params);

// Generic particle-coverage buffer strategy (see CoverageGridBufferParams).
// `h` must be aligned by index with `points` (i.e. h[i] is the smoothing
// length of points.id[i]/x[i]/y[i]/z[i]) -- this is only meaningful before
// any Lloyd relaxation moves the sites away from their originating
// particles, so it should be called right after extract_points(), same as
// the graded_shells call site. Appends buffer points (new sequential ids
// starting at points.size()+1, matching add_graded_shells_buffer's
// convention) and returns the new set.
PointSet add_coverage_grid_buffer(const PointSet& points, const std::vector<double>& h,
                                   const Domain& domain, const CoverageGridBufferParams& params);

// The cell-centered Cartesian candidate-point grid shared by
// add_coverage_grid_buffer() (which filters it down to just the
// SPH-uncovered candidates) and, standalone, by the lattice-based
// Voronoi regression test (tools/test_lattice_grid.cpp): lay spacing =
// spacing_fraction * min(domain extents) apart on all 3 axes, first/last
// point spacing/2 in from each wall so nothing sits exactly on a domain
// boundary. `nx`/`ny`/`nz` = round(extent/spacing) per axis (>= 1). Ids are
// sequential starting at id_offset+1 in the same i-outer/j-middle/k-inner
// nesting order used to place the points -- point (i,j,k) (0-based) gets id
// id_offset + 1 + i*ny*nz + j*nz + k, so a caller that knows nx/ny/nz can
// decode a returned id back to its grid index.
//
// NOTE for exact-cube-cell use (the lattice regression test): this only
// produces perfect cubes at every cell -- including the boundary layer,
// where voro++ clips against the domain wall instead of a neighbor's
// bisector plane -- when spacing evenly divides every domain extent (e.g.
// spacing_fraction = 1/N on a cubic domain). Otherwise round()'s residual
// leaves the outermost layer very slightly off-cube. The coverage-grid
// buffer use doesn't care about this (it only ever adds candidates, never
// needs an exact wall-to-wall tiling), but a caller that does should choose
// spacing_fraction accordingly.
struct CartesianLattice {
    PointSet points;
    double spacing = 0;
    int nx = 0, ny = 0, nz = 0;
};
CartesianLattice generate_cartesian_lattice(const Domain& domain, double spacing_fraction,
                                             int id_offset = 0);

// One Lloyd relaxation step: build a Voronoi tessellation of `points` inside
// `domain` (using `blocks` for voro++'s spatial hashing), then move every
// site to its cell's centroid.
PointSet lloyd_step(const PointSet& points, const Domain& domain, const GridBlocks& blocks);

// Extracts points from cfg.input_snapshot, then runs buffer (if enabled)
// and N Lloyd iterations per cfg, writing intermediate files under
// cfg.run_dir() (and every iteration if cfg.lloyd.save_each_iteration is
// set). Returns the final point set. This is the single entry point for
// the point-set stage -- no other code should reimplement this
// sequencing.
//
// If cfg.stage_path("points-final.dat") already exists on disk, it is
// loaded directly instead of being recomputed -- Lloyd relaxation is slow
// (~1 min/iteration at ~450k points), so re-running an unchanged stage on
// every invocation of a later stage would make iterating on the pipeline
// painfully slow. Delete the file (or the whole run directory) to force a
// recompute.
PointSet run_point_set_stage(const Config& cfg);

} // namespace kiwi
