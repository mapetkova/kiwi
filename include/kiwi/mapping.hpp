// Maps SPH particle density and mass-weighted velocity onto the
// final Voronoi grid (post buffer + Lloyd).
//
// A single breadth-first walk over the tessellation evaluates the expensive
// SPH kernel line integral once per (particle, cell) pair and accumulates
// both outputs -- mass for density, mass-weighted velocity -- from that one
// value. No smoothing is applied to the result.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"

namespace kiwi {

struct CellMapping {
    int id = 0;
    double x = 0, y = 0, z = 0;
    double volume = 0;
    double density = 0;
    double vx = 0, vy = 0, vz = 0;
};

struct MappingResult {
    std::vector<CellMapping> cells;
};

void save_mapping(const MappingResult& result, const std::string& path);
MappingResult load_mapping(const std::string& path);

// Runs the point-set stage (cached) to get the final grid, tessellates it
// while retaining full per-cell face geometry (needed for the kernel
// integral -- unlike tessellation.hpp's CellGeometry, which only keeps
// volume + neighbor ids), then walks every SPH particle in cfg.input_snapshot
// and BFS-accumulates density and mass-weighted velocity together.
//
// Cached like the other stages: if cfg.stage_path("mapping.dat") already
// exists, it's loaded directly instead of recomputed.
MappingResult run_mapping_stage(const Config& cfg);

// Optional validation step for a future Voronoi-grid -> SPH-particle inverse
// mapping, adapted from CMacIonize's Petkova ionic-fraction inverse mapping
// (SPHArrayInterface::InverseMappingFunction, SPHARRAY_MAPPING_PETKOVA
// branch; Petkova et al. 2021) to particle mass instead. Restricted to the
// "testing step" case: recovering each particle's mass
// from the SAME grid that particle distribution already forward-mapped onto
// via run_mapping_stage() -- NOT a grid-agnostic standalone inverse mapping
// (that needs a second BFS pass; see mapping.cpp for why, and the "second
// pass" discussion is deliberately deferred to a later, separate feature).
//
// One particle is one row: m_original (from the snapshot), m_recovered (the
// per-cell-normalized weighted sum defined below), and coverage_fraction =
// m_recovered / m_original (undefined, reported as 0, if m_original <= 0).
//
// The math, and an important caveat about what this specific test can and
// cannot catch: for each cell i, run_mapping_stage() already computed
// M_i = density_i * volume_i = sum over ALL particles b of K_{b,i} * m_b,
// where K_{b,i} is the same normalized SPH kernel wall-integral this file's
// compute_cell_integral()/compute_face_integral() already evaluates walking
// outward from particle b (0 if cell i is outside particle b's kernel
// support). CMacIonize's actual shipped code (not quite what its own paper's
// eq. 8 states -- see the design discussion this stage grew out of) weights
// particle a's recovered share of cell i's property P_i by
// weight_{a,i} = K_{a,i} * m_a / M_i, i.e. recovered_a = sum_i weight_{a,i} * P_i.
// This stage sets P_i = M_i itself (recovering mass from the mass/density
// grid) -- and M_i then cancels exactly: recovered_a = sum_i K_{a,i} * m_a =
// m_a * coverage_fraction_a. So this test is mathematically guaranteed to
// reduce to "how much of particle a's kernel support did the BFS actually
// reach" (1.0 for a fully-covered interior particle; less than 1.0 near a
// domain boundary or where the grid/BFS truncates the kernel early) --
// it CANNOT detect a bug in run_mapping_stage()'s own density values,
// because whatever M_i (right or wrong) is read back cancels out of the
// weight*property product identically. It is still a genuine, useful check
// of kernel/BFS support completeness and of the aggregate mass identity
// (sum_a recovered_a == sum_i M_i exactly, by construction, for any masses)
// -- just not a check that the forward mapping computed the right density.
//
// Requires run_mapping_stage() to have already produced (or caches here
// trigger it to produce) cfg.stage_path("mapping.dat"); does its own fresh
// voro++ tessellation (geometry isn't persisted between stages) via the same
// build_cell_geometry() helper run_mapping_stage() uses, and looks up each
// cell's M_i by its voro++ site id (not by array position) so it never has
// to assume this run's tessellation enumerates cells in the same order the
// original forward-mapping run did.
//
// Cached like every other stage: if cfg.stage_path("inverse-mapping-test.dat")
// already exists, it's loaded directly instead of recomputed.
struct InverseMappingTestRow {
    int id = 0;
    double m_original = 0;
    double m_recovered = 0;
    double coverage_fraction = 0;
};

struct InverseMappingTestResult {
    std::vector<InverseMappingTestRow> particles;
    double total_original_mass = 0;   // sum of m_original over all particles
    double total_grid_mass = 0;       // sum of M_i over all cells (from the forward mapping)
    double total_recovered_mass = 0;  // sum of m_recovered over all particles -- == total_grid_mass exactly, by construction
};

void save_inverse_mapping_test(const InverseMappingTestResult& result, const std::string& path);
InverseMappingTestResult load_inverse_mapping_test(const std::string& path);

InverseMappingTestResult run_inverse_mapping_test_stage(const Config& cfg);

} // namespace kiwi
