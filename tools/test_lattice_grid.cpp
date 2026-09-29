// Permanent regression test: a Voronoi tessellation of a simple cubic
// lattice of generating points is, exactly, a partition of the domain into
// cubes -- each cell bounded by the perpendicular bisector planes to its
// +/-x/+/-y/+/-z neighbors (interior cells) or by the domain wall itself
// (boundary cells, which come out as perfect cubes too as long as the
// domain wall sits exactly spacing/2 beyond the outermost lattice point on
// every axis -- see generate_cartesian_lattice()'s header comment).
//
// This exercises exactly the tessellation/block-count/voro++-container
// machinery every real run depends on, against a case with a fully known
// analytic answer: every cell has volume == spacing^3, and every cell has
// exactly 6 neighbor-list entries (a mix of real positive-id neighbors and
// negative domain-wall ids, per voro++'s convention -- see CellGeometry),
// with each real neighbor's center offset from this cell's center by
// exactly `spacing` along exactly one axis. No SPH data, no buffer, no
// Lloyd relaxation -- the lattice points themselves are the full site set.
//
// This is also the ground-truth reference fixture for validating
// "development one" (the dedicated Cartesian-mapping rewrite, see the
// project notes): the same lattice, fed through this tool's existing
// (already-validated) kernel-integral mapping machinery, gives known-correct
// density/velocity values to check the new fast path against.
//
// Usage: test_lattice_grid [N] [extent] [block_scale_factor]
//   N       -- lattice points per axis (domain is an NxNxN cube). Default 8.
//   extent  -- domain side length. Default 10.0.
//   block_scale_factor -- passed through to guess_scaled_blocks(); doesn't
//                         affect correctness, only voro++'s search speed.
//                         Default 3.0 (this tool's own default).
//
// spacing_fraction is deliberately fixed at 1.0/N (not a CLI knob) --
// spacing = extent/N divides the domain exactly, which is the whole point:
// anything else leaves the outermost layer of cells very slightly
// off-cube (see generate_cartesian_lattice()'s header comment), which
// would make this test's own analytic answer wrong, not the code under
// test.
//
// Exit code 0 on a fully clean pass, 1 on any assertion failure -- suitable
// as an actual regression gate, not just an eyeballed report.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>

#include "kiwi/config.hpp"
#include "kiwi/point_set.hpp"
#include "kiwi/tessellation.hpp"

using namespace kiwi;

namespace {

int decode_i(int id0based, int ny, int nz) { return id0based / (ny * nz); }
int decode_j(int id0based, int ny, int nz) { return (id0based / nz) % ny; }
int decode_k(int id0based, int nz) { return id0based % nz; }

} // namespace

int main(int argc, char** argv) {
    int N = argc >= 2 ? std::atoi(argv[1]) : 8;
    double extent = argc >= 3 ? std::atof(argv[2]) : 10.0;
    double block_scale_factor = argc >= 4 ? std::atof(argv[3]) : 3.0;

    if (N < 2) {
        std::fprintf(stderr, "test_lattice_grid: N must be >= 2 (need at least one interior-vs-boundary "
                              "distinction to be meaningful)\n");
        return 1;
    }

    Domain d;
    d.x_min = 0.0; d.x_max = extent;
    d.y_min = 0.0; d.y_max = extent;
    d.z_min = 0.0; d.z_max = extent;

    const double spacing_fraction = 1.0 / N;
    CartesianLattice lat = generate_cartesian_lattice(d, spacing_fraction);

    std::fprintf(stderr, "test_lattice_grid: N=%d extent=%.6g -> nx=%d ny=%d nz=%d spacing=%.10g\n",
                 N, extent, lat.nx, lat.ny, lat.nz, lat.spacing);

    const long expected_count = static_cast<long>(lat.nx) * lat.ny * lat.nz;
    if (lat.nx != N || lat.ny != N || lat.nz != N) {
        std::fprintf(stderr, "test_lattice_grid: FAIL -- expected nx=ny=nz=%d, got %d/%d/%d "
                              "(spacing_fraction=1/N on a cube domain should always give this exactly)\n",
                     N, lat.nx, lat.ny, lat.nz);
        return 1;
    }
    if (static_cast<long>(lat.points.size()) != expected_count) {
        std::fprintf(stderr, "test_lattice_grid: FAIL -- generated %zu lattice points, expected %ld\n",
                     lat.points.size(), expected_count);
        return 1;
    }

    GridBlocks blocks;
    blocks.block_scale_factor = block_scale_factor;
    Tessellation tess = tessellate_voronoi(lat.points, d, blocks);

    int n_failures = 0;
    auto fail = [&](const std::string& msg) {
        std::fprintf(stderr, "test_lattice_grid: FAIL -- %s\n", msg.c_str());
        ++n_failures;
    };

    if (static_cast<long>(tess.cells.size()) != expected_count) {
        fail("tessellation produced " + std::to_string(tess.cells.size()) + " cells, expected " +
             std::to_string(expected_count) + " (every lattice point should yield exactly one cell -- "
             "no degenerate/missing cells are possible on an exact grid)");
    }

    // id -> cell, for the neighbor-position cross-check below. Ids are
    // 1-based sequential in i-outer/j-middle/k-inner order (see
    // generate_cartesian_lattice()).
    std::map<int, const CellGeometry*> by_id;
    for (const auto& c : tess.cells) by_id[c.id] = &c;

    const double expected_volume = lat.spacing * lat.spacing * lat.spacing;
    const double vol_tol = 1e-9 * expected_volume; // relative, generous vs. observed ~1e-13 fp noise
    const double pos_tol = 1e-9 * lat.spacing;

    double max_vol_rel_err = 0.0;
    long n_wall_total = 0, n_real_neighbor_total = 0;

    for (const auto& c : tess.cells) {
        const double vol_err = std::abs(c.volume - expected_volume);
        max_vol_rel_err = std::max(max_vol_rel_err, vol_err / expected_volume);
        if (vol_err > vol_tol) {
            fail("cell id=" + std::to_string(c.id) + " volume=" + std::to_string(c.volume) +
                 " != expected cube volume " + std::to_string(expected_volume));
        }

        if (c.neighbor_ids.size() != 6) {
            fail("cell id=" + std::to_string(c.id) + " has " + std::to_string(c.neighbor_ids.size()) +
                 " neighbor-list entries, expected exactly 6 (a cube always has 6 faces, each either a "
                 "real neighbor or a domain wall)");
            continue;
        }

        // Expected wall-face count from this cell's (i,j,k) grid index: one
        // per axis where it sits on the domain boundary (index 0 or N-1).
        const int idx0 = c.id - 1;
        const int i = decode_i(idx0, lat.ny, lat.nz);
        const int j = decode_j(idx0, lat.ny, lat.nz);
        const int k = decode_k(idx0, lat.nz);
        int expected_walls = 0;
        expected_walls += (i == 0 || i == lat.nx - 1) ? 1 : 0;
        expected_walls += (j == 0 || j == lat.ny - 1) ? 1 : 0;
        expected_walls += (k == 0 || k == lat.nz - 1) ? 1 : 0;

        int actual_walls = 0, actual_real = 0;
        for (int nid : c.neighbor_ids) {
            if (nid < 0) {
                ++actual_walls;
            } else {
                ++actual_real;
                auto it = by_id.find(nid);
                if (it == by_id.end()) {
                    fail("cell id=" + std::to_string(c.id) + " neighbor id=" + std::to_string(nid) +
                         " is not in the lattice at all");
                    continue;
                }
                const CellGeometry& nb = *it->second;
                const double ddx = nb.x - c.x, ddy = nb.y - c.y, ddz = nb.z - c.z;
                // Exactly one axis should differ by +/-spacing, the other two by ~0.
                int axes_off = 0;
                double total_dist = std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
                for (double comp : {ddx, ddy, ddz}) {
                    if (std::abs(std::abs(comp) - lat.spacing) < pos_tol) ++axes_off;
                }
                if (axes_off != 1 || std::abs(total_dist - lat.spacing) > pos_tol) {
                    fail("cell id=" + std::to_string(c.id) + " neighbor id=" + std::to_string(nid) +
                         " center offset (" + std::to_string(ddx) + "," + std::to_string(ddy) + "," +
                         std::to_string(ddz) + ") is not exactly one lattice spacing along one axis");
                }
            }
        }
        n_wall_total += actual_walls;
        n_real_neighbor_total += actual_real;

        if (actual_walls != expected_walls) {
            fail("cell id=" + std::to_string(c.id) + " (i,j,k)=(" + std::to_string(i) + "," +
                 std::to_string(j) + "," + std::to_string(k) + ") has " + std::to_string(actual_walls) +
                 " domain-wall neighbor entries, expected " + std::to_string(expected_walls) +
                 " from its grid position");
        }
    }

    std::fprintf(stderr, "test_lattice_grid: %zu cells, max cell-volume relative error %.3e, "
                          "%ld total wall-face entries, %ld total real-neighbor entries\n",
                 tess.cells.size(), max_vol_rel_err, n_wall_total, n_real_neighbor_total);

    if (n_failures == 0) {
        std::fprintf(stderr, "test_lattice_grid: PASS (%d x %d x %d lattice, spacing=%.10g, "
                              "every cell an exact cube, full neighbor topology verified)\n",
                     lat.nx, lat.ny, lat.nz, lat.spacing);
        return 0;
    } else {
        std::fprintf(stderr, "test_lattice_grid: %d FAILURE(S)\n", n_failures);
        return 1;
    }
}
