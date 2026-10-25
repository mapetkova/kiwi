// Benchmark harness for the grid_blocks setting (config.hpp, GridBlocks).
// Not part of the kiwi CLI itself.
//
// guess_optimal(nx,ny,nz) (third_party/voro/src/pre_container.cc) completely
// overwrites whatever nx/ny/nz it's given -- see GridBlocks in config.hpp.
// This harness bypasses that to
// actually test whether voro++'s own heuristic (a domain-wide average of
// optimal_particles=5.6 per block, config.hh) is well matched to our workload,
// which is NOT domain-uniform: a real SPH cloud, often clumpy/filamentary, plus
// (optionally) a uniform coverage_grid buffer lattice filling the rest of the
// domain.
//
// Usage: bench_grid_blocks <config.json> [repeats] [comma,separated,scales]
// Loads the config's cached points-final.dat (via run_point_set_stage, so this
// does NOT recompute buffer/Lloyd if already cached) and domain, then times a
// full tessellation pass (identical inner loop to tessellate_voronoi) across a
// range of block counts: several multiples of guess_optimal's own choice.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>

#include "voro++.hh"
#include "kiwi/config.hpp"
#include "kiwi/point_set.hpp"

using namespace voro;
using Clock = std::chrono::steady_clock;

namespace {

// One timed tessellation pass with a FORCED block count (no guess_optimal
// call) -- otherwise byte-for-byte the same container setup/loop as
// tessellate_voronoi() in tessellation.cpp.
double time_tessellation(const kiwi::PointSet& points, const kiwi::Domain& d,
                          int nx, int ny, int nz) {
    pre_container pcon(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i) {
        pcon.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    }

    container con(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, nx, ny, nz, false, false, false, 8);
    pcon.setup(con);

    auto t0 = Clock::now();

    voronoicell_neighbor c;
    std::vector<int> n_list;
    c_loop_all clo(con);
    double x0, y0, z0;
    size_t n_cells = 0;
    if (clo.start()) do {
        if (con.compute_cell(c, clo)) {
            clo.pos(x0, y0, z0);
            c.neighbors(n_list);
            ++n_cells;
        }
    } while (clo.inc());

    auto t1 = Clock::now();
    (void)n_cells;
    return std::chrono::duration<double>(t1 - t0).count();
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: bench_grid_blocks <config.json> [repeats]\n");
        return 1;
    }
    int repeats = argc >= 3 ? std::atoi(argv[2]) : 3;

    kiwi::Config cfg = kiwi::Config::load_from_file(argv[1]);
    kiwi::PointSet points = kiwi::run_point_set_stage(cfg); // cached load if present
    std::fprintf(stderr, "loaded %zu points from %s\n", points.size(), cfg.run_dir().c_str());

    // Reproduce guess_optimal's own choice for reference.
    pre_container pcon_ref(cfg.domain.x_min, cfg.domain.x_max, cfg.domain.y_min, cfg.domain.y_max,
                            cfg.domain.z_min, cfg.domain.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i)
        pcon_ref.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    int gx = 0, gy = 0, gz = 0;
    pcon_ref.guess_optimal(gx, gy, gz);
    std::fprintf(stderr, "guess_optimal() chose nx=%d ny=%d nz=%d (%d blocks, %.2f particles/block avg)\n",
                 gx, gy, gz, gx * gy * gz, double(points.size()) / (gx * gy * gz));

    std::vector<double> scales = {0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0};
    if (argc >= 4) {
        scales.clear();
        std::stringstream ss(argv[3]);
        std::string tok;
        while (std::getline(ss, tok, ',')) scales.push_back(std::atof(tok.c_str()));
    }
    std::printf("scale\tnx\tny\tnz\tblocks\tavg_pts/block\tmin_s\tmedian_s\n");
    for (double s : scales) {
        int nx = std::max(1, int(std::round(gx * s)));
        int ny = std::max(1, int(std::round(gy * s)));
        int nz = std::max(1, int(std::round(gz * s)));
        std::vector<double> times;
        for (int r = 0; r < repeats; ++r) {
            times.push_back(time_tessellation(points, cfg.domain, nx, ny, nz));
        }
        std::sort(times.begin(), times.end());
        double min_t = times.front();
        double med_t = times[times.size() / 2];
        long blocks = long(nx) * ny * nz;
        std::printf("%.2f\t%d\t%d\t%d\t%ld\t%.3f\t%.4f\t%.4f\n",
                    s, nx, ny, nz, blocks, double(points.size()) / blocks, min_t, med_t);
        std::fflush(stdout);
    }
    return 0;
}
