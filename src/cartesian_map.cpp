#include "kiwi/cartesian_map.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "voro++.hh"
#include "kiwi/mapping.hpp"
#include "kiwi/point_set.hpp"
#include "kiwi/voro_grid.hpp"

namespace {
bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}
} // namespace

namespace kiwi {

namespace {

// Builds a fresh voro++ container from `points`, identical to
// tessellate_voronoi()'s own container-build step (pre_container -> scaled
// block count -> container::setup) but without walking cells afterward --
// callers here only ever do point-location queries (find_voronoi_cell), not
// a full cell computation pass.
//
// Deliberately NOT shared as a single object across OpenMP threads.
// container::find_voronoi_cell() delegates to voro_compute::find_voronoi_cell()
// (v_compute.cc), which reads and writes a `mask` scratch array owned by the
// container's single internal voro_compute instance -- allocated once in the
// constructor and reused across calls (to avoid a per-call allocation) to
// track which blocks of the masked search have already been tested during
// THIS call. That reuse is exactly what makes concurrent calls from
// multiple threads a genuine data race: two threads searching at the same
// time would read and write the same mask array, corrupting each other's
// in-progress search. Unlike the mapping stage's BFS (which already gives
// every thread its own accumulator/scratch arrays, see mapping.cpp), there
// is no way to make a single container's find_voronoi_cell() calls
// thread-safe without either locking around the whole search (which would
// serialize it and defeat the point of parallelizing) or giving every
// thread its own copy of the mutable search state.
//
// This takes the simpler of those two options: each thread gets its own
// FULL container, built from the same points. That costs more memory than
// the mapping stage's per-thread scratch (each copy holds the whole point
// set, not just a small accumulator), but in exchange the search itself is
// then genuinely lock-free -- and because it's a pure nearest-site lookup
// with no floating-point summation order to vary, every thread's container
// finds the same answer for the same query point regardless of how many
// threads are running.
std::unique_ptr<voro::container> build_container(const PointSet& points, const Domain& d,
                                                   const GridBlocks& blocks) {
    using namespace voro;
    pre_container pcon(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i) {
        pcon.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    }
    int nx, ny, nz;
    guess_scaled_blocks(pcon, blocks.block_scale_factor, nx, ny, nz);
    auto con = std::make_unique<container>(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max,
                                            nx, ny, nz, false, false, false, 8);
    pcon.setup(*con);
    return con;
}

constexpr char kMagic[8] = {'V', 'G', 'T', 'C', 'M', 'A', 'P', '1'};

} // namespace

CartesianMap run_cartesian_map_stage(const Config& cfg) {
    std::string final_path = cfg.stage_path("cartesian_map.dat");
    if (file_exists(final_path)) {
        return load_cartesian_map(final_path);
    }

    const auto& d = cfg.domain;
    {
        // Not a hard requirement the way it is for the POLARIS writer (this
        // output has no external format baking in a single side length) --
        // but a non-cube domain with a single per-axis resolution means the
        // Cartesian cells themselves aren't cubes either, worth flagging
        // since that's easy to not notice when just eyeballing a plot.
        const double dx = d.x_max - d.x_min;
        const double dy = d.y_max - d.y_min;
        const double dz = d.z_max - d.z_min;
        const double max_extent = std::max({dx, dy, dz});
        const double tol = 1e-6 * max_extent;
        if (std::abs(dx - dy) > tol || std::abs(dx - dz) > tol || std::abs(dy - dz) > tol) {
            std::cerr << "[kiwi] warning: domain is not a cube (x extent=" << dx << ", y extent="
                      << dy << ", z extent=" << dz << ") -- the Cartesian map uses a single "
                         "resolution for all 3 axes, so its cells will not be cubes (density/"
                         "velocity values themselves are still correct, only the implied cell "
                         "shape is skewed).\n";
        }
    }

    int res = cfg.cartesian_map.resolution;
    if (res <= 0) throw std::runtime_error("cartesian_map.resolution must be > 0");

    PointSet points = run_point_set_stage(cfg);
    MappingResult mapping = run_mapping_stage(cfg);

    std::unordered_map<int, const CellMapping*> mapping_by_id;
    mapping_by_id.reserve(mapping.cells.size() * 2);
    for (const auto& m : mapping.cells) mapping_by_id[m.id] = &m;

    CartesianMap map;
    map.resolution = res;
    map.domain = d;
    const size_t n = static_cast<size_t>(res) * static_cast<size_t>(res) * static_cast<size_t>(res);
    map.density.assign(n, 0.0);
    map.vx.assign(n, 0.0);
    map.vy.assign(n, 0.0);
    map.vz.assign(n, 0.0);
    map.voronoi_cell_id.assign(n, 0);

    const double dpx = (d.x_max - d.x_min) / res;
    const double dpy = (d.y_max - d.y_min) / res;
    const double dpz = (d.z_max - d.z_min) / res;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    // One full container per thread -- see build_container()'s comment above
    // for why this is necessary rather than sharing one container.
    std::vector<std::unique_ptr<voro::container>> containers;
    containers.reserve(static_cast<size_t>(n_threads));
    for (int t = 0; t < n_threads; ++t) {
        containers.push_back(build_container(points, d, cfg.grid_blocks));
    }

    long unresolved = 0; // find_voronoi_cell() found no enclosing cell at all
    long unmapped = 0;   // found a cell id, but that id has no mapping.dat entry

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic) reduction(+ : unresolved, unmapped)
#endif
    for (int ix = 0; ix < res; ++ix) {
        for (int iy = 0; iy < res; ++iy) {
            int tid = 0;
#ifdef _OPENMP
            tid = omp_get_thread_num();
#endif
            voro::container& con = *containers[static_cast<size_t>(tid)];
            const double x = d.x_min + (ix + 0.5) * dpx;
            const double y = d.y_min + (iy + 0.5) * dpy;
            for (int iz = 0; iz < res; ++iz) {
                const double z = d.z_min + (iz + 0.5) * dpz;
                const size_t idx =
                    (static_cast<size_t>(ix) * res + static_cast<size_t>(iy)) * res + static_cast<size_t>(iz);

                double rx, ry, rz;
                int pid;
                if (!con.find_voronoi_cell(x, y, z, rx, ry, rz, pid)) {
                    ++unresolved;
                    continue;
                }
                map.voronoi_cell_id[idx] = pid;

                auto it = mapping_by_id.find(pid);
                if (it == mapping_by_id.end()) {
                    ++unmapped;
                    continue;
                }
                map.density[idx] = it->second->density;
                map.vx[idx] = it->second->vx;
                map.vy[idx] = it->second->vy;
                map.vz[idx] = it->second->vz;
            }
        }
    }

    if (unresolved > 0) {
        std::cerr << "[kiwi] warning: " << unresolved << " of " << n
                  << " Cartesian grid cell(s) found no enclosing Voronoi cell "
                     "(left at density/velocity 0, voronoi_cell_id 0)\n";
    }
    if (unmapped > 0) {
        std::cerr << "[kiwi] warning: " << unmapped << " of " << n
                  << " Cartesian grid cell(s) resolved to a Voronoi cell id with no mapping "
                     "data (left at density/velocity 0)\n";
    }

    save_cartesian_map(map, final_path);
    return map;
}

void save_cartesian_map(const CartesianMap& map, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("could not open output file: " + path);

    f.write(kMagic, sizeof(kMagic));

    int32_t res32 = map.resolution;
    f.write(reinterpret_cast<const char*>(&res32), sizeof(res32));

    double bounds[6] = {map.domain.x_min, map.domain.x_max, map.domain.y_min,
                         map.domain.y_max, map.domain.z_min, map.domain.z_max};
    f.write(reinterpret_cast<const char*>(bounds), sizeof(bounds));

    auto write_doubles = [&](const std::vector<double>& v) {
        f.write(reinterpret_cast<const char*>(v.data()),
                static_cast<std::streamsize>(v.size() * sizeof(double)));
    };
    write_doubles(map.density);
    write_doubles(map.vx);
    write_doubles(map.vy);
    write_doubles(map.vz);

    f.write(reinterpret_cast<const char*>(map.voronoi_cell_id.data()),
            static_cast<std::streamsize>(map.voronoi_cell_id.size() * sizeof(int)));

    if (!f) throw std::runtime_error("failed writing output file (disk full?): " + path);
}

CartesianMap load_cartesian_map(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("could not open cartesian map file: " + path);

    char magic[8];
    f.read(magic, sizeof(magic));
    if (!f || std::memcmp(magic, kMagic, sizeof(magic)) != 0)
        throw std::runtime_error("cartesian map file has an unrecognized header (wrong format "
                                  "or corrupt): " + path);

    int32_t res32 = 0;
    f.read(reinterpret_cast<char*>(&res32), sizeof(res32));
    if (!f || res32 <= 0)
        throw std::runtime_error("cartesian map file has an invalid resolution: " + path);

    double bounds[6];
    f.read(reinterpret_cast<char*>(bounds), sizeof(bounds));
    if (!f) throw std::runtime_error("cartesian map file is truncated (domain bounds): " + path);

    CartesianMap map;
    map.resolution = res32;
    map.domain.x_min = bounds[0];
    map.domain.x_max = bounds[1];
    map.domain.y_min = bounds[2];
    map.domain.y_max = bounds[3];
    map.domain.z_min = bounds[4];
    map.domain.z_max = bounds[5];

    const size_t n = static_cast<size_t>(map.resolution) * static_cast<size_t>(map.resolution) *
                      static_cast<size_t>(map.resolution);

    auto read_doubles = [&](std::vector<double>& v) {
        v.resize(n);
        f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(double)));
    };
    read_doubles(map.density);
    read_doubles(map.vx);
    read_doubles(map.vy);
    read_doubles(map.vz);

    map.voronoi_cell_id.resize(n);
    f.read(reinterpret_cast<char*>(map.voronoi_cell_id.data()),
           static_cast<std::streamsize>(n * sizeof(int)));

    if (!f) throw std::runtime_error("cartesian map file is truncated/corrupt: " + path);
    return map;
}

} // namespace kiwi
