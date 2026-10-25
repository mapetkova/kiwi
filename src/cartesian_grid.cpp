#include "kiwi/cartesian_grid.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "kiwi/kernel_integral.hpp"
#include "kiwi/kernel_integrator.hpp"
#include "kiwi/snapshot.hpp"

namespace {
bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}
} // namespace

namespace kiwi {

namespace {

constexpr char kMagic[8] = {'V', 'G', 'T', 'C', 'G', 'R', 'D', '2'};

// Builds the 6 axis-aligned rectangular faces (each 4 vertices, flattened
// as x,y,z triples) of a box cell centered at (cx,cy,cz) with half-extents
// (hx,hy,hz), wound to match compute_face_integral()'s sign convention --
// see kernel_integral.hpp's doc comment on that function for what this
// winding empirically needs to be and how it was validated (NOT the naive
// "(v1-v0) x (v2-v0) points outward" reading of voro++.hh's own documented
// convention -- that gives the wrong sign here; see the note below).
//
// This is deliberately fixed geometry, not something to hand-tune per call
// site: getting a face's winding backward silently flips that face's sign
// contribution to Msum rather than erroring. Validated two independent ways
// before settling on this winding: (1) a synthetic cube face-integral test
// (a large box fully enclosing an h=1 particle's kernel support must sum to
// Msum == 1.0, the kernel's own normalization) and (2) an end-to-end ground
// truth -- a single SPH particle centered in a cubic domain has exactly one
// Voronoi cell (the whole domain box), so running the `kiwi map` pipeline on
// that exact setup gives an independent density*volume == 1.0 to compare
// against. The naive outward-normal
// winding gave Msum == -1.0 on test (1) (wrong sign, right magnitude) until
// reversed to match test (2)'s independently-produced ground truth.
void box_cell_faces(double cx, double cy, double cz, double hx, double hy, double hz,
                     std::vector<double> faces[6]) {
    const double xm = cx - hx, xp = cx + hx;
    const double ym = cy - hy, yp = cy + hy;
    const double zm = cz - hz, zp = cz + hz;

    // +x face (in-plane axes u=y, v=z): (u-,v+),(u+,v+),(u+,v-),(u-,v-)
    faces[0] = {xp, ym, zp, xp, yp, zp, xp, yp, zm, xp, ym, zm};
    // -x face
    faces[1] = {xm, yp, zm, xm, yp, zp, xm, ym, zp, xm, ym, zm};
    // +y face (u=z, v=x): (u-,v-),(u+,v-),(u+,v+),(u-,v+) -- see note above
    faces[2] = {xp, yp, zm, xp, yp, zp, xm, yp, zp, xm, yp, zm};
    // -y face
    faces[3] = {xm, ym, zp, xp, ym, zp, xp, ym, zm, xm, ym, zm};
    // +z face
    faces[4] = {xm, yp, zp, xp, yp, zp, xp, ym, zp, xm, ym, zp};
    // -z face
    faces[5] = {xp, ym, zm, xp, yp, zm, xm, yp, zm, xm, ym, zm};
}

} // namespace

CartesianGridResult run_sph_to_grid_stage(const Config& cfg) {
    std::string final_path = cfg.stage_path("cartesian_grid.dat");
    if (file_exists(final_path)) {
        return load_cartesian_grid(final_path);
    }

    const auto& d = cfg.domain;

    // Resolve per-axis cell counts: an explicit resolution_x/y/z overrides
    // the shared `resolution` for that axis only (0 = unset, per
    // CartesianGridConfig's doc comment) -- lets a non-cube domain be given
    // deliberately-chosen, independently-sized axes instead of forcing one
    // shared count everywhere.
    const int nx = cfg.cartesian_grid.resolution_x > 0 ? cfg.cartesian_grid.resolution_x : cfg.cartesian_grid.resolution;
    const int ny = cfg.cartesian_grid.resolution_y > 0 ? cfg.cartesian_grid.resolution_y : cfg.cartesian_grid.resolution;
    const int nz = cfg.cartesian_grid.resolution_z > 0 ? cfg.cartesian_grid.resolution_z : cfg.cartesian_grid.resolution;
    if (nx <= 0 || ny <= 0 || nz <= 0)
        throw std::runtime_error("cartesian_grid resolution (resolved nx/ny/nz) must be > 0");

    const double dpx = (d.x_max - d.x_min) / nx;
    const double dpy = (d.y_max - d.y_min) / ny;
    const double dpz = (d.z_max - d.z_min) / nz;
    {
        // Checks the RESOLVED cell dimensions, not the domain's own shape --
        // per-axis resolution can make cells cubes even on a non-cube
        // domain (that's the point of resolution_x/y/z), so warning on the
        // domain shape alone would be a false positive in that case.
        const double max_extent = std::max({dpx, dpy, dpz});
        const double tol = 1e-6 * max_extent;
        if (std::abs(dpx - dpy) > tol || std::abs(dpx - dpz) > tol || std::abs(dpy - dpz) > tol) {
            std::cerr << "[kiwi] warning: resolved Cartesian grid cells are not cubes (dx=" << dpx
                      << ", dy=" << dpy << ", dz=" << dpz << ", from nx=" << nx << " ny=" << ny
                      << " nz=" << nz << " over domain extents " << (d.x_max - d.x_min) << " x "
                      << (d.y_max - d.y_min) << " x " << (d.z_max - d.z_min) << ") -- density/velocity "
                         "values are still correct, only the implied cell shape is skewed. Set "
                         "cartesian_grid.resolution_x/y/z to match the domain's own aspect ratio if cube "
                         "cells matter for your use case.\n";
        }
    }

    DensityLookupTable table;
    if (cfg.cartesian_grid.use_tabulated_integral) {
        // Same load-or-generate-in-memory fallback as run_mapping_stage()
        // (see MappingConfig's doc comment in config.hpp) -- these are two
        // independent config knob sets by design, but this one behavior is
        // worth keeping identical between them.
        if (!cfg.cartesian_grid.density_table.empty() && file_exists(cfg.cartesian_grid.density_table)) {
            std::cerr << "[kiwi] loading density lookup table from " << cfg.cartesian_grid.density_table
                      << "...\n";
            table = load_density_table(cfg.cartesian_grid.density_table);
        } else {
            if (cfg.cartesian_grid.density_table.empty()) {
                std::cerr << "[kiwi] cartesian_grid.density_table not set -- generating the "
                             "kernel-integral lookup table in memory instead of loading one from disk...\n";
            } else {
                std::cerr << "[kiwi] warning: cartesian_grid.density_table (" << cfg.cartesian_grid.density_table
                          << ") not found on disk -- generating it in memory instead...\n";
            }
            table = generate_density_table();
        }
        std::cerr << "[kiwi] density lookup table ready (nr0=" << table.nr0 << " nR_0=" << table.nR_0
                   << " nphi=" << table.nphi << ")\n";
    } else {
        std::cerr << "[kiwi] using full_integral_new (analytic) -- no lookup table needed\n";
    }

    const KernelIntegrator kernel = make_kernel_integrator(cfg, cfg.cartesian_grid.use_tabulated_integral, &table);

    SPHSnapshot snap = load_sph_snapshot(cfg.input_snapshot, cfg.input_snapshot_format, cfg.kernel_hfact(), cfg.phantom_units,
                                          cfg.allow_missing_velocity);
    std::cerr << "[kiwi] sph-to-grid: loaded " << snap.size() << " SPH particles (no point-set/"
                 "tessellation/mapping stage involved -- this path reads the snapshot directly)\n";

    CartesianGridResult grid;
    grid.nx = nx;
    grid.ny = ny;
    grid.nz = nz;
    grid.domain = d;
    const size_t n = static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
    grid.density.assign(n, 0.0);
    grid.vx.assign(n, 0.0);
    grid.vy.assign(n, 0.0);
    grid.vz.assign(n, 0.0);

    const double hpx = 0.5 * dpx, hpy = 0.5 * dpy, hpz = 0.5 * dpz;
    const double cell_volume = dpx * dpy * dpz;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    std::cerr << "[kiwi] sph-to-grid: running with " << n_threads << " thread(s)\n";

    // Per-thread accumulator grids, reduced after the parallel region --
    // same reasoning as run_mapping_stage()'s dens_cells_t/v{x,y,z}_cells_t
    // (mapping.cpp): nearby particles routinely deposit into the same
    // handful of cells, so a shared accumulator would need locking that
    // erases most of the parallel speedup. Cost here: O(resolution^3 *
    // n_threads) extra memory -- e.g. res=128, 8 threads is ~512MB (4
    // arrays x 8 threads x 2,097,152 cells x 8 bytes). Trivial next to a
    // real SPH snapshot's own memory footprint at usual
    // resolutions, but worth knowing before cranking resolution and thread
    // count both up at once.
    std::vector<std::vector<double>> dens_t(n_threads, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> vx_t(n_threads, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> vy_t(n_threads, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> vz_t(n_threads, std::vector<double>(n, 0.0));

    double vx_mean = 0, vy_mean = 0, vz_mean = 0;
    long skipped_outside = 0; // particle's kernel support doesn't overlap the grid at all
    const long n_particles = static_cast<long>(snap.size());

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256) reduction(+:vx_mean,vy_mean,vz_mean,skipped_outside)
#endif
    for (long idx = 0; idx < n_particles; ++idx) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<double>& dens_cells = dens_t[tid];
        std::vector<double>& vx_cells = vx_t[tid];
        std::vector<double>& vy_cells = vy_t[tid];
        std::vector<double>& vz_cells = vz_t[tid];

        const double x0 = snap.x[idx], y0 = snap.y[idx], z0 = snap.z[idx];
        const double vx = snap.vx[idx], vy = snap.vy[idx], vz = snap.vz[idx];
        const double m = snap.m[idx], h = snap.h[idx];

        vx_mean += vx;
        vy_mean += vy;
        vz_mean += vz;

        if (h <= 0.0) continue; // degenerate/invalid smoothing length -- no contribution

        const double radius = kernel.support(h); // kernel support radius

        // Arithmetic cell-index range this particle's kernel could reach --
        // no discovery/search at all, unlike run_mapping_stage()'s BFS (see
        // this file's header comment). A bounding BOX around the kernel's
        // support sphere, not trimmed to the sphere itself -- same as
        // SPLASH's interpolate3D.f90: cells in the box's corners the sphere
        // doesn't actually reach just evaluate to ~0 contribution below
        // rather than being explicitly excluded up front.
        int ixmin = static_cast<int>(std::floor((x0 - radius - d.x_min) / dpx));
        int ixmax = static_cast<int>(std::floor((x0 + radius - d.x_min) / dpx));
        int iymin = static_cast<int>(std::floor((y0 - radius - d.y_min) / dpy));
        int iymax = static_cast<int>(std::floor((y0 + radius - d.y_min) / dpy));
        int izmin = static_cast<int>(std::floor((z0 - radius - d.z_min) / dpz));
        int izmax = static_cast<int>(std::floor((z0 + radius - d.z_min) / dpz));

        ixmin = std::max(ixmin, 0);
        ixmax = std::min(ixmax, nx - 1);
        iymin = std::max(iymin, 0);
        iymax = std::min(iymax, ny - 1);
        izmin = std::max(izmin, 0);
        izmax = std::min(izmax, nz - 1);

        if (ixmin > ixmax || iymin > iymax || izmin > izmax) {
            ++skipped_outside; // kernel support doesn't overlap the grid at all
            continue;
        }

        for (int ix = ixmin; ix <= ixmax; ++ix) {
            const double cx = d.x_min + (ix + 0.5) * dpx;
            for (int iy = iymin; iy <= iymax; ++iy) {
                const double cy = d.y_min + (iy + 0.5) * dpy;
                for (int iz = izmin; iz <= izmax; ++iz) {
                    const double cz = d.z_min + (iz + 0.5) * dpz;

                    std::vector<double> faces[6];
                    box_cell_faces(cx, cy, cz, hpx, hpy, hpz, faces);

                    // Sum the kernel-overlap integral across this cell's 6
                    // faces -- exactly the same "closed cell boundary flux
                    // sums to the volume overlap integral" principle
                    // run_mapping_stage()'s compute_cell_integral() relies
                    // on for an arbitrary Voronoi cell, just applied here to
                    // a box's fixed 6 faces instead of walking however many
                    // faces a Voronoi cell happens to have.
                    double Msum = 0.0;
                    for (int f = 0; f < 6; ++f) {
                        Msum += kernel.face(x0, y0, z0, h, faces[f]);
                    }

                    // Same clamp as compute_cell_integral() applies to its
                    // own Msum -- a small negative value is a face-marching
                    // numerical artefact, not a real negative mass.
                    if (Msum < 0.0) {
                        if (Msum <= -0.01)
                            std::cerr << "[kiwi] warning: negative mass fraction " << Msum << " (particle "
                                      << snap.id[idx] << ", cell (" << ix << "," << iy << "," << iz << "))\n";
                        Msum = 0.0;
                    }
                    if (Msum == 0.0) continue; // common for the bounding box's non-overlapping corner cells

                    const size_t idx3 =
                        (static_cast<size_t>(ix) * ny + static_cast<size_t>(iy)) * nz + static_cast<size_t>(iz);
                    dens_cells[idx3] += Msum * m;
                    vx_cells[idx3] += Msum * m * vx;
                    vy_cells[idx3] += Msum * m * vy;
                    vz_cells[idx3] += Msum * m * vz;
                }
            }
        }
    }

    std::cerr << "[kiwi] sph-to-grid: reducing per-thread accumulators...\n";
    std::vector<double> dens_cells(n, 0.0), vx_cells(n, 0.0), vy_cells(n, 0.0), vz_cells(n, 0.0);
    for (int t = 0; t < n_threads; ++t) {
        for (size_t i = 0; i < n; ++i) {
            dens_cells[i] += dens_t[t][i];
            vx_cells[i] += vx_t[t][i];
            vy_cells[i] += vy_t[t][i];
            vz_cells[i] += vz_t[t][i];
        }
    }

    if (snap.size() > 0) {
        vx_mean /= static_cast<double>(snap.size());
        vy_mean /= static_cast<double>(snap.size());
        vz_mean /= static_cast<double>(snap.size());
    }

    std::cerr << "[kiwi] sph-to-grid: " << skipped_outside
              << " particle(s)' kernel support did not overlap the grid at all\n";

    for (size_t i = 0; i < n; ++i) {
        grid.density[i] = dens_cells[i] / cell_volume;
        if (dens_cells[i] > 0.0) {
            grid.vx[i] = vx_cells[i] / dens_cells[i];
            grid.vy[i] = vy_cells[i] / dens_cells[i];
            grid.vz[i] = vz_cells[i] / dens_cells[i];
        } else {
            // No particle's kernel reached this cell at all -- fall back to
            // the mean velocity across all particles, same convention
            // run_mapping_stage() uses for its own unreached cells.
            grid.vx[i] = vx_mean;
            grid.vy[i] = vy_mean;
            grid.vz[i] = vz_mean;
        }
    }

    save_cartesian_grid(grid, final_path);
    return grid;
}

void save_cartesian_grid(const CartesianGridResult& grid, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("could not open output file: " + path);

    f.write(kMagic, sizeof(kMagic));

    int32_t nx32 = grid.nx, ny32 = grid.ny, nz32 = grid.nz;
    f.write(reinterpret_cast<const char*>(&nx32), sizeof(nx32));
    f.write(reinterpret_cast<const char*>(&ny32), sizeof(ny32));
    f.write(reinterpret_cast<const char*>(&nz32), sizeof(nz32));

    double bounds[6] = {grid.domain.x_min, grid.domain.x_max, grid.domain.y_min,
                         grid.domain.y_max, grid.domain.z_min, grid.domain.z_max};
    f.write(reinterpret_cast<const char*>(bounds), sizeof(bounds));

    auto write_doubles = [&](const std::vector<double>& v) {
        f.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(double)));
    };
    write_doubles(grid.density);
    write_doubles(grid.vx);
    write_doubles(grid.vy);
    write_doubles(grid.vz);

    if (!f) throw std::runtime_error("failed writing output file (disk full?): " + path);
}

CartesianGridResult load_cartesian_grid(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("could not open cartesian grid file: " + path);

    char magic[8];
    f.read(magic, sizeof(magic));
    if (!f || std::memcmp(magic, kMagic, sizeof(magic)) != 0)
        throw std::runtime_error("cartesian grid file has an unrecognized header (wrong format or "
                                  "corrupt): " + path);

    int32_t nx32 = 0, ny32 = 0, nz32 = 0;
    f.read(reinterpret_cast<char*>(&nx32), sizeof(nx32));
    f.read(reinterpret_cast<char*>(&ny32), sizeof(ny32));
    f.read(reinterpret_cast<char*>(&nz32), sizeof(nz32));
    if (!f || nx32 <= 0 || ny32 <= 0 || nz32 <= 0)
        throw std::runtime_error("cartesian grid file has an invalid resolution (nx/ny/nz): " + path);

    double bounds[6];
    f.read(reinterpret_cast<char*>(bounds), sizeof(bounds));
    if (!f) throw std::runtime_error("cartesian grid file is truncated (domain bounds): " + path);

    CartesianGridResult grid;
    grid.nx = nx32;
    grid.ny = ny32;
    grid.nz = nz32;
    grid.domain.x_min = bounds[0];
    grid.domain.x_max = bounds[1];
    grid.domain.y_min = bounds[2];
    grid.domain.y_max = bounds[3];
    grid.domain.z_min = bounds[4];
    grid.domain.z_max = bounds[5];

    const size_t n = static_cast<size_t>(grid.nx) * static_cast<size_t>(grid.ny) * static_cast<size_t>(grid.nz);

    auto read_doubles = [&](std::vector<double>& v) {
        v.resize(n);
        f.read(reinterpret_cast<char*>(v.data()), static_cast<std::streamsize>(n * sizeof(double)));
    };
    read_doubles(grid.density);
    read_doubles(grid.vx);
    read_doubles(grid.vy);
    read_doubles(grid.vz);

    if (!f) throw std::runtime_error("cartesian grid file is truncated/corrupt: " + path);
    return grid;
}

} // namespace kiwi
