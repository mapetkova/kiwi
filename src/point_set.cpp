#include "kiwi/point_set.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <iostream>
#include <sys/stat.h>
#include <cmath>
#include <algorithm>
#include <array>
#include <unordered_map>

#include "voro++.hh"
#include "kiwi/diagnostics.hpp"
#include "kiwi/snapshot.hpp"
#include "kiwi/voro_grid.hpp"

namespace {
bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}

// Packs three bin-grid indices into one 64-bit key for the spatial hash
// below. Each component is offset by a large bias and clamped to 20 bits
// (+/- ~500k bins per axis), which is far more than any real domain/spacing
// combination in this pipeline needs.
int64_t bin_key(int ix, int iy, int iz) {
    auto biased = [](int v) -> int64_t {
        int64_t b = static_cast<int64_t>(v) + (1 << 19);
        return std::max<int64_t>(0, std::min<int64_t>((1 << 20) - 1, b));
    };
    return (biased(ix) << 40) | (biased(iy) << 20) | biased(iz);
}
} // namespace

namespace kiwi {

void PointSet::append(int pid, double px, double py, double pz) {
    id.push_back(pid);
    x.push_back(px);
    y.push_back(py);
    z.push_back(pz);
}

PointSet load_point_set(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open point-set file: " + path);

    PointSet points;
    int pid;
    double px, py, pz;
    while (f >> pid >> px >> py >> pz) {
        points.append(pid, px, py, pz);
    }
    if (points.size() == 0)
        throw std::runtime_error("point-set file contained no rows (wrong format?): " + path);
    return points;
}

void save_point_set(const PointSet& points, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("could not open point-set file for writing: " + path);
    f << std::setprecision(15);
    for (size_t i = 0; i < points.size(); ++i) {
        f << points.id[i] << " " << points.x[i] << " " << points.y[i] << " " << points.z[i] << "\n";
    }
}

// Exact port of add-buffer-cells.py. Kept deliberately close to the
// original's variable names/structure (rather than "cleaned up") so it can
// be diffed against the source it was ported from.
PointSet add_legacy_brick_buffer(const PointSet& points, const Domain& d,
                                  const LegacyBrickBufferParams& p) {
    PointSet out = points;
    long l = static_cast<long>(points.size());
    long num = 0;

    const double x_min = d.x_min, x_max = d.x_max;
    const double y_min = d.y_min, y_max = d.y_max;
    const double z_min = d.z_min, z_max = d.z_max;
    const double dx = p.dx, dy = p.dy, dz = p.dz;

    // Block 1: two shells on the y_min/y_max faces (spacing_1), spanning a
    // z sub-range narrowed by (dx - dz) / 2 on each side.
    {
        double a = p.shell_spacing_1;
        double y = y_min + a / 2.0;
        double x = x_min + a / 2.0;
        while (x < x_max) {
            double z = (z_min + (dx - dz) / 2.0) + a / 2.0;
            while (z < z_max - (dx - dz) / 2.0) {
                num++;
                out.append(static_cast<int>(num + l), x, y, z);
                num++;
                out.append(static_cast<int>(num + l), x, y_max - a / 2.0, z);
                z += a;
            }
            x += a;
        }
    }

    // Block 2: first shell on the z_min/z_max faces (spacing_1 again),
    // spanning the full y range.
    double z_hi, z_lo;
    {
        double a = p.shell_spacing_1;
        z_hi = z_max - (dx - dz) / 2.0 + a / 2.0;
        z_lo = z_min + (dx - dz) / 2.0 - a / 2.0;
        double x = x_min + a / 2.0;
        while (x < x_max) {
            double y = y_min + a / 2.0;
            while (y < y_max) {
                num++;
                out.append(static_cast<int>(num + l), x, y, z_hi);
                num++;
                out.append(static_cast<int>(num + l), x, y, z_lo);
                y += a;
            }
            x += a;
        }
    }

    // Block 3: second z-face shell, stepped outward (spacing_2).
    {
        double a = p.shell_spacing_2;
        z_hi = z_hi + (a - 1.0) / 2.0 + a / 2.0;
        z_lo = z_lo - (a - 1.0) / 2.0 - a / 2.0;
        double x = x_min + a / 2.0;
        while (x < x_max) {
            double y = y_min + a / 2.0;
            while (y < y_max) {
                num++;
                out.append(static_cast<int>(num + l), x, y, z_hi);
                num++;
                out.append(static_cast<int>(num + l), x, y, z_lo);
                y += a;
            }
            x += a;
        }
    }

    // Block 4: third z-face shell, stepped outward again (spacing_3).
    {
        double a = p.shell_spacing_3;
        z_hi = z_hi + (a - 1.0) / 2.0 + a / 2.0;
        z_lo = z_lo - (a - 1.0) / 2.0 - a / 2.0;
        double x = x_min + a / 2.0;
        while (x < x_max) {
            double y = y_min + a / 2.0;
            while (y < y_max) {
                num++;
                out.append(static_cast<int>(num + l), x, y, z_hi);
                num++;
                out.append(static_cast<int>(num + l), x, y, z_lo);
                y += a;
            }
            x += a;
        }
    }

    return out;
}

// See point_set.hpp for the full algorithm description.
CartesianLattice generate_cartesian_lattice(const Domain& d, double spacing_fraction, int id_offset) {
    const double extent_x = d.x_max - d.x_min;
    const double extent_y = d.y_max - d.y_min;
    const double extent_z = d.z_max - d.z_min;
    const double spacing = spacing_fraction * std::min({extent_x, extent_y, extent_z});

    // Cell-centered: first candidate is spacing/2 in from x_min, last is
    // spacing/2 in from x_max, so no candidate sits exactly on a domain
    // wall.
    CartesianLattice lat;
    lat.spacing = spacing;
    lat.nx = std::max(1, static_cast<int>(std::round(extent_x / spacing)));
    lat.ny = std::max(1, static_cast<int>(std::round(extent_y / spacing)));
    lat.nz = std::max(1, static_cast<int>(std::round(extent_z / spacing)));

    long num = 0;
    for (int i = 0; i < lat.nx; ++i) {
        double x = d.x_min + spacing * (i + 0.5);
        for (int j = 0; j < lat.ny; ++j) {
            double y = d.y_min + spacing * (j + 0.5);
            for (int k = 0; k < lat.nz; ++k) {
                double z = d.z_min + spacing * (k + 0.5);
                ++num;
                lat.points.append(static_cast<int>(num + id_offset), x, y, z);
            }
        }
    }

    return lat;
}

// Generic particle-coverage buffer. See CoverageGridBufferParams (config.hpp)
// and point_set.hpp for the algorithm description.
//
// Coverage test uses a uniform spatial hash over the real particles, cell
// size = coverage_radius_factor * h_max (the largest smoothing length in the
// input), so that checking a candidate's own bin plus its 3x3x3 neighborhood
// is guaranteed to find every particle whose coverage sphere could possibly
// reach it -- no false negatives. Known trade-off, deliberately deferred: a
// single large-h outlier coarsens the hash for every candidate, not just the
// ones near that outlier; revisit if this shows up as a real bottleneck.
PointSet add_coverage_grid_buffer(const PointSet& points, const std::vector<double>& h,
                                   const Domain& d, const CoverageGridBufferParams& params) {
    if (h.size() != points.size())
        throw std::runtime_error("add_coverage_grid_buffer: h.size() must match points.size()");

    double h_max = 0.0;
    for (double hi : h) h_max = std::max(h_max, hi);
    if (h_max <= 0.0)
        throw std::runtime_error("add_coverage_grid_buffer: all smoothing lengths are <= 0");

    const double bin_size = params.coverage_radius_factor * h_max;

    auto bin_of = [&](double x, double y, double z) {
        return std::array<int, 3>{
            static_cast<int>(std::floor(x / bin_size)),
            static_cast<int>(std::floor(y / bin_size)),
            static_cast<int>(std::floor(z / bin_size))};
    };

    std::unordered_map<int64_t, std::vector<size_t>> hash;
    hash.reserve(points.size() * 2);
    for (size_t i = 0; i < points.size(); ++i) {
        auto b = bin_of(points.x[i], points.y[i], points.z[i]);
        hash[bin_key(b[0], b[1], b[2])].push_back(i);
    }

    auto is_covered = [&](double x, double y, double z) {
        auto b = bin_of(x, y, z);
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = hash.find(bin_key(b[0] + dx, b[1] + dy, b[2] + dz));
                    if (it == hash.end()) continue;
                    for (size_t i : it->second) {
                        double ddx = points.x[i] - x, ddy = points.y[i] - y, ddz = points.z[i] - z;
                        double r = params.coverage_radius_factor * h[i];
                        if (ddx * ddx + ddy * ddy + ddz * ddz <= r * r) return true;
                    }
                }
            }
        }
        return false;
    };

    // Candidate grid generation itself is shared with the lattice-based
    // Voronoi regression test (tools/test_lattice_grid.cpp) -- see
    // generate_cartesian_lattice() above. Ids get renumbered here to just the
    // surviving (uncovered) candidates, sequential starting at
    // points.size()+1, matching add_legacy_brick_buffer's convention -- the
    // lattice's own ids (which number every candidate, covered or not) are
    // discarded.
    CartesianLattice lat = generate_cartesian_lattice(d, params.spacing_fraction);

    PointSet out = points;
    long l = static_cast<long>(points.size());
    long num = 0;
    for (size_t idx = 0; idx < lat.points.size(); ++idx) {
        double x = lat.points.x[idx], y = lat.points.y[idx], z = lat.points.z[idx];
        if (!is_covered(x, y, z)) {
            ++num;
            out.append(static_cast<int>(num + l), x, y, z);
        }
    }

    return out;
}

PointSet lloyd_step(const PointSet& points, const Domain& d, const GridBlocks& blocks) {
    using namespace voro;

    pre_container pcon(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i) {
        pcon.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    }

    int nx, ny, nz;
    guess_scaled_blocks(pcon, blocks.block_scale_factor, nx, ny, nz);

    container con(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, nx, ny, nz, false, false, false, 8);
    pcon.setup(con);

    PointSet out;
    voronoicell_neighbor c;
    c_loop_all clo(con);
    double x0, y0, z0, cx, cy, cz;
    if (clo.start()) do {
        if (con.compute_cell(c, clo)) {
            clo.pos(x0, y0, z0);
            c.centroid(cx, cy, cz);
            out.append(clo.pid(), x0 + cx, y0 + cy, z0 + cz);
        }
    } while (clo.inc());

    return out;
}

PointSet run_point_set_stage(const Config& cfg) {
    std::string final_path = cfg.stage_path("points-final.dat");
    if (file_exists(final_path)) {
        return load_point_set(final_path);
    }

    SPHSnapshot snap = load_sph_snapshot(cfg.input_snapshot, cfg.input_snapshot_format, cfg.kernel_hfact(), cfg.phantom_units,
                                          cfg.allow_missing_velocity);
    PointSet points = extract_points(snap);
    save_point_set(points, cfg.stage_path("points-raw.dat"));

    // Real-particle count, fixed at this point regardless of whether a
    // buffer runs below -- every buffer strategy appends new points with ids
    // starting at n_real+1 (see add_legacy_brick_buffer/
    // add_coverage_grid_buffer), so "id <= n_real" identifies a real
    // (SPH-derived) site at any later stage, including after Lloyd
    // relaxation reorders rows. Diagnostics uses this split; nothing else
    // currently needs it.
    const size_t n_real = points.size();

    if (cfg.buffer.enabled) {
        size_t n_before = points.size();
        if (cfg.buffer.strategy == "legacy_brick") {
            points = add_legacy_brick_buffer(points, cfg.domain, cfg.buffer.legacy_brick);
        } else if (cfg.buffer.strategy == "coverage_grid") {
            // snap.h is aligned by index with the points extract_points()
            // just produced from it -- must run before Lloyd relaxation
            // moves sites away from their originating particles.
            points = add_coverage_grid_buffer(points, snap.h, cfg.domain, cfg.buffer.coverage_grid);
        } else {
            throw std::runtime_error(
                "buffer.strategy '" + cfg.buffer.strategy +
                "' is not recognized (expected 'legacy_brick' or 'coverage_grid')");
        }
        std::cerr << "[kiwi] buffer (" << cfg.buffer.strategy << "): added "
                  << (points.size() - n_before) << " buffer point(s), " << points.size()
                  << " total\n";
        std::string buff_path = cfg.stage_path("points-buff.dat");
        save_point_set(points, buff_path);

        if (cfg.diagnostics.enabled) {
            std::ostringstream title;
            title << "kiwi diagnostic -- " << cfg.run_name << "\nafter buffer ("
                  << cfg.buffer.strategy << "): " << n_real << " real, "
                  << (points.size() - n_real) << " buffer";
            generate_diagnostic_plot(cfg, buff_path, n_real, title.str(), "buffer_diagnostic.png");
        }
    }

    for (int it = 1; it <= cfg.lloyd.iterations; ++it) {
        points = lloyd_step(points, cfg.domain, cfg.grid_blocks);
        if (cfg.lloyd.save_each_iteration) {
            save_point_set(points, cfg.stage_path("points-buff-iter-" + std::to_string(it) + ".dat"));
        }
    }

    save_point_set(points, final_path);

    if (cfg.diagnostics.enabled && cfg.lloyd.iterations > 0) {
        std::ostringstream title;
        title << "kiwi diagnostic -- " << cfg.run_name << "\nafter Lloyd relaxation ("
              << cfg.lloyd.iterations << " iter): " << n_real << " real, "
              << (points.size() - n_real) << " buffer";
        generate_diagnostic_plot(cfg, final_path, n_real, title.str(), "lloyd_diagnostic.png");
    }

    return points;
}

} // namespace kiwi
