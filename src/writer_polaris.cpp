#include "kiwi/writer_polaris.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "kiwi/mapping.hpp"
#include "kiwi/tessellation.hpp"

namespace kiwi {

namespace {

// POLARIS GRID data-quantity IDs, fixed order: [gas temperature, vx, vy, vz,
// gas density]. Not derived from anything -- these are the values the
// target POLARIS setup expects, treated as opaque POLARIS-internal
// constants. See
// config.hpp's PolarisWriterConfig comment for why these aren't config yet.
constexpr unsigned short kParams = 5;
constexpr unsigned short kPlist[kParams] = {3, 7, 8, 9, 28};

template <typename T>
void write_raw(std::ofstream& f, const T& v) {
    f.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

} // namespace

PolarisWriteStats write_polaris_grid(const Config& cfg) {
    Tessellation tess = run_tessellation_stage(cfg);
    MappingResult mapping = run_mapping_stage(cfg);

    std::unordered_map<int, const CellMapping*> mapping_by_id;
    mapping_by_id.reserve(mapping.cells.size() * 2);
    for (const auto& m : mapping.cells) mapping_by_id[m.id] = &m;

    // Cells actually written: every tessellated cell that also has mapping
    // data, in ascending id order (i.e. the order of the original particle
    // ids), so rows can be compared one-to-one with a grid written from the
    // same particles.
    std::vector<const CellGeometry*> ordered;
    ordered.reserve(tess.cells.size());
    for (const auto& c : tess.cells) {
        if (mapping_by_id.count(c.id)) ordered.push_back(&c);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const CellGeometry* a, const CellGeometry* b) { return a->id < b->id; });

    PolarisWriteStats stats;
    stats.cells_written = ordered.size();
    stats.missing_cells = tess.cells.size() - ordered.size();

    // Explicit id -> output-row-index map, built from the cells we're
    // actually about to write. Assuming row = id - 1 would only be valid when
    // every id in [1, max_id] produced a cell (no gaps); this mapping stays
    // correct even when some cells are missing.
    std::unordered_map<int, int> row_of_id;
    row_of_id.reserve(ordered.size() * 2);
    for (size_t row = 0; row < ordered.size(); ++row) {
        row_of_id[ordered[row]->id] = static_cast<int>(row);
    }

    const auto& d = cfg.domain;
    const auto& pw = cfg.polaris_writer;
    const double midx = (d.x_min + d.x_max) / 2.0;
    const double midy = (d.y_min + d.y_max) / 2.0;
    const double midz = (d.z_min + d.z_max) / 2.0;

    // POLARIS's own Voronoi grid loader (GridVoronoi::loadGridFromBinaryFile)
    // stores a single scalar max_len in the header and treats the domain as
    // a CUBE spanning [-0.5*max_len, +0.5*max_len] on all three axes -- not
    // a convention of ours, baked into POLARIS itself. `len` below is
    // derived from the x-extent alone, which
    // is only correct if the y/z extents equal it too. A non-cubic domain
    // here would silently write a `len` that doesn't match the actual y/z
    // extent, producing a wrong grid with no error -- warn, don't block,
    // same pattern as the OpenMP-not-found/manifest-mismatch cases elsewhere.
    {
        const double dx = d.x_max - d.x_min;
        const double dy = d.y_max - d.y_min;
        const double dz = d.z_max - d.z_min;
        const double max_extent = std::max({dx, dy, dz});
        const double tol = 1e-6 * max_extent; // relative tolerance, scale-invariant
        if (std::abs(dx - dy) > tol || std::abs(dx - dz) > tol || std::abs(dy - dz) > tol) {
            std::cerr << "[kiwi] warning: domain is not a cube (x extent=" << dx
                      << ", y extent=" << dy << ", z extent=" << dz
                      << ") -- POLARIS's Voronoi grid format assumes a single cubic side length "
                         "(derived here from the x extent only); the written grid's y/z scale "
                         "will be wrong. Adjust domain so x_max-x_min == y_max-y_min == z_max-z_min.\n";
        }
    }

    std::string out_path = cfg.stage_path("grid-polaris.dat");
    std::ofstream f(out_path, std::ios::binary);
    if (!f) throw std::runtime_error("could not open output file: " + out_path);

    unsigned short sim = pw.sim_id;
    write_raw(f, sim);
    write_raw(f, kParams);
    f.write(reinterpret_cast<const char*>(kPlist), sizeof(kPlist));

    double cells_d = static_cast<double>(ordered.size());
    double len = (d.x_max - d.x_min) * pw.length_scale_m;
    write_raw(f, cells_d);
    write_raw(f, len);

    for (const CellGeometry* cell : ordered) {
        const CellMapping& m = *mapping_by_id.at(cell->id);

        float x0 = static_cast<float>((cell->x - midx) * pw.length_scale_m);
        float y0 = static_cast<float>((cell->y - midy) * pw.length_scale_m);
        float z0 = static_cast<float>((cell->z - midz) * pw.length_scale_m);
        write_raw(f, x0);
        write_raw(f, y0);
        write_raw(f, z0);

        double vol = cell->volume * pw.length_scale_m * pw.length_scale_m * pw.length_scale_m;
        write_raw(f, vol);

        float Tg = static_cast<float>(pw.gas_temperature_k);
        float vx0 = static_cast<float>(m.vx * pw.velocity_scale);
        float vy0 = static_cast<float>(m.vy * pw.velocity_scale);
        float vz0 = static_cast<float>(m.vz * pw.velocity_scale);
        float dens_f = static_cast<float>(m.density * pw.density_scale);
        write_raw(f, Tg);
        write_raw(f, vx0);
        write_raw(f, vy0);
        write_raw(f, vz0);
        write_raw(f, dens_f);

        std::vector<int> translated;
        translated.reserve(cell->neighbor_ids.size());
        for (int n : cell->neighbor_ids) {
            if (n > 0) {
                auto it = row_of_id.find(n);
                if (it == row_of_id.end()) {
                    ++stats.dropped_neighbor_refs;
                    continue;
                }
                translated.push_back(it->second);
            } else {
                translated.push_back(n); // voro++ domain-wall marker, unchanged
            }
        }

        int neigh = static_cast<int>(translated.size());
        write_raw(f, neigh);
        for (int nid : translated) write_raw(f, nid);
    }

    f.close();

    std::cerr << "[kiwi] wrote " << stats.cells_written << " cell(s) to " << out_path << "\n";
    if (stats.missing_cells > 0) {
        std::cerr << "[kiwi] warning: " << stats.missing_cells
                  << " tessellated cell(s) had no mapping data and were skipped\n";
    }
    if (stats.dropped_neighbor_refs > 0) {
        std::cerr << "[kiwi] warning: " << stats.dropped_neighbor_refs
                  << " neighbor reference(s) pointed at a missing cell and were dropped\n";
    }

    return stats;
}

} // namespace kiwi
