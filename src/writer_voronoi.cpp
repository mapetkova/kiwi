#include "kiwi/writer_voronoi.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "kiwi/mapping.hpp"
#include "kiwi/tessellation.hpp"

namespace kiwi {

VoronoiWriteStats write_voronoi_text(const Config& cfg) {
    Tessellation tess = run_tessellation_stage(cfg);
    MappingResult mapping = run_mapping_stage(cfg);

    std::unordered_map<int, const CellMapping*> mapping_by_id;
    mapping_by_id.reserve(mapping.cells.size() * 2);
    for (const auto& m : mapping.cells) mapping_by_id[m.id] = &m;

    std::vector<const CellGeometry*> ordered;
    ordered.reserve(tess.cells.size());
    for (const auto& c : tess.cells) {
        if (mapping_by_id.count(c.id)) ordered.push_back(&c);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const CellGeometry* a, const CellGeometry* b) { return a->id < b->id; });

    std::unordered_set<int> written_ids;
    written_ids.reserve(ordered.size() * 2);
    for (const auto* c : ordered) written_ids.insert(c->id);

    VoronoiWriteStats stats;
    stats.cells_written = ordered.size();
    stats.missing_cells = tess.cells.size() - ordered.size();

    std::string out_path = cfg.stage_path("grid-voronoi.txt");
    std::ofstream f(out_path);
    if (!f) throw std::runtime_error("could not open output file: " + out_path);
    f << "# id x y z volume density vx vy vz n_neighbors neighbor_ids...\n";
    f << std::setprecision(15);

    for (const CellGeometry* cell : ordered) {
        const CellMapping& m = *mapping_by_id.at(cell->id);

        std::vector<int> neighbors;
        neighbors.reserve(cell->neighbor_ids.size());
        for (int n : cell->neighbor_ids) {
            if (n > 0) {
                if (written_ids.count(n)) {
                    neighbors.push_back(n);
                } else {
                    ++stats.dropped_neighbor_refs;
                }
            } else {
                neighbors.push_back(n);
            }
        }

        f << cell->id << " " << cell->x << " " << cell->y << " " << cell->z << " "
          << cell->volume << " " << m.density << " " << m.vx << " " << m.vy << " " << m.vz
          << " " << neighbors.size();
        for (int nid : neighbors) f << " " << nid;
        f << "\n";
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
