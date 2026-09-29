#include "kiwi/tessellation.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_set>
#include <iomanip>

#include "voro++.hh"
#include "kiwi/point_set.hpp"
#include "kiwi/voro_grid.hpp"

namespace {
bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}
} // namespace

namespace kiwi {

Tessellation tessellate_voronoi(const PointSet& points, const Domain& d, const GridBlocks& blocks) {
    using namespace voro;

    pre_container pcon(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i) {
        pcon.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    }

    int nx, ny, nz;
    guess_scaled_blocks(pcon, blocks.block_scale_factor, nx, ny, nz);

    container con(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, nx, ny, nz, false, false, false, 8);
    pcon.setup(con);

    Tessellation tess;
    tess.cells.reserve(points.size());

    voronoicell_neighbor c;
    std::vector<int> n_list;
    c_loop_all clo(con);
    double x0, y0, z0;
    if (clo.start()) do {
        if (con.compute_cell(c, clo)) {
            clo.pos(x0, y0, z0);
            c.neighbors(n_list);

            CellGeometry cell;
            cell.id = clo.pid();
            cell.x = x0;
            cell.y = y0;
            cell.z = z0;
            cell.volume = c.volume();
            cell.neighbor_ids = n_list;
            tess.cells.push_back(std::move(cell));
        }
    } while (clo.inc());

    return tess;
}

TessellationDiagnostics diagnose(const PointSet& input, const Tessellation& tess) {
    TessellationDiagnostics diag;

    std::unordered_set<int> produced;
    produced.reserve(tess.cells.size() * 2);
    for (const auto& c : tess.cells) {
        produced.insert(c.id);
        if (c.volume <= 0.0) diag.zero_or_negative_volume_ids.push_back(c.id);
    }
    for (int id : input.id) {
        if (produced.find(id) == produced.end()) diag.missing_cell_ids.push_back(id);
    }
    return diag;
}

void save_tessellation_summary(const Tessellation& tess, const std::string& path) {
    std::ofstream f(path);
    f << std::setprecision(15);
    for (const auto& c : tess.cells) {
        f << c.id << " " << c.x << " " << c.y << " " << c.z << " " << c.volume
          << " " << c.neighbor_ids.size();
        for (int n : c.neighbor_ids) f << " " << n;
        f << "\n";
    }
}

Tessellation load_tessellation(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open tessellation file: " + path);

    Tessellation tess;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ls(line);
        CellGeometry cell;
        size_t n_neighbors;
        if (!(ls >> cell.id >> cell.x >> cell.y >> cell.z >> cell.volume >> n_neighbors))
            throw std::runtime_error("tessellation file is malformed: " + path);
        cell.neighbor_ids.resize(n_neighbors);
        for (size_t i = 0; i < n_neighbors; ++i) {
            if (!(ls >> cell.neighbor_ids[i]))
                throw std::runtime_error("tessellation file is malformed (short neighbor list): " + path);
        }
        tess.cells.push_back(std::move(cell));
    }
    if (tess.cells.empty())
        throw std::runtime_error("tessellation file contained no rows: " + path);
    return tess;
}

Tessellation run_tessellation_stage(const Config& cfg) {
    std::string final_path = cfg.stage_path("tessellation.dat");
    if (file_exists(final_path)) {
        return load_tessellation(final_path);
    }

    PointSet points = run_point_set_stage(cfg);
    Tessellation tess = tessellate_voronoi(points, cfg.domain, cfg.grid_blocks);
    save_tessellation_summary(tess, final_path);
    return tess;
}

} // namespace kiwi
