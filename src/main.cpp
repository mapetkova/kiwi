#include <iostream>

#include "kiwi/cartesian_grid.hpp"
#include "kiwi/cartesian_map.hpp"
#include "kiwi/config.hpp"
#include "kiwi/manifest.hpp"
#include "kiwi/mapping.hpp"
#include "kiwi/point_set.hpp"
#include "kiwi/tessellation.hpp"
#include "kiwi/writer_polaris.hpp"
#include "kiwi/writer_voronoi.hpp"

namespace {

int run_point_set(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);
    kiwi::PointSet points = kiwi::run_point_set_stage(cfg);
    std::cerr << "[kiwi] final point set size: " << points.size() << "\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_tessellate(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    // Chains from the point-set stage. run_point_set_stage() loads its own
    // cached "points-final.dat" instead of recomputing if it already exists
    // under this run's directory (see point_set.hpp) -- so re-running
    // `tessellate` after a `point-set` run doesn't redo buffer/Lloyd from
    // scratch.
    kiwi::PointSet points = kiwi::run_point_set_stage(cfg);
    std::cerr << "[kiwi] loaded " << points.size() << " points\n";

    // Cached like every other stage: run_tessellation_stage() loads
    // "tessellation.dat" directly if it already exists instead of rebuilding
    // the voro++ container from scratch -- this used to be the one stage
    // command that always recomputed even when nothing had changed (and the
    // write-polaris/write-voronoi commands each did their own separate
    // from-scratch tessellation pass on top of that). Delete
    // tessellation.dat (or the run directory) to force a recompute.
    kiwi::Tessellation tess = kiwi::run_tessellation_stage(cfg);
    std::cerr << "[kiwi] computed " << tess.cells.size() << " cells\n";

    auto diag = kiwi::diagnose(points, tess);
    if (!diag.missing_cell_ids.empty()) {
        std::cerr << "[kiwi] warning: " << diag.missing_cell_ids.size()
                  << " input point(s) did not produce a cell (ids:";
        for (size_t i = 0; i < diag.missing_cell_ids.size() && i < 10; ++i)
            std::cerr << " " << diag.missing_cell_ids[i];
        if (diag.missing_cell_ids.size() > 10) std::cerr << " ...";
        std::cerr << ")\n";
    }
    if (!diag.zero_or_negative_volume_ids.empty()) {
        std::cerr << "[kiwi] warning: " << diag.zero_or_negative_volume_ids.size()
                  << " cell(s) with volume <= 0\n";
    }

    kiwi::write_manifest(cfg);
    return 0;
}

int run_map(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    // Chains from point-set (cached) internally -- see run_mapping_stage.
    kiwi::MappingResult result = kiwi::run_mapping_stage(cfg);
    std::cerr << "[kiwi] mapped " << result.cells.size() << " cells\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_write_polaris(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    // Chains from mapping (cached) and tessellation (also cached, see
    // run_tessellation_stage) for neighbor topology -- see writer_polaris.hpp.
    kiwi::PolarisWriteStats stats = kiwi::write_polaris_grid(cfg);
    std::cerr << "[kiwi] wrote " << stats.cells_written << " cells\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_write_voronoi(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    kiwi::VoronoiWriteStats stats = kiwi::write_voronoi_text(cfg);
    std::cerr << "[kiwi] wrote " << stats.cells_written << " cells\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_cartesian_map(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    // Chains from point-set and mapping (both cached) -- see cartesian_map.hpp.
    kiwi::CartesianMap map = kiwi::run_cartesian_map_stage(cfg);
    std::cerr << "[kiwi] resampled onto a " << map.resolution << "^3 Cartesian grid\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_sph_to_grid(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);
    // Note: unlike every other command above, this stage (see
    // cartesian_grid.hpp) doesn't chain from point-set/tessellation/mapping
    // at all -- the manifest check here is only about catching a stale
    // config in a reused run directory, same as everywhere else.

    kiwi::CartesianGridResult grid = kiwi::run_sph_to_grid_stage(cfg);
    std::cerr << "[kiwi] mapped directly onto a " << grid.nx << "x" << grid.ny << "x" << grid.nz
              << " Cartesian grid\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_inverse_mapping_test(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    kiwi::check_manifest_consistency(cfg);

    // Chains from mapping (cached) -- see run_inverse_mapping_test_stage()
    // for why this "testing step" case needs only one extra BFS pass, and
    // for the important caveat on what it can and can't actually validate.
    kiwi::InverseMappingTestResult result = kiwi::run_inverse_mapping_test_stage(cfg);
    std::cerr << "[kiwi] inverse-mapping-test: " << result.particles.size() << " particle(s), "
              << "total original mass " << result.total_original_mass
              << ", total recovered mass " << result.total_recovered_mass << "\n";
    kiwi::write_manifest(cfg);
    return 0;
}

int run_manifest(const std::string& config_path) {
    kiwi::Config cfg = kiwi::Config::load_from_file(config_path);
    std::cerr << "[kiwi] run: " << cfg.run_name << "\n";
    // No stage-consistency check here -- this command's whole purpose is to
    // (re)write the manifest to match the current config and whatever
    // output files already exist, e.g. after hand-editing a config that's
    // known to already match the cached outputs.
    kiwi::write_manifest(cfg);
    std::cerr << "[kiwi] wrote " << cfg.stage_path("manifest.json") << "\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: kiwi <point-set|tessellate|map|write-polaris|write-voronoi|"
                     "cartesian-map|sph-to-grid|inverse-mapping-test|manifest> <config.json>\n";
        return 1;
    }
    std::string stage = argv[1];
    std::string config_path = argv[2];

    try {
        if (stage == "point-set") return run_point_set(config_path);
        if (stage == "tessellate") return run_tessellate(config_path);
        if (stage == "map") return run_map(config_path);
        if (stage == "write-polaris") return run_write_polaris(config_path);
        if (stage == "write-voronoi") return run_write_voronoi(config_path);
        if (stage == "cartesian-map") return run_cartesian_map(config_path);
        if (stage == "sph-to-grid") return run_sph_to_grid(config_path);
        if (stage == "inverse-mapping-test") return run_inverse_mapping_test(config_path);
        if (stage == "manifest") return run_manifest(config_path);
        std::cerr << "unknown stage '" << stage
                  << "' (expected point-set, tessellate, map, write-polaris, write-voronoi, "
                     "cartesian-map, sph-to-grid, inverse-mapping-test, or manifest)\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "[kiwi] error: " << e.what() << "\n";
        return 1;
    }
}
