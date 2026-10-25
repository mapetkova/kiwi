#include <iostream>
#include <stdexcept>
#include <string>

#include "kiwi/cartesian_grid.hpp"
#include "kiwi/cartesian_map.hpp"
#include "kiwi/config.hpp"
#include "kiwi/manifest.hpp"
#include "kiwi/mapping.hpp"
#include "kiwi/point_set.hpp"
#include "kiwi/snapshot.hpp"
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
    // the voro++ container from scratch. Delete tessellation.dat (or the run
    // directory) to force a recompute.
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
    std::cerr << "[kiwi] note: write-polaris is experimental -- the written quantities (gas temperature, "
                 "velocity, density) are fixed and the interface may change\n";

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

// `kiwi convert-snapshot <snapshot> <output> [options]`: read a snapshot in any supported format and write it
// as plain_ascii (gas particles only; Phantom dumps in pc, Msun and km/s unless overridden). It lets
// tools that don't read the binary formats (tools/suggest_config.py) work on a Phantom dump.
static int run_convert_snapshot(int argc, char** argv) {
    const char* usage =
        "usage: kiwi convert-snapshot <snapshot> <output.txt> [--format auto|plain_ascii|splash_ascii|phantom_binary]\n"
        "                             [--allow-missing-velocity] [--length-cm X] [--mass-g X] [--velocity-cms X]\n";
    if (argc < 4) {
        std::cerr << usage;
        return 1;
    }
    std::string input = argv[2], output = argv[3], format = "auto";
    bool allow_missing_velocity = false;
    kiwi::PhantomUnits units;
    for (int i = 4; i < argc; ++i) {
        std::string a = argv[i];
        auto value = [&](const char* name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error(std::string(name) + " needs a value");
            return argv[++i];
        };
        if (a == "--format") format = value("--format");
        else if (a == "--allow-missing-velocity") allow_missing_velocity = true;
        else if (a == "--length-cm") units.length_cm = std::stod(value("--length-cm"));
        else if (a == "--mass-g") units.mass_g = std::stod(value("--mass-g"));
        else if (a == "--velocity-cms") units.velocity_cms = std::stod(value("--velocity-cms"));
        else {
            std::cerr << "unknown option '" << a << "'\n" << usage;
            return 1;
        }
    }
    kiwi::Config defaults;  // only for the default kernel's hfact (used by splash_ascii files that lack h or rho)
    kiwi::SPHSnapshot snap =
        kiwi::load_sph_snapshot(input, format, defaults.kernel_hfact(), units, allow_missing_velocity);
    kiwi::write_plain_ascii(snap, output);
    std::cerr << "[kiwi] wrote " << snap.size() << " particle(s) to " << output << "\n";
    return 0;
}

#ifndef KIWI_VERSION
#error "KIWI_VERSION must be defined by the build -- see CMakeLists.txt"
#endif

int main(int argc, char** argv) {
    if (argc == 2 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-V")) {
        std::cout << "kiwi " << KIWI_VERSION << "\n";
        return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "convert-snapshot") {
        try {
            return run_convert_snapshot(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << "[kiwi] error: " << e.what() << "\n";
            return 1;
        }
    }
    if (argc != 3) {
        std::cerr << "usage: kiwi <point-set|tessellate|map|write-polaris|write-voronoi|"
                     "cartesian-map|sph-to-grid|inverse-mapping-test|manifest> <config.json>\n"
                     "       kiwi convert-snapshot <snapshot> <output.txt> [options]\n"
                     "       kiwi --version\n";
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
