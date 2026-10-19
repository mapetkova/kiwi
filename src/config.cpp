#include "kiwi/config.hpp"

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

namespace kiwi {

// Strips directory components and the last extension only (e.g.
// "/a/b/ascii_CMZ_..._00147.out" -> "ascii_CMZ_..._00147"). Deliberately
// simple -- multi-dot filenames keep everything before the *last* dot, which
// matches how every snapshot file we've seen is named.
static std::string trimmed_basename(const std::string& path) {
    std::string base = path;
    auto slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
    auto dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return base;
}

static void mkdir_p(const std::string& path) {
    std::string partial;
    for (size_t i = 0; i < path.size(); ++i) {
        if (path[i] == '/') {
            if (!partial.empty()) mkdir(partial.c_str(), 0755);
        }
        partial += path[i];
    }
    mkdir(partial.c_str(), 0755);
}

using nlohmann::json;

void Domain::validate() const {
    if (!(x_min < x_max)) throw ConfigError("domain.x_min must be < domain.x_max");
    if (!(y_min < y_max)) throw ConfigError("domain.y_min must be < domain.y_max");
    if (!(z_min < z_max)) throw ConfigError("domain.z_min must be < domain.z_max");
}

static Domain domain_from_json(const json& j) {
    Domain d;
    d.x_min = j.at("x_min").get<double>();
    d.x_max = j.at("x_max").get<double>();
    d.y_min = j.at("y_min").get<double>();
    d.y_max = j.at("y_max").get<double>();
    d.z_min = j.at("z_min").get<double>();
    d.z_max = j.at("z_max").get<double>();
    return d;
}

static json domain_to_json(const Domain& d) {
    return json{{"x_min", d.x_min}, {"x_max", d.x_max},
                {"y_min", d.y_min}, {"y_max", d.y_max},
                {"z_min", d.z_min}, {"z_max", d.z_max}};
}

Config Config::load_from_file(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw ConfigError("could not open config file: " + path);

    json j;
    try {
        f >> j;
    } catch (const json::parse_error& e) {
        throw ConfigError("config file is not valid JSON (" + path + "): " + e.what());
    }

    Config c;
    try {
        c.run_name = j.at("run_name").get<std::string>();
        c.domain = domain_from_json(j.at("domain"));
        c.domain.validate();

        if (j.contains("grid_blocks")) {
            c.grid_blocks.block_scale_factor = j["grid_blocks"].value("block_scale_factor", 3.0);
            if (c.grid_blocks.block_scale_factor <= 0.0)
                throw ConfigError("grid_blocks.block_scale_factor must be > 0");
        }

        c.input_snapshot = j.at("input_snapshot").get<std::string>();
        c.input_snapshot_format = j.value("input_snapshot_format", std::string("auto"));
        if (c.input_snapshot_format != "auto" && c.input_snapshot_format != "dale_phantom" &&
            c.input_snapshot_format != "splash_ascii" && c.input_snapshot_format != "phantom_binary")
            throw ConfigError("input_snapshot_format must be 'auto', 'dale_phantom', 'splash_ascii', or "
                               "'phantom_binary' (got '" + c.input_snapshot_format + "')");
        c.work_folder = j.at("work_folder").get<std::string>();

        if (j.contains("buffer")) {
            const auto& b = j["buffer"];
            c.buffer.enabled = b.value("enabled", false);
            c.buffer.strategy = b.value("strategy", std::string("legacy_brick"));
            if (c.buffer.enabled && c.buffer.strategy == "legacy_brick") {
                if (!b.contains("legacy_brick"))
                    throw ConfigError("buffer.strategy is legacy_brick but buffer.legacy_brick params are missing");
                const auto& lb = b.at("legacy_brick");
                c.buffer.legacy_brick.dx = lb.at("dx").get<double>();
                c.buffer.legacy_brick.dy = lb.at("dy").get<double>();
                c.buffer.legacy_brick.dz = lb.at("dz").get<double>();
                c.buffer.legacy_brick.shell_spacing_1 = lb.value("shell_spacing_1", 3.0);
                c.buffer.legacy_brick.shell_spacing_2 = lb.value("shell_spacing_2", 4.0);
                c.buffer.legacy_brick.shell_spacing_3 = lb.value("shell_spacing_3", 5.0);
            }
            if (c.buffer.enabled && c.buffer.strategy == "coverage_grid") {
                const auto& cg = b.value("coverage_grid", json::object());
                c.buffer.coverage_grid.spacing_fraction = cg.value("spacing_fraction", 0.05);
                c.buffer.coverage_grid.coverage_radius_factor = cg.value("coverage_radius_factor", 2.0);
                if (c.buffer.coverage_grid.spacing_fraction <= 0.0)
                    throw ConfigError("buffer.coverage_grid.spacing_fraction must be > 0");
                if (c.buffer.coverage_grid.coverage_radius_factor <= 0.0)
                    throw ConfigError("buffer.coverage_grid.coverage_radius_factor must be > 0");
            }
        }

        if (j.contains("lloyd")) {
            c.lloyd.iterations = j["lloyd"].value("iterations", 0);
            c.lloyd.save_each_iteration = j["lloyd"].value("save_each_iteration", false);
            if (c.lloyd.iterations < 0)
                throw ConfigError("lloyd.iterations must be >= 0");
        }

        if (j.contains("mapping")) {
            c.mapping.use_tabulated_integral = j["mapping"].value("use_tabulated_integral", false);
            c.mapping.density_table = j["mapping"].value("density_table", std::string());
        }

        if (j.contains("diagnostics")) {
            c.diagnostics.enabled = j["diagnostics"].value("enabled", false);
            c.diagnostics.python_executable =
                j["diagnostics"].value("python_executable", std::string("python3"));
            if (c.diagnostics.python_executable.empty())
                throw ConfigError("diagnostics.python_executable must not be empty");
        }

        if (j.contains("cartesian_map")) {
            c.cartesian_map.resolution = j["cartesian_map"].value("resolution", 128);
            if (c.cartesian_map.resolution <= 0)
                throw ConfigError("cartesian_map.resolution must be > 0");
        }

        if (j.contains("cartesian_grid")) {
            const auto& cgc = j["cartesian_grid"];
            c.cartesian_grid.resolution = cgc.value("resolution", 128);
            c.cartesian_grid.resolution_x = cgc.value("resolution_x", 0);
            c.cartesian_grid.resolution_y = cgc.value("resolution_y", 0);
            c.cartesian_grid.resolution_z = cgc.value("resolution_z", 0);
            c.cartesian_grid.use_tabulated_integral = cgc.value("use_tabulated_integral", false);
            c.cartesian_grid.density_table = cgc.value("density_table", std::string());
            if (c.cartesian_grid.resolution <= 0)
                throw ConfigError("cartesian_grid.resolution must be > 0");
            if (c.cartesian_grid.resolution_x < 0 || c.cartesian_grid.resolution_y < 0 ||
                c.cartesian_grid.resolution_z < 0)
                throw ConfigError("cartesian_grid.resolution_x/resolution_y/resolution_z must be >= 0 "
                                   "(0 means: fall back to resolution)");
        }

        if (j.contains("kernel")) {
            const auto& k = j["kernel"];
            c.kernel.name = k.value("name", std::string("cubic_spline"));
            c.kernel.support_radius_in_h = k.value("support_radius_in_h", 0.0);
            c.kernel.hfact = k.value("hfact", 0.0);
            c.kernel.generic_path = k.value("generic_path", false);
            try {
                (void)get_kernel(c.kernel.name);
            } catch (const std::invalid_argument& e) {
                throw ConfigError(std::string("kernel.name: ") + e.what());
            }
            if (c.kernel.support_radius_in_h < 0.0)
                throw ConfigError("kernel.support_radius_in_h must be > 0 (or 0 / unset for the kernel's default)");
            if (c.kernel.hfact < 0.0)
                throw ConfigError("kernel.hfact must be > 0 (or 0 / unset for the kernel's default)");
        }
        if (c.use_generic_kernel_path() &&
            (c.mapping.use_tabulated_integral || c.cartesian_grid.use_tabulated_integral))
            throw ConfigError("use_tabulated_integral is only available for the cubic spline kernel with its "
                              "default support radius (and kernel.generic_path unset); the selected kernel "
                              "settings need the analytic integral -- set use_tabulated_integral to false");

        if (j.contains("phantom_units")) {
            const auto& u = j["phantom_units"];
            c.phantom_units.length_cm = u.value("length_cm", c.phantom_units.length_cm);
            c.phantom_units.mass_g = u.value("mass_g", c.phantom_units.mass_g);
            c.phantom_units.velocity_cms = u.value("velocity_cms", c.phantom_units.velocity_cms);
            if (c.phantom_units.length_cm <= 0.0 || c.phantom_units.mass_g <= 0.0 || c.phantom_units.velocity_cms <= 0.0)
                throw ConfigError("phantom_units.length_cm / mass_g / velocity_cms must all be > 0");
        }

        if (j.contains("polaris_writer")) {
            const auto& pw = j["polaris_writer"];
            c.polaris_writer.sim_id = pw.value("sim_id", (unsigned short)50);
            c.polaris_writer.length_scale_m = pw.value("length_scale_m", 3.08567758e16);
            c.polaris_writer.velocity_scale = pw.value("velocity_scale", 1000.0);
            c.polaris_writer.density_scale = pw.value("density_scale",
                1.9891e30 / (3.08567758e16 * 3.08567758e16 * 3.08567758e16));
            c.polaris_writer.gas_temperature_k = pw.value("gas_temperature_k", 65.0);
        }

    } catch (const json::out_of_range& e) {
        throw ConfigError(std::string("config missing required field: ") + e.what());
    }

    return c;
}

json Config::to_json() const {
    json j;
    j["run_name"] = run_name;
    j["domain"] = domain_to_json(domain);
    j["grid_blocks"] = {{"block_scale_factor", grid_blocks.block_scale_factor}};
    j["input_snapshot"] = input_snapshot;
    j["input_snapshot_format"] = input_snapshot_format;
    j["work_folder"] = work_folder;

    json b;
    b["enabled"] = buffer.enabled;
    b["strategy"] = buffer.strategy;
    if (buffer.strategy == "legacy_brick") {
        b["legacy_brick"] = {{"dx", buffer.legacy_brick.dx},
                              {"dy", buffer.legacy_brick.dy},
                              {"dz", buffer.legacy_brick.dz},
                              {"shell_spacing_1", buffer.legacy_brick.shell_spacing_1},
                              {"shell_spacing_2", buffer.legacy_brick.shell_spacing_2},
                              {"shell_spacing_3", buffer.legacy_brick.shell_spacing_3}};
    } else if (buffer.strategy == "coverage_grid") {
        b["coverage_grid"] = {{"spacing_fraction", buffer.coverage_grid.spacing_fraction},
                               {"coverage_radius_factor", buffer.coverage_grid.coverage_radius_factor}};
    }
    j["buffer"] = b;

    j["lloyd"] = {{"iterations", lloyd.iterations}, {"save_each_iteration", lloyd.save_each_iteration}};
    j["mapping"] = {{"use_tabulated_integral", mapping.use_tabulated_integral},
                     {"density_table", mapping.density_table}};
    j["polaris_writer"] = {{"sim_id", polaris_writer.sim_id},
                            {"length_scale_m", polaris_writer.length_scale_m},
                            {"velocity_scale", polaris_writer.velocity_scale},
                            {"density_scale", polaris_writer.density_scale},
                            {"gas_temperature_k", polaris_writer.gas_temperature_k}};
    j["diagnostics"] = {{"enabled", diagnostics.enabled},
                         {"python_executable", diagnostics.python_executable}};
    j["cartesian_map"] = {{"resolution", cartesian_map.resolution}};
    j["cartesian_grid"] = {{"resolution", cartesian_grid.resolution},
                            {"resolution_x", cartesian_grid.resolution_x},
                            {"resolution_y", cartesian_grid.resolution_y},
                            {"resolution_z", cartesian_grid.resolution_z},
                            {"use_tabulated_integral", cartesian_grid.use_tabulated_integral},
                            {"density_table", cartesian_grid.density_table}};
    // Only written when something differs from the default, so the config (and the provenance manifest
    // built from it) of a pre-existing cubic-spline run is unchanged by the existence of this section.
    if (!kernel.is_default()) {
        j["kernel"] = {{"name", kernel.name},
                       {"support_radius_in_h", kernel.support_radius_in_h},
                       {"hfact", kernel.hfact},
                       {"generic_path", kernel.generic_path}};
    }
    if (!phantom_units.is_default()) {
        j["phantom_units"] = {{"length_cm", phantom_units.length_cm},
                              {"mass_g", phantom_units.mass_g},
                              {"velocity_cms", phantom_units.velocity_cms}};
    }
    return j;
}

const KernelSpec& Config::kernel_spec() const { return get_kernel(kernel.name); }

double Config::kernel_support_in_h() const {
    return kernel.support_radius_in_h > 0.0 ? kernel.support_radius_in_h : kernel_spec().radkern;
}

double Config::kernel_hfact() const { return kernel.hfact > 0.0 ? kernel.hfact : kernel_spec().hfact; }

bool Config::use_generic_kernel_path() const {
    return kernel.generic_path || std::string(kernel_spec().key) != "M4" || kernel_support_in_h() != 2.0;
}

std::string Config::run_dir() const {
    return work_folder + "/" + trimmed_basename(input_snapshot);
}

std::string Config::stage_path(const std::string& filename) const {
    std::string dir = run_dir();
    mkdir_p(dir);
    return dir + "/" + filename;
}

} // namespace kiwi
