#include "kiwi/manifest.hpp"

#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <vector>

namespace kiwi {

namespace {

constexpr size_t kFullHashLimit = 64ull * 1024 * 1024;  // 64 MB
constexpr size_t kSampleSize = 32ull * 1024 * 1024;      // 32 MB from each end

uint64_t fnv1a_update(uint64_t h, const char* data, size_t n) {
    constexpr uint64_t prime = 1099511628211ull;
    for (size_t i = 0; i < n; ++i) {
        h ^= static_cast<unsigned char>(data[i]);
        h *= prime;
    }
    return h;
}

std::string to_hex(uint64_t v) {
    std::ostringstream oss;
    oss << std::hex << v;
    return oss.str();
}

// Known output filenames every stage may have cached under cfg.run_dir() --
// see point_set.hpp/tessellation.hpp/mapping.hpp/writer_polaris.hpp/
// writer_voronoi.hpp/cartesian_grid.hpp for where each is written.
const std::vector<std::string>& known_output_files() {
    static const std::vector<std::string> names = {
        "points-final.dat", "tessellation.dat", "mapping.dat",
        "grid-polaris.dat", "grid-voronoi.txt", "cartesian_map.dat",
        "cartesian_grid.dat", "buffer_diagnostic.png", "lloyd_diagnostic.png"};
    return names;
}

} // namespace

FileFingerprint fingerprint_file(const std::string& path) {
    FileFingerprint fp;
    if (path.empty()) return fp;

    struct stat st;
    if (stat(path.c_str(), &st) != 0) return fp; // doesn't exist -- leave exists=false

    fp.exists = true;
    fp.size_bytes = static_cast<size_t>(st.st_size);
    fp.mtime_unix = static_cast<long>(st.st_mtime);

    std::ifstream f(path, std::ios::binary);
    if (!f) return fp; // stat succeeded but couldn't open -- report metadata, no hash

    uint64_t h = 1469598103934665603ull; // FNV offset basis
    std::vector<char> buf(1 << 20); // 1 MB read buffer

    if (fp.size_bytes <= kFullHashLimit) {
        while (f) {
            f.read(buf.data(), static_cast<std::streamsize>(buf.size()));
            std::streamsize got = f.gcount();
            if (got > 0) h = fnv1a_update(h, buf.data(), static_cast<size_t>(got));
        }
    } else {
        // Head sample.
        size_t remaining = kSampleSize;
        while (remaining > 0 && f) {
            size_t chunk = std::min(remaining, buf.size());
            f.read(buf.data(), static_cast<std::streamsize>(chunk));
            std::streamsize got = f.gcount();
            if (got <= 0) break;
            h = fnv1a_update(h, buf.data(), static_cast<size_t>(got));
            remaining -= static_cast<size_t>(got);
        }
        // Tail sample.
        f.seekg(static_cast<std::streamoff>(fp.size_bytes - kSampleSize), std::ios::beg);
        remaining = kSampleSize;
        while (remaining > 0 && f) {
            size_t chunk = std::min(remaining, buf.size());
            f.read(buf.data(), static_cast<std::streamsize>(chunk));
            std::streamsize got = f.gcount();
            if (got <= 0) break;
            h = fnv1a_update(h, buf.data(), static_cast<size_t>(got));
            remaining -= static_cast<size_t>(got);
        }
        // Mix the total size in too, so a same-sized-sample-but-different-middle
        // file (rare, but possible) at least changes the fingerprint via size,
        // which is already recorded separately -- this just keeps the hash
        // itself size-sensitive as well.
        h = fnv1a_update(h, reinterpret_cast<const char*>(&fp.size_bytes), sizeof(fp.size_bytes));
    }

    fp.fingerprint = to_hex(h);
    return fp;
}

namespace {
nlohmann::json fingerprint_to_json(const FileFingerprint& fp) {
    nlohmann::json j;
    j["exists"] = fp.exists;
    if (fp.exists) {
        j["size_bytes"] = fp.size_bytes;
        j["mtime_unix"] = fp.mtime_unix;
        j["fingerprint"] = fp.fingerprint;
    }
    return j;
}
} // namespace

nlohmann::json build_manifest(const Config& cfg) {
    nlohmann::json j;
    j["manifest_schema_version"] = 1;
    j["written_at_unix"] = static_cast<long>(std::time(nullptr));
    j["config"] = cfg.to_json();

    nlohmann::json inputs;
    inputs["input_snapshot"] = fingerprint_to_json(fingerprint_file(cfg.input_snapshot));
    if (cfg.mapping.use_tabulated_integral && !cfg.mapping.density_table.empty()) {
        // fingerprint_file() already degrades to {"exists": false} rather
        // than throwing if this path isn't there -- which is a
        // meaningful, non-error state: run_mapping_stage() auto-generates
        // the lookup table in memory in that case (see mapping.cpp) rather
        // than failing, so an "exists": false entry here records that this
        // run used the generated table, not a broken config.
        inputs["density_table"] = fingerprint_to_json(fingerprint_file(cfg.mapping.density_table));
    }
    if (cfg.cartesian_grid.use_tabulated_integral && !cfg.cartesian_grid.density_table.empty()) {
        // Same auto-generate-in-memory fallback as mapping.density_table
        // above -- see run_sph_to_grid_stage() (cartesian_grid.cpp).
        inputs["cartesian_grid_density_table"] =
            fingerprint_to_json(fingerprint_file(cfg.cartesian_grid.density_table));
    }
    j["inputs"] = inputs;

    nlohmann::json outputs;
    for (const auto& name : known_output_files()) {
        FileFingerprint fp = fingerprint_file(cfg.stage_path(name));
        if (fp.exists) outputs[name] = fingerprint_to_json(fp);
    }
    j["outputs"] = outputs;

    return j;
}

void write_manifest(const Config& cfg) {
    std::string path = cfg.stage_path("manifest.json");
    std::ofstream f(path);
    if (!f) {
        std::cerr << "[kiwi] warning: could not write manifest to " << path << "\n";
        return;
    }
    f << build_manifest(cfg).dump(2) << "\n";
}

bool check_manifest_consistency(const Config& cfg) {
    std::string path = cfg.stage_path("manifest.json");
    std::ifstream f(path);
    if (!f) return false; // no prior manifest -- first run in this directory

    nlohmann::json prior;
    try {
        f >> prior;
    } catch (const nlohmann::json::parse_error&) {
        std::cerr << "[kiwi] warning: existing manifest.json is not valid JSON, ignoring: " << path << "\n";
        return false;
    }

    if (!prior.contains("config")) return false;
    nlohmann::json prior_cfg = prior["config"];
    nlohmann::json cur_cfg = cfg.to_json();

    std::vector<std::string> mismatched;
    for (auto it = cur_cfg.begin(); it != cur_cfg.end(); ++it) {
        const std::string& key = it.key();
        if (!prior_cfg.contains(key) || prior_cfg[key] != it.value()) {
            mismatched.push_back(key);
        }
    }

    if (mismatched.empty()) return false;

    std::cerr << "[kiwi] WARNING: this run directory (" << cfg.run_dir() << ") was previously "
              << "used with a DIFFERENT config. Cached stage output here may have been computed "
              << "under the old settings. Changed section(s): ";
    for (size_t i = 0; i < mismatched.size(); ++i) {
        std::cerr << mismatched[i] << (i + 1 < mismatched.size() ? ", " : "");
    }
    std::cerr << "\n[kiwi]   old: " << prior_cfg.dump() << "\n";
    std::cerr << "[kiwi]   new: " << cur_cfg.dump() << "\n";
    std::cerr << "[kiwi] Delete the stale cached file(s) under that directory (or use a "
                 "different work_folder/run_name) if you intended this change to take effect.\n";

    return true;
}

} // namespace kiwi
