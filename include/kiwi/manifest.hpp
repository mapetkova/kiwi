// Provenance manifest -- records which config and which input/output files
// produced a run directory, so that stale cached output is detected instead
// of silently reused.
//
// Every stage is cached: if its output file already exists under
// cfg.run_dir(), it's loaded instead of recomputed (see point_set.hpp,
// mapping.hpp, etc.). That's fast, but a cached mapping.dat sitting in a run
// directory tells you nothing about which config (tabulated vs analytic
// integral, which density table, ...) produced it, so a config edit followed
// by a re-run could silently reuse stale output computed under the OLD config.
//
// The manifest closes that gap: cfg.to_json() is the single source of truth
// for a run (see config.hpp), so it's stamped into manifest.json alongside
// fingerprints of every configured input and produced output file. Every
// CLI subcommand calls check_manifest_consistency() before running its
// stage and write_manifest() after -- a config change becomes a loud stderr
// warning naming exactly which section changed, instead of a silent stale
// read.
#pragma once

#include <string>

#include "nlohmann/json.hpp"
#include "kiwi/config.hpp"

namespace kiwi {

// Cheap, non-cryptographic file fingerprint: size + mtime + an FNV-1a hash.
// Full-file for anything <= 64 MB; head+tail sampling (32 MB each end) for
// anything larger, so this stays fast on multi-hundred-MB SPH snapshots.
// Purpose is catching "wrong file wired up" bugs, not tamper-proofing.
struct FileFingerprint {
    bool exists = false;
    size_t size_bytes = 0;
    long mtime_unix = 0;
    std::string fingerprint; // hex FNV-1a-64, empty if !exists
};
FileFingerprint fingerprint_file(const std::string& path);

// The run's resolved config (every default made explicit via to_json()),
// fingerprints of the configured inputs (input_snapshot, and
// mapping.density_table when use_tabulated_integral is set), a manifest
// write timestamp, and fingerprints of every known stage-output file
// present in cfg.run_dir() (points-final.dat, tessellation.dat, mapping.dat,
// grid-polaris.dat, grid-voronoi.txt -- whichever exist).
nlohmann::json build_manifest(const Config& cfg);

// Writes/overwrites cfg.stage_path("manifest.json"). Call after a stage
// finishes so the manifest reflects the config that most recently touched
// this run directory.
void write_manifest(const Config& cfg);

// Reads any existing manifest.json in cfg.run_dir() and compares its stored
// config against cfg's current config, top-level section by section.
// Differences are printed as stderr warnings (not thrown -- a stage caching
// decision is still made by that stage's own cache check, this only makes
// the risk visible). Returns true if a mismatch was found; false if the
// configs match or no manifest exists yet (first run in this directory).
bool check_manifest_consistency(const Config& cfg);

} // namespace kiwi
