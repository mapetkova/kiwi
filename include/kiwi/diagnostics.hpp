// Optional diagnostic plots -- a quick visual sanity check of where buffer
// points landed and how Lloyd relaxation moved the sites, without needing to
// open points-buff.dat/points-final.dat by hand. Off by default
// (diagnostics.enabled in the config); when on, the point-set stage calls
// generate_diagnostic_plot() at two checkpoints: right after the buffer
// stage (if buffer.enabled) and right after Lloyd relaxation finishes (if
// lloyd.iterations > 0).
//
// Plotting itself is delegated to a bundled Python/matplotlib script
// (tools/plot_diagnostic.py) rather than reimplemented in C++ -- this is
// best-effort output, not part of the pipeline's data contract, so pulling a
// C++ plotting library into the build isn't worth it. If python3 or
// matplotlib isn't available, generate_diagnostic_plot() prints a warning
// and returns; it never throws, so a missing plotting toolchain can't break
// a real run (same "warn, don't block" pattern as the OpenMP-not-found and
// manifest-mismatch cases elsewhere).
#pragma once

#include <cstddef>
#include <string>

#include "kiwi/config.hpp"

namespace kiwi {

// Reads "id x y z" rows from `points_file` (any file written by
// save_point_set -- points-buff.dat or points-final.dat), splits them into
// "real" (id <= n_real) and "buffer" (id > n_real) using the same id
// convention every buffer strategy in point_set.cpp relies on (buffer ids
// start at n_real + 1), and writes an x-y scatter plot to
// cfg.stage_path(out_filename). `title` is passed through to the plot
// as-is (may contain '\n' for a multi-line title).
void generate_diagnostic_plot(const Config& cfg, const std::string& points_file,
                               size_t n_real, const std::string& title,
                               const std::string& out_filename);

} // namespace kiwi
