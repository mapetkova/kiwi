#include "kiwi/diagnostics.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>

#ifndef KIWI_SOURCE_DIR
#error "KIWI_SOURCE_DIR must be defined by the build -- see CMakeLists.txt"
#endif

namespace kiwi {

namespace {
// Wraps `s` in double quotes for use as a single shell argument, escaping
// any embedded '"' or '\'. Defensive only -- every caller in this codebase
// passes a path or a title string we generated ourselves, never raw user
// input, but std::system() still goes through a shell so this is cheap
// insurance against a stray space or quote in a run_name/path.
std::string shell_quote(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    out += "\"";
    return out;
}
} // namespace

void generate_diagnostic_plot(const Config& cfg, const std::string& points_file,
                               size_t n_real, const std::string& title,
                               const std::string& out_filename) {
    std::string script = std::string(KIWI_SOURCE_DIR) + "/tools/plot_diagnostic.py";
    std::string out_path = cfg.stage_path(out_filename);
    const std::string& python = cfg.diagnostics.python_executable;

    std::ostringstream cmd;
    cmd << shell_quote(python) << " " << shell_quote(script)
        << " " << shell_quote(points_file)
        << " " << n_real
        << " " << shell_quote(out_path)
        << " " << shell_quote(title);

    int rc = std::system(cmd.str().c_str());
    if (rc != 0) {
        std::cerr << "[kiwi] warning: diagnostic plot generation failed (\"" << python
                  << "\" exited with code " << rc << "). Continuing without writing " << out_path
                  << ".\n[kiwi]   Try, in a terminal: " << python << " -m pip install --user matplotlib\n"
                  << "[kiwi]   (using \"" << python << " -m pip\" rather than a bare pip3 install "
                  << "avoids installing into a different python3 than the one this ran.) If you "
                  << "have more than one python3 on this machine, set "
                  << "diagnostics.python_executable in the config to the exact interpreter path "
                  << "(e.g. the output of `which python3`).\n";
        return;
    }
    std::cerr << "[kiwi] diagnostic plot written: " << out_path << "\n";
}

} // namespace kiwi
