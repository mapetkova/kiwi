// Regression test: every config file given on the command line must load through
// Config::load_from_file() (schema, types, kernel name, ...) and survive a to_json() round trip.
// Run by ctest on everything in configs/, so a shipped example can't silently go stale.
//
// Usage: test_configs <config.json> [<config.json> ...]
#include <iostream>

#include "kiwi/config.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: test_configs <config.json> [...]\n";
        return 2;
    }
    int failures = 0;
    for (int i = 1; i < argc; ++i) {
        try {
            kiwi::Config cfg = kiwi::Config::load_from_file(argv[i]);
            (void)cfg.to_json();
            (void)cfg.kernel_spec();
            std::cout << "ok    " << argv[i] << "\n";
        } catch (const std::exception& e) {
            std::cout << "FAIL  " << argv[i] << ": " << e.what() << "\n";
            ++failures;
        }
    }
    return failures == 0 ? 0 : 1;
}
