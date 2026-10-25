# Contributing to KIWI

Thanks for your interest. KIWI is pre-beta, so for anything bigger than a small fix please open an issue
first to talk it through; that avoids work on something that does not fit.

## Reporting a bug

Please include:

- the output of `kiwi --version`, your OS and compiler (and whether OpenMP was found);
- the command you ran and the **complete** config JSON;
- the snapshot format, and, if you can share it, a small snapshot that reproduces the problem (a cutout of
  a few thousand particles is ideal, like `tests/data/brick_03927_cutout.ascii`, which keeps only the rows in
  a small x-y box);
- the full stderr output. Errors start with `[kiwi] error:`; warnings with `[kiwi] warning:`.

## Building and testing

You need CMake >= 3.13 and a C++17 compiler. Python 3 (standard library only) is needed for most tests;
the plot and config helpers in `tools/` also need `numpy` and `matplotlib`.

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure -LE slow
```

`-LE slow` skips `real_data`, which needs a ~350 MB snapshot that is not in the repository (it is skipped
automatically when the file is absent); the same checks run on a small cutout as `real_data_cutout`. CI
(`.github/workflows/ci.yml`) runs the tests on Ubuntu and macOS, with OpenMP at 1, 2 and 4 threads, and
with AddressSanitizer and UBSan, so please check that your change builds without new warnings and passes
on at least one platform before opening a pull request.

## Layout

| Path | Contents |
|---|---|
| `src/`, `include/kiwi/` | the library (one stage per file pair) and `main.cpp`, the command-line front end |
| `tests/` | the ctest tests (Python drivers around the `kiwi` executable) and `CMakeLists.txt` listing them |
| `tools/` | helper scripts, the kernel-table generator (`kernel_tables/`) and C++ test programs |
| `configs/` | an example config |
| `third_party/` | bundled Voro++ and nlohmann/json. Please do not edit these; see `THIRD_PARTY.md` |

## Making a change

- **Add a test with every behaviour change.** Numerical results are checked against known answers
  (mass and momentum conservation, exact lattices, brute-force quadrature), not stored output. If you fix
  a bug, add a test that fails without the fix.
- **Errors**: invalid input should exit with status 1 and one `[kiwi] error: ...` line that names the
  problem and, for a bad row in a snapshot, the line number. `tests/bad_input_test.py` lists the cases
  covered; add yours there.
- **Config options**: add the field in `include/kiwi/config.hpp`, parse and validate it in
  `src/config.cpp`, include it in `to_json()` (the manifest compares it), document it in the README and
  add a case to `bad_input_test.py` if it can be invalid. Every file in `configs/` must load
  (`example_configs` test).
- **Comments** explain why, not how. Describe what the code does now; the history of a change belongs in
  the commit message and `CHANGELOG.md`.
- Follow the style of the surrounding code (4-space indentation in C++ and Python, snake_case names).
  There is no formatter configuration.
- Keep new dependencies out of the C++ build. Python helpers should stay within the standard library,
  except where `numpy` or `matplotlib` is clearly needed.

### Kernels

The kernel-integral coefficients in `src/kernel_spec_data.cpp` are generated, not hand-written. To change
or add a kernel, edit the definitions in `tools/kernel_tables/` (see its README), run

```bash
python3 tools/kernel_tables/gen_cpp.py
```

and commit the regenerated file. The `kernel_spec_data_current` test fails if the committed file does not
match the generator (`gen_cpp.py --check`), and `kernel_tables` and `generic_kernel` check the
coefficients against brute-force quadrature.

## Releases

Bump the version in `CMakeLists.txt` (`project(kiwi VERSION ...)`) and in `CITATION.cff` together (the
`version_consistent` test checks that they agree), add the date to `CITATION.cff` (`date-released`), move
the entries in `CHANGELOG.md` from Unreleased to the new version, and tag `vX.Y.Z`.

## Licence

KIWI is distributed under the GNU General Public License, version 3 or later (see `LICENSE`). By
contributing you agree that your contribution is licensed under the same terms.
