# KIWI

**KIWI** stands for **Kernel-Integral Weighted Interpolation**.

KIWI maps SPH (smoothed-particle hydrodynamics) simulation snapshots onto grids. Instead of
sampling the density at a point, it integrates each particle's smoothing kernel over the volume of
every grid cell, so mass is conserved cell by cell and the result doesn't depend on how the
particles happen to fall relative to the cells. The mapping methods follow
Petkova et al. 2018 ([doi:10.1016/j.jcp.2017.10.024](https://doi.org/10.1016/j.jcp.2017.10.024))
and Petkova et al. 2021 ([doi:10.1093/mnras/stab2178](https://doi.org/10.1093/mnras/stab2178)).

It can produce:

- **Density and velocity on a Voronoi grid.** KIWI builds the Voronoi tessellation from the
  (relaxed) particle positions with [Voro++](https://math.lbl.gov/voro++/) and maps the SPH data
  onto its cells.
- **Density and velocity on a Cartesian grid**, either resampled from the Voronoi grid
  (`cartesian-map`) or computed directly from the particles (`sph-to-grid`).
- **As an optional add-on**, the Voronoi grid written in the input format of the radiative-transfer
  code [POLARIS](https://www1.astrophysik.uni-kiel.de/~polaris/), or as plain text.

One JSON config per run drives every stage, and each stage caches its output, so re-running a
config only recomputes what is missing.

> **Status: pre-beta.** The mapping stages have been validated against earlier tools, including
> SPLASH, but the interface may still change. Bug reports are welcome.

## Build

Requires CMake >= 3.13 and a C++17 compiler. Voro++ and nlohmann/json are bundled in
`third_party/`, so there is nothing else to install.

```
mkdir build && cd build
cmake ..
cmake --build . -j
```

This produces `build/kiwi`. OpenMP is optional and speeds up the mapping stage. On macOS the
system compiler has no OpenMP, so use Homebrew's:

```
brew install llvm libomp
cmake .. -DCMAKE_CXX_COMPILER=$(brew --prefix llvm)/bin/clang++ -DOpenMP_ROOT=$(brew --prefix libomp)
```

## Quick start

1. Copy one of the configs in `configs/` and edit the paths. `hdens-00282-smoke-test.json` is
   the fastest. Set `input_snapshot` to your snapshot and `domain` to the box you want to map.
2. Run a stage, from the repo root:

```
build/kiwi map <config.json>
```

Each stage runs the stages before it if their output is missing:

| Command | What it does | Output (in `work_folder/<snapshot>/`) |
|---|---|---|
| `kiwi point-set <cfg>` | Load the snapshot, optionally add buffer points, relax with Lloyd iterations | `points-final.dat` |
| `kiwi tessellate <cfg>` | Voronoi tessellation of the points | `tessellation.dat` |
| `kiwi map <cfg>` | Kernel-integral density and velocity per Voronoi cell | `mapping.dat` |
| `kiwi cartesian-map <cfg>` | Resample the Voronoi result onto a Cartesian grid | `cartesian_map.dat` |
| `kiwi sph-to-grid <cfg>` | Map particles directly onto a Cartesian grid (no Voronoi step) | `cartesian_grid.dat` |
| `kiwi write-voronoi <cfg>` | Plain-text Voronoi grid | `grid-voronoi.txt` |
| `kiwi write-polaris <cfg>` | Optional: POLARIS binary grid | `grid-polaris.dat` |
| `kiwi inverse-mapping-test <cfg>` | Validation: map the grid back onto the particles | see log |
| `kiwi manifest <cfg>` | Write `manifest.json` recording the settings behind each output | `manifest.json` |

The output of `sph-to-grid` is currently in KIWI's own simple binary format, which is still a
placeholder, and `tools/plot_cartesian_map.py` cannot read it yet. `tools/plot_cartesian_map.py`
plots the `cartesian-map` output.

## Minimal config

```json
{
  "run_name": "my-run",
  "domain": {"x_min": -10, "x_max": 10, "y_min": -10, "y_max": 10, "z_min": -10, "z_max": 10},
  "input_snapshot": "path/to/snapshot.out",
  "work_folder": "runs",
  "buffer": {"enabled": false},
  "lloyd": {"iterations": 1, "save_each_iteration": false},
  "cartesian_map": {"resolution": 128}
}
```

Useful options:

- `buffer.strategy`: `"coverage_grid"` (general purpose) or `"legacy_brick"` (specific to one
  cloud geometry). A buffer of extra points keeps the cells at the domain edge finite.
- `kernel`: the SPH kernel the snapshot was run with. Default is the cubic spline; the others are
  `quartic_spline`, `quintic_spline`, `wendland_c2`, `wendland_c4` and `wendland_c6` (all analytic;
  the tabulated integral exists only for the cubic spline):

  ```json
  "kernel": {"name": "wendland_c4"}
  ```

  The kernel's support radius is taken in units of the snapshot's smoothing length `h` following SPLASH
  and Phantom (2h for the cubic spline and the Wendland kernels, 2.5h for the quartic, 3h for the
  quintic); set `kernel.support_radius_in_h` if your `h` follows another convention. `kernel.hfact`
  (default: Phantom's value for the chosen kernel) is only used to derive `h` or the density when a
  `splash_ascii` snapshot lacks one of them. With a kernel wider than 2h, consider raising
  `buffer.coverage_grid.coverage_radius_factor` (default 2) to match. `kernel.generic_path: true` runs
  even the cubic spline through the generic integral, which should agree with the default path to
  about 1e-11.
- `mapping.use_tabulated_integral`: `false` (default) uses the analytic kernel integral. `true`
  uses a lookup table, which is faster but approximate.
- `polaris_writer`: unit scales and temperature for the optional POLARIS output. The domain must
  be a cube if you use it.
- `diagnostics.enabled`: writes diagnostic plots (needs Python with matplotlib).

## Input snapshots

`input_snapshot_format` is `"auto"` by default. Two ASCII formats are supported:

- `splash_ascii`: the output of `splash to ascii`. Columns are read by name from the header, so
  the exact column set and order can vary. Required: `x y z m vx vy vz` and one of `h` or density.
- `dale_phantom`: headerless, 9 columns per particle: `x y z vx vy vz m h rho`.

`tools/suggest_config.py` (below) reads either format and proposes a config for a snapshot.

## Setting up a run: `suggest_config.py`

`tools/suggest_config.py` builds a ready-to-run config from a raw snapshot, so you don't have to
pick a domain by hand or copy an old config. It needs Python 3 with `numpy` (and `matplotlib` for
the interactive zoom).

```
python3 tools/suggest_config.py path/to/snapshot.out
```

It reads the snapshot (format auto-detected), suggests a **cubic** domain that covers the particles
(trimming a tiny fraction of far-flung outliers, and reporting how many particles fall outside),
asks you for the few things it can't guess (`run_name`, `work_folder`, and whether to use a
tabulated kernel integral), and writes `config.json`. Everything else gets sensible defaults, with a
`coverage_grid` buffer, no Lloyd iterations and a 128-cell Cartesian map, which you can edit
afterwards.

**Zooming in on part of a snapshot.** Use `--mode zoom` to map a smaller cubic region instead of the
whole snapshot. The centre can be the geometric centre, the centre of mass, the densest region, or
coordinates you give:

```
python3 tools/suggest_config.py snapshot.out --mode zoom --center-mode densest --box-side 8 \
    --run-name my-zoom --work-folder runs --non-interactive
```

If you leave out `--center-mode` and `--box-side`, an **interactive** window opens instead. It
plots the particles in two projections (x-y and x-z) with the candidate centres marked and the
current box drawn on top. At the terminal you can switch the centre (`g` geometric, `c` centre of
mass, `d` densest, `u` your own coordinates) and change the box size, and the plot redraws. Press
`a` to accept, or `q` to quit without writing anything. On a machine with no display it writes a
PNG each round instead. The accepted selection is also saved as `<snapshot>_zoom_selection.png`.

Other useful options: `--shift-snapshot` (also write a copy of the snapshot centred on the origin),
`--output` (where to write the config), `--non-interactive` (no prompts, for scripting; give
`--run-name` and `--work-folder`), and `--help` for the full list.

## Tools

`tools/` holds helper scripts for plotting (`plot_cartesian_map.py`, `plot_diagnostic.py`,
`plot_h_comparison.py`), the config builder described above (`suggest_config.py`), and a
lattice regression test (`test_lattice_grid.cpp`, built with the project).

## Tests

From the build directory, after building:

```
ctest --output-on-failure
```

- `lattice_tessellation`: a Voronoi tessellation of a cubic lattice must be exact cubes.
- `synthetic_pipeline`: builds a small synthetic SPH cloud with a known total mass and bulk velocity,
  runs `sph-to-grid` and the full Voronoi pipeline, and checks that mass is conserved and the
  velocity is recovered, for every kernel. It needs Python 3 and nothing else.
- `real_data` (about a minute): `sph-to-grid` on a real Phantom snapshot, `tests/data/brick_03927.ascii`
  (a `splash to ascii` dump of a cubic-spline run; ~350 MB, so it is not tracked by git and the test is
  skipped when the file is absent). Checks that mass and momentum on the grid equal the particle sums
  for the cubic, quintic and Wendland C6 kernels, that the generic integral reproduces the cubic fast
  path, and that the tabulated integral stays close. Run it alone with `ctest -R real_data`.
- `generic_kernel`: the generic kernel integral against the hand-written cubic-spline one and against
  brute-force quadrature, plus whole-cell mass checks, for all six kernels.
- `kernel_tables`, `kernel_spec_data_current`: the kernel-integral coefficients (derived exactly in
  `tools/kernel_tables/`, see its README) against brute-force quadrature, and a check that the
  generated `src/kernel_spec_data.cpp` is up to date.

The test data is generated on the fly, so nothing needs downloading.

## Known limitations

- The POLARIS quantity list is fixed and not configurable.
- `sph-to-grid` has no native POLARIS or ASCII writer yet.
- The Voronoi pipeline (`map`, `cartesian-map`) is only tested on synthetic data; the real-snapshot
  test covers `sph-to-grid`.

## Citation

If you use KIWI, please cite the two papers behind its mapping methods:

- Petkova, M. A., Laibe, G., Bonnell, I. A., 2018, *Fast and accurate Voronoi density gridding from
  Lagrangian hydrodynamics data*, Journal of Computational Physics, 353, 300.
  [doi:10.1016/j.jcp.2017.10.024](https://doi.org/10.1016/j.jcp.2017.10.024)
- Petkova, M. A., Vandenbroucke, B., Bonnell, I. A., Kruijssen, J. M. D., 2021, *Modelling of
  ionizing feedback with smoothed particle hydrodynamics and Monte Carlo radiative transfer on a
  Voronoi grid*, MNRAS, 507, 858.
  [doi:10.1093/mnras/stab2178](https://doi.org/10.1093/mnras/stab2178)

GitHub's "Cite this repository" button uses `CITATION.cff`.

## License

Copyright (C) 2026 Maya A. Petkova.

KIWI is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. See [LICENSE](LICENSE).

KIWI bundles third-party code under its own licenses; see [THIRD_PARTY.md](THIRD_PARTY.md).
