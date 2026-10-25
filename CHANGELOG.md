# Changelog

All notable changes to KIWI are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/) (while the version is 0.x, the interface may change between
releases).

## [Unreleased]

## [0.1.0] - not yet released

First public version (pre-beta).

### Added

- **Command-line stages** driven by one JSON config per run: `point-set`, `tessellate`, `map`,
  `cartesian-map`, `sph-to-grid`, `write-voronoi`, `write-polaris`, `inverse-mapping-test` and `manifest`.
  Each stage caches its output and runs the stages before it if they are missing. `kiwi --version` prints
  the version.
- **Kernel-integral mapping** of SPH density and velocity onto a Voronoi grid (built with the bundled
  Voro++), following Petkova et al. (2018, 2021): the density in each cell is the integral of the particle
  smoothing kernels over the cell, so mass is conserved cell by cell.
- **Six SPH kernels**: cubic, quartic and quintic splines and Wendland C2, C4 and C6, selected with the
  `kernel` config block, through one generic analytic integral. The cubic spline also has a hand-written
  closed-form fast path and an optional lookup-table integral (`mapping.use_tabulated_integral`). The
  coefficients are generated and checked by `tools/kernel_tables/`.
- **Direct SPH-to-Cartesian mapping** (`sph-to-grid`), with per-axis resolution, and resampling of the
  Voronoi result onto a Cartesian grid (`cartesian-map`).
- **Input formats**: Phantom binary dumps read directly (`phantom_binary`), `splash to ascii` output
  (`splash_ascii`) and a headerless 9-column text format (`plain_ascii`), with automatic detection.
- **Buffer points and Lloyd relaxation** for the Voronoi grid: `coverage_grid` (general purpose, the
  default) and `graded_shells` (for one slab-shaped cloud geometry).
- **Provenance manifest** (`manifest.json`) recording the config and files behind each output, so a changed
  config is reported instead of silently reusing stale cached output.
- **Plain-text Voronoi grid writer** (`write-voronoi`) and an **experimental POLARIS grid writer**
  (`write-polaris`: fixed quantity list of gas temperature, velocity and density; the interface may change).
- **Helper scripts** in `tools/`: `suggest_config.py` (propose a config from a raw snapshot, including a
  zoom assistant), `plot_cartesian_map.py` (slices and projections of either Cartesian output, with units)
  and diagnostic plots.
- **Tests and CI**: lattice tessellation, generic-kernel, synthetic end-to-end, Phantom-dump, POLARIS
  writer, bad-input, example-config, version and kernel-table tests, plus a real-data test on a small
  cutout of a Phantom snapshot (`tests/data/`) and an optional full-snapshot test. GitHub Actions builds
  and tests on Ubuntu and macOS, with an OpenMP thread-count matrix and an ASan/UBSan job.
- **Documentation**: a worked example (`configs/example-cutout.json`), `CITATION.cff` (the two method
  papers and the thesis with the extra kernel coefficients) and `CONTRIBUTING.md`.

### Known limitations

- The output file of `sph-to-grid` is a placeholder binary format; there is no native ASCII or POLARIS
  writer for it yet.
- The POLARIS writer's quantity list is fixed.
- `cartesian-map` samples the Voronoi cells at the Cartesian cell centres, so it does not conserve mass
  when the Cartesian cells are larger than the Voronoi cells. The `map` output itself does conserve mass.
- The full-snapshot test covers `sph-to-grid` only; the Voronoi `map` stage is tested on synthetic data and
  on the small real cutout.

[Unreleased]: https://github.com/mapetkova/kiwi/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/mapetkova/kiwi/releases/tag/v0.1.0
