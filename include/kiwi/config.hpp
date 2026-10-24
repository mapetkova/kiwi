// Config schema for the KIWI pipeline.
//
// Design intent: this struct is the SINGLE SOURCE OF TRUTH for a run. Every
// stage reads only from a Config (or a narrower struct built from it) --
// never from a hardcoded path or a compiled-in bool. A Config loaded from
// disk round-trips back to JSON unchanged, so it can be stamped into a
// provenance manifest alongside the run's outputs (see manifest.hpp).
#pragma once

#include <string>
#include <vector>
#include <optional>
#include <stdexcept>

#include "nlohmann/json.hpp"
#include "kiwi/kernel_spec.hpp"

namespace kiwi {

struct ConfigError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

struct Domain {
    double x_min = 0, x_max = 0;
    double y_min = 0, y_max = 0;
    double z_min = 0, z_max = 0;

    void validate() const;
};

// voro++'s own pre_container_base::guess_optimal() sizes its internal
// spatial hash for a DOMAIN-UNIFORM particle density (its optimal_particles
// constant targets ~5.6 points/block on average) -- and completely ignores
// whatever nx/ny/nz values are passed to it (ilscale is derived purely from
// total_particles() and domain volume; see pre_container.cc). KIWI's
// point sets are usually NOT domain-uniform: a real SPH cloud (often
// clumpy/filamentary) plus, when buffer.enabled, a much sparser buffer
// lattice filling the rest of the domain. guess_optimal()'s domain-average
// sizing under-resolves the dense cloud region, where the actual local
// particle density is far above the nominal 5.6/block target.
//
// block_scale_factor multiplies guess_optimal()'s own nx/ny/nz by this
// factor (rounded, clamped to >= 1) before the container is built. 1.0
// keeps guess_optimal()'s own sizing -- block count only changes voro++'s
// internal search granularity, never which cells/volumes/neighbors it
// finds, so this can never change output values, only speed.
//
// Benchmarked (bench_grid_blocks, tools/) on two real point sets -- 140,955
// points (coverage_grid buffer + Lloyd) and 448,619 points (graded_shells buffer +
// Lloyd):
// tessellation time as a function of scale forms a clear U-shape, with
// guess_optimal's own choice (scale=1.0) 2-3.5x SLOWER than the minimum,
// which fell consistently in the 2.5-3.5x-finer range on both datasets
// (going further, e.g. 6x, starts getting slower again -- too many
// near-empty blocks to search). Default here (3.0) sits at or near that
// minimum on both; re-run bench_grid_blocks against your own data if it
// looks very different (e.g. much less clumpy) to check whether a different
// value serves it better.
struct GridBlocks {
    double block_scale_factor = 3.0;
};

// The graded-shells buffer strategy places three nested shells of synthetic
// points around the faces of the domain, with spacing growing outward
// (shell_spacing_1/2/3), offset from a smaller inner region defined by
// (dx, dy, dz). It was written for one specific slab-shaped cloud geometry
// and is NOT a generic strategy -- see BufferConfig::strategy for the general
// alternative.
struct GradedShellsBufferParams {
    double dx = 0, dy = 0, dz = 0;
    double shell_spacing_1 = 3.0; // first two shells use this spacing
    double shell_spacing_2 = 4.0;
    double shell_spacing_3 = 5.0;
};

// The generic buffer strategy: lay a Cartesian candidate grid across the
// whole domain (spacing = spacing_fraction * the domain's smallest extent,
// cell-centered so no candidate sits exactly on the box boundary), then keep
// only the candidates NOT already covered by any SPH particle's smoothing
// sphere (radius = coverage_radius_factor * that particle's h). Applies
// domain-wide, not just near the walls -- interior voids in an irregular
// cloud get buffer cells too, not only the outskirts.
struct CoverageGridBufferParams {
    double spacing_fraction = 0.05;      // dx=dy=dz = this fraction of the domain's smallest extent
    double coverage_radius_factor = 2.0; // a candidate is "covered" if within factor*h of any particle
                                         // (2 = the cubic spline's support; raise it for wider kernels)
};

struct BufferConfig {
    bool enabled = false;
    // "graded_shells": three nested shells of growing spacing, for one specific slab-shaped cloud geometry (needs graded_shells params)
    // "coverage_grid": generic particle-coverage-based buffer (needs coverage_grid params)
    std::string strategy = "coverage_grid";
    GradedShellsBufferParams graded_shells;
    CoverageGridBufferParams coverage_grid;
};

struct LloydConfig {
    int iterations = 0;
    bool save_each_iteration = false;
};

// Optional visual sanity-check plots (see diagnostics.hpp). Off by default --
// when enabled, the point-set stage writes an x-y scatter plot right after
// the buffer stage (if buffer.enabled) and right after Lloyd relaxation
// finishes (if lloyd.iterations > 0), so buffer placement and the effect of
// relaxation can be eyeballed without opening points-buff.dat/points-final.dat
// by hand. Best-effort: requires python3 + matplotlib, but a missing
// toolchain only prints a warning, never fails the run.
struct DiagnosticsConfig {
    bool enabled = false;
    // Which python3 to invoke. Defaults to "python3" (resolved via PATH at
    // run time, same as before this field existed). Override with a full
    // interpreter path (e.g. "/usr/bin/python3" or a venv's
    // ".../bin/python3") if `pip3`/`python3` on PATH turn out to be two
    // different installs -- a real failure mode on macOS, where several
    // Python installs (system, Homebrew, a venv, ...) can coexist and a
    // package installed via one isn't visible to another.
    std::string python_executable = "python3";
};

// Stage 4: density + mass-weighted velocity mapping. Deliberately minimal:
// there is no smoothing step to configure.
//
// use_tabulated_integral selects which form of the SPH kernel line integral
// (see kernel_integral.hpp) computes density/velocity: false (the default)
// uses full_integral_new(), the analytic closed-form "ground truth"; true
// uses gridded_integral(), a fast trilinear-interpolated approximation of
// it, normally reading its lookup table from density_table (a precomputed
// table file). If density_table is
// left unset, or points at a file that isn't actually there,
// run_mapping_stage() generates the same canonical table in memory instead
// of failing (see generate_density_table() in kernel_integral.hpp/cpp) --
// so this field is a preference, not a hard requirement, when
// use_tabulated_integral is true.
struct MappingConfig {
    bool use_tabulated_integral = false;
    std::string density_table; // optional: loaded if set and present, else auto-generated in memory
};

// Stage 5a: POLARIS binary grid writer (EXPERIMENTAL: fixed quantity list, interface may change).
// Writes POLARIS's Voronoi grid file (header: sim/params/plist/cells_d/len;
// per cell: position, volume, [gas temperature, vx, vy, vz, density],
// neighbor list -- see writer_polaris.cpp for the byte layout). The file
// stores SI quantities, so unit-conversion scale factors and a uniform gas
// temperature are config; the defaults convert pc, km/s and Msun/pc^3 to SI
// and use 65 K.
//
// sim_id/plist (the 5 POLARIS GRID data-quantity IDs, in fixed order:
// [gas temperature, vx, vy, vz, gas density]) are NOT configurable yet --
// this writer targets exactly the quantity set one particular POLARIS setup expects.
// A fully generic quantity list is natural future work once there's a
// second use case (different species, different fixed fields) to design it
// against.
struct PolarisWriterConfig {
    unsigned short sim_id = 50;
    double length_scale_m = 3.08567758e16;    // meters per input coordinate unit (default: 1 pc)
    double velocity_scale = 1000.0;           // (m/s) per input velocity unit (default: km/s -> m/s)
    double density_scale = 1.9891e30 / (3.08567758e16 * 3.08567758e16 * 3.08567758e16);
                                               // (kg/m^3) per input density unit (default: Msun/pc^3 -> kg/m^3)
    double gas_temperature_k = 65.0;          // uniform gas/dust temperature written to every cell (Kelvin)
};

// Stage: resample the already-computed Voronoi grid (density/velocity per
// cell, from the mapping stage) onto a uniform Cartesian grid -- a quick,
// simple-to-visualize output, NOT part of the pipeline's data contract. See
// cartesian_map.hpp for the algorithm: for each Cartesian grid-cell center,
// find the enclosing Voronoi cell via voro++'s own nearest-site search and
// copy that cell's density/velocity across.
//
// resolution is cells per axis (the domain must be a cube, same requirement
// write_polaris_grid() already has). The default is 128; it is a config
// field so that a quick look and a final high-resolution pass are both
// cheap to ask for (500^3 = 125M lookups is a lot).
struct CartesianMapConfig {
    int resolution = 128;
};

// Stage: sph-to-grid -- a wholly independent SPH -> Cartesian grid
// mapping, computed directly (no voro++, no buffer, no Lloyd, no
// tessellation) rather than resampling the Voronoi pipeline's output like
// CartesianMapConfig/cartesian_map.hpp does. See cartesian_grid.hpp for the
// algorithm.
//
// Deliberately a separate, self-contained knob set from MappingConfig
// rather than reusing it -- these two mapping paths are meant to stay
// wholly independent of one another (sph-to-grid is a separate path through
// the tool, without buffer or Lloyd iterations), so a
// change to one's config shouldn't have to reason about the other.
//
// resolution: cells per axis, used for all 3 axes UNLESS overridden below --
// on a cube domain this gives cube cells directly, same as
// CartesianMapConfig::resolution.
//
// resolution_x/resolution_y/resolution_z: optional per-axis overrides (0 =
// unset, falls back to `resolution` for that axis). They let a non-cube
// domain be given cube -- or at least deliberately chosen -- cells by
// picking per-axis counts that match the domain's own aspect ratio, and let
// the grid line up EXACTLY, cell-for-cell, with an external tool's
// differently-shaped grid (e.g. SPLASH's `splash to grid`, which routinely
// picks a different nx/ny/nz per axis to keep voxels roughly cubic on a
// non-cube domain) for a direct cross-check. run_sph_to_grid_stage()'s
// cell-shape warning checks the actual resolved cell dimensions (dpx/dpy/dpz)
// rather than the domain's own shape, since per-axis resolution can make
// cells cubes even on a non-cube domain.
//
// use_tabulated_integral / density_table: same semantics as MappingConfig's
// fields of the same name (see kernel_integral.hpp) -- false (the default)
// uses the analytic full_integral_new()-derived tiers directly; true uses
// the tabulated fast-path approximation.
struct CartesianGridConfig {
    int resolution = 128;
    int resolution_x = 0;
    int resolution_y = 0;
    int resolution_z = 0;
    bool use_tabulated_integral = false;
    std::string density_table;
};

// Which SPH kernel the snapshot was run with, and therefore which kernel the grid mapping integrates.
// Everything here is optional; the defaults reproduce the cubic spline (M4) exactly as before.
//
// name: "cubic_spline" (default), "quartic_spline", "quintic_spline", "wendland_c2", "wendland_c4",
//   "wendland_c6" (see kernel_spec.hpp for aliases).
// support_radius_in_h: the kernel's support radius in units of the snapshot's smoothing length h; 0
//   (default) means the SPLASH/Phantom convention for the chosen kernel (2 for the cubic spline and
//   the Wendland kernels, 2.5 for the quartic, 3 for the quintic). Set it if the snapshot's h follows
//   another convention (e.g. a code whose kernel vanishes at r = h would need 1).
// hfact: only used when a splash_ascii snapshot lacks either h or density and one is derived from the
//   other via h = hfact * (m/rho)^(1/3); 0 (default) means Phantom's default for the chosen kernel.
// generic_path: false (default) runs the cubic spline through its original hand-written, optionally
//   tabulated integral and every other kernel through the generic analytic one; true forces the
//   generic analytic path even for the cubic spline (mainly for testing; incompatible with
//   use_tabulated_integral, which only exists for the cubic spline).
struct KernelConfig {
    std::string name = "cubic_spline";
    double support_radius_in_h = 0.0;
    double hfact = 0.0;
    bool generic_path = false;

    bool is_default() const {
        return name == "cubic_spline" && support_radius_in_h == 0.0 && hfact == 0.0 && !generic_path;
    }
};

// Units a "phantom_binary" snapshot is converted to when read (see phantom_dump.hpp): the dump stores
// code units plus the cgs factors udist (length), umass (mass) and utime (time). The defaults are the
// conventions of SPLASH's `to ascii` export (1 pc = 3.086e18 cm, 1 M_Sun = 1.989e33 g, 1 km/s = 1e5 cm/s),
// so a dump read directly agrees with its own SPLASH ascii export. Ignored for ascii snapshots.
struct PhantomUnits {
    double length_cm = 3.086e18;
    double mass_g = 1.989e33;
    double velocity_cms = 1.0e5;

    bool is_default() const {
        return length_cm == 3.086e18 && mass_g == 1.989e33 && velocity_cms == 1.0e5;
    }
};

struct Config {
    std::string run_name;
    Domain domain;
    GridBlocks grid_blocks;
    std::string input_snapshot; // raw SPH snapshot file -- see input_snapshot_format
    // "auto" (default): detect between the formats below by sniffing the
    // file (see detect_snapshot_format() in snapshot.hpp). "phantom_binary":
    // Phantom's native binary dump, read directly (see phantom_dump.hpp;
    // converted to the units in phantom_units). "plain_ascii":
    // headerless text, 9 columns per particle (x y z vx vy vz m h rho). "splash_ascii": output of `splash to ascii` (SPLASH repo) --
    // commented header, then per-particle data whose column count/order is
    // read from the header line itself (not hardcoded -- it varies between
    // files depending on the source Phantom run), non-gas (itype != 1) rows
    // dropped when an itype column is present. See snapshot.hpp for the
    // full format documentation.
    std::string input_snapshot_format = "auto";
    std::string work_folder;    // base output directory
    BufferConfig buffer;
    LloydConfig lloyd;
    MappingConfig mapping;
    PolarisWriterConfig polaris_writer;
    DiagnosticsConfig diagnostics;
    CartesianMapConfig cartesian_map;
    CartesianGridConfig cartesian_grid;
    KernelConfig kernel;
    PhantomUnits phantom_units;
    // phantom_binary only: accept a small dump that has no velocity arrays, setting all velocities to
    // zero (density is unaffected, velocity outputs are meaningless). See phantom_dump.hpp.
    bool allow_missing_velocity = false;

    // Resolved kernel settings (throw for an unknown kernel name).
    const KernelSpec& kernel_spec() const;
    double kernel_support_in_h() const; // support radius / h
    double kernel_hfact() const;
    // True when the generic analytic integral must be used instead of the cubic-spline fast path.
    bool use_generic_kernel_path() const;

    static Config load_from_file(const std::string& path);
    nlohmann::json to_json() const;

    // <work_folder>/<basename of input_snapshot, extension stripped>/ --
    // every stage of a run writes its outputs here. Created on first use.
    std::string run_dir() const;
    // run_dir() + "/" + filename -- the single place stage output paths are
    // computed, so no stage hardcodes its own naming scheme.
    std::string stage_path(const std::string& filename) const;
};

} // namespace kiwi
