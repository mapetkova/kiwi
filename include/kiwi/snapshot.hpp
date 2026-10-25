// Raw SPH snapshot ingestion. This is the entry point of the pipeline:
// everything downstream (buffer, Lloyd, tessellation, and the
// density/velocity mapping stage) starts from what's read here.
//
// Three file formats are recognized (the third, "phantom_binary", is Phantom's native binary dump, read
// directly -- see phantom_dump.hpp):
//
//   "plain_ascii": a headerless text file (the format of the early KIWI runs, not
//   tied to any simulation code), 9 whitespace-separated columns per particle --
//     x  y  z  vx  vy  vz  m  h  rho
//   (m is constant across particles for an equal-mass SPH run; h is the
//   smoothing length; rho is density.) Particle id is the 1-based row
//   number -- there is no explicit id column, and the per-cell stage outputs
//   are aligned to this same row order, so that convention is kept rather
//   than inventing IDs.
//
//   "splash_ascii": output of `splash to ascii` (see the SPLASH repo,
//   github.com/danieljprice/splash). A commented preamble (a signature
//   line naming SPLASH, then `time`/`npart`/`units` blocks, each a `#`
//   comment line followed by one uncommented line of values, then a final
//   `#` column-header line) precedes the per-particle data. The number,
//   order, and identity of the data columns are NOT hardcoded: `splash to
//   ascii` writes whatever columns the source Phantom dump actually has,
//   in whatever order (see write_sphdump(), SPLASH's write_sphdata.f90 --
//   both ncolumns and the column order are dynamic there, driven by which
//   physics modules a given run had switched on -- e.g. Bx/By/Bz only for
//   an MHD run), and even whether an `itype` column is present at all
//   depends on whether the dump has more than one particle type. So this
//   reader parses the column-header comment line itself (fixed-width
//   fields -- see split_splash_ascii_labels() in snapshot.cpp for why a
//   plain whitespace split doesn't work) into a name -> column-index map,
//   and looks columns up by name rather than assuming a fixed layout.
//   Only x, y, z, m (particle mass), vx, vy, vz, and at least one of h
//   (smoothing length) or density are required; anything else present
//   (pressure, alpha, div_v, poten, dt, or an unrecognized physics-module
//   column) is read past but not kept -- `pressure` in particular is not
//   needed by anything downstream, even though real files carry it. If only
//   one of h/density is present, the other is derived via
//   h = hfact*(m/rho)^(1/3) (Price et al. 2018, eq. 6; hfact = 1.2 for the
//   cubic spline). The relation needs m, h and rho in one consistent unit
//   system, but `splash to ascii` writes them in three different ones
//   (M_Sun, pc, g/cm^3), so the values are converted through Msun/pc^3
//   first; if the header's own unit brackets say anything else, the reader
//   refuses with a clear error rather than silently misconverting.
//   If no `itype` column is present at all, every row is treated as gas
//   (matching SPLASH's own single-particle-type convention); otherwise
//   only itype == 1 (gas) rows are kept, and sink/dust/star/etc. rows are
//   dropped with a warning rather than silently treated as gas. Particle
//   id is the 1-based row number among the *kept* rows.
//
//   Column identity is always read from the header, never assumed from a
//   fixed layout: a hard-coded column order silently misreads (for example
//   comparing `itype` against the wrong column and dropping every particle
//   as "non-gas") as soon as a file has a column more or fewer than expected.
//
// "auto" (the default -- see Config::input_snapshot_format) picks the format
// from the file itself: a Phantom capture-pattern record means
// phantom_binary; otherwise the first line decides -- the SPLASH signature
// comment means splash_ascii, anything else is assumed to be plain_ascii.
#pragma once

#include <string>
#include <vector>

#include "kiwi/config.hpp"
#include "kiwi/point_set.hpp"

namespace kiwi {

struct SPHSnapshot {
    std::vector<int> id; // 1-based row number (among rows actually kept)
    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<double> m, h, rho;

    size_t size() const { return id.size(); }
};

// Sniffs `path` to tell the formats apart: a Phantom binary dump is recognised by its capture-pattern
// record, otherwise the first line decides between plain_ascii and splash_ascii.
// Returns "phantom_binary", "plain_ascii" or "splash_ascii" -- never "auto". Exposed
// separately from load_sph_snapshot() so callers (suggest_config.py has a cheaper
// version of this rule) can report or act on the detected format
// before/without loading the whole file.
std::string detect_snapshot_format(const std::string& path);

// `format` is one of "auto" (default -- resolved via detect_snapshot_format),
// "plain_ascii", "splash_ascii", or "phantom_binary". Throws std::runtime_error for any other
// value, an unreadable file, or a file that doesn't parse as the resolved
// format.
//
// `hfact` is the factor in h = hfact * (m/rho)^(1/3), used only for a splash_ascii snapshot that has a
// density or an h column but not both (default 1.2, the Phantom cubic-spline value; Config::kernel_hfact()
// gives the per-kernel value).
//
// `units` and `allow_missing_velocity` are only used for "phantom_binary" snapshots (see phantom_dump.hpp).
SPHSnapshot load_sph_snapshot(const std::string& path, const std::string& format = "auto", double hfact = 1.2,
                              const PhantomUnits& units = PhantomUnits(), bool allow_missing_velocity = false);

// Writes `snap` as a plain_ascii file (headerless, 9 columns per particle: x y z vx vy vz m h rho, full
// double precision so reading it back gives the same values). Throws std::runtime_error if the file can't
// be written. The rho column is written as loaded: g/cm^3 for Phantom and splash_ascii input, as read for
// plain_ascii input. Used by `kiwi convert-snapshot`.
void write_plain_ascii(const SPHSnapshot& snap, const std::string& path);

// Just the positions, as a PointSet ready for the buffer/Lloyd/tessellation
// stages (which don't need to know anything about velocity/mass/density).
PointSet extract_points(const SPHSnapshot& snap);

} // namespace kiwi
