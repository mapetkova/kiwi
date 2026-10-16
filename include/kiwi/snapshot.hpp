// Raw SPH snapshot ingestion. This is the real entry point of the pipeline:
// everything downstream (buffer, Lloyd, tessellation, and eventually the
// density/velocity mapping stage) starts from what's read here.
//
// Two file formats are recognized:
//
//   "dale_phantom": the original Dale-et-al.-2019-style Phantom ASCII
//   dump, no header line, 9 whitespace-separated columns per particle --
//     x  y  z  vx  vy  vz  m  h  rho
//   (m is constant across particles for an equal-mass SPH run; h is the
//   smoothing length; rho is density.) Particle id is the 1-based row
//   number -- there is no explicit id column, and every other file in the
//   legacy pipeline (tabulated density/velocity) is implicitly aligned to
//   this same row order, so we preserve that convention rather than
//   inventing IDs.
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
//   column) is read past but not kept -- `pressure` in particular is NOT
//   needed by anything downstream, even though real files from the
//   original pipeline carry it. If only one of h/density is present, the other is
//   derived via h = 1.2*(m/rho)^(1/3) (Price et al. 2018, eq. 6; the sign
//   matters -- (m/rho)^(-1/3) would be backwards), converting through
//   Msun/pc^3 first since m/h/density are written in three different unit
//   systems (M_Sun/pc/g_cm^3) that the relation can't be applied to
//   directly -- this refuses with a clear error, rather than silently
//   misconverting, if the header's own unit brackets say anything else.
//   If no `itype` column is present at all, every row is treated as gas
//   (matching SPLASH's own single-particle-type convention); otherwise
//   only itype == 1 (gas) rows are kept, and sink/dust/star/etc. rows are
//   dropped with a warning rather than silently treated as gas. Particle
//   id is the 1-based row number among the *kept* rows.
//
//   (An earlier version of this reader assumed a fixed 15-column layout --
//   x y z m h rho vx vy vz pressure alpha div_v poten dt itype -- and an
//   even earlier version of that was missing the `pressure` column
//   entirely, which silently shifted every field from alpha onward by
//   one, so itype ended up compared against what was actually the `dt`
//   column and every real gas particle got dropped as "non-gas": see the
//   "Nineteenth delivery" section of the project progress doc. Both were
//   real bugs; reading column identity from the header instead of
//   assuming it fixes the whole class, not just that one instance.)
//
// "auto" (the default -- see Config::input_snapshot_format) picks between
// the two by sniffing the file's first line for the SPLASH signature
// comment; a file with no such first line is assumed to be dale_phantom.
#pragma once

#include <string>
#include <vector>

#include "kiwi/point_set.hpp"

namespace kiwi {

struct SPHSnapshot {
    std::vector<int> id; // 1-based row number (among rows actually kept)
    std::vector<double> x, y, z;
    std::vector<double> vx, vy, vz;
    std::vector<double> m, h, rho;

    size_t size() const { return id.size(); }
};

// Sniffs `path`'s first line to tell dale_phantom and splash_ascii apart.
// Returns "dale_phantom" or "splash_ascii" -- never "auto". Exposed
// separately from load_sph_snapshot() so callers (and suggest_config.py's
// Python port of this same rule) can report or act on the detected format
// before/without loading the whole file.
std::string detect_snapshot_format(const std::string& path);

// `format` is one of "auto" (default -- resolved via detect_snapshot_format),
// "dale_phantom", or "splash_ascii". Throws std::runtime_error for any other
// value, an unreadable file, or a file that doesn't parse as the resolved
// format.
SPHSnapshot load_sph_snapshot(const std::string& path, const std::string& format = "auto");

// Just the positions, as a PointSet ready for the buffer/Lloyd/tessellation
// stages (which don't need to know anything about velocity/mass/density).
PointSet extract_points(const SPHSnapshot& snap);

} // namespace kiwi
