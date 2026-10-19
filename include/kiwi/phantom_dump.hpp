// Direct reader for Phantom's native binary dump files ("phantom_binary" snapshot format), so a dump
// can be used without first converting it with `splash to ascii` (a ~4x larger text file).
//
// Format (documented at https://phantomsph.readthedocs.io/en/latest/user-guide/dumpfile.html):
// Fortran unformatted sequential records (each wrapped by a 4-byte length marker), in this order:
//   1. a capture-pattern record  i1 r1 i2 iversion i3  (60769, 60878.0, 60878, version, 690706) whose
//      size and values reveal the byte order and the widths of the default integer and real;
//   2. a 100-character file identifier ("FT:Phantom:..." = full dump, tagged format);
//   3. the global header: eight typed sections (default int, int*1, int*2, int*4, int*8, default real,
//      real*4, real*8), each = int32 count, 16-character tags, values;
//   4. int32 nblocks, then per block int64 n + 8 x int32 array counts (one per type);
//   5. per block, for each type in order, its arrays as 16-character tag + n values.
// Non-MPI dumps have two blocks: the gas (and other non-sink) particles, then the sink particles; MPI
// dumps repeat that pair once per MPI task (header "nblocks").
//
// What KIWI takes from it (everything is looked up by tag, so extra arrays are skipped unread):
//   x y z vx vy vz h   particle arrays (real*4 or real*8);   itype (any int width) if present
//   massoftype, hfact, udist, umass, utime   from the header
// Phantom stores neither particle mass nor density: the mass of a type-1 (gas) particle is
// massoftype(1) and the density is massoftype(1) * (hfact / h)^3. Particles with h <= 0 are dead
// (accreted onto sinks) and are dropped, as are particles with itype != 1 and the sink block.
//
// Written from the format description above and checked against a real dump and its SPLASH ascii
// export; Sarracen's `read_phantom` (Harris & Tricco 2023, JOSS 8(86), 5263) served as a cross-reference
// for the layout -- no code was copied from it.
#pragma once

#include <string>

#include "kiwi/snapshot.hpp"

namespace kiwi {

// True iff `path` starts with a Phantom capture-pattern record (any byte order / int / real width).
// Never throws; false for unreadable or non-Phantom files.
bool is_phantom_dump(const std::string& path);

// Loads the gas particles of a Phantom dump, converting code units to the units in `units`
// (default: pc, M_Sun, km/s) using the dump's own udist/umass/utime. Throws std::runtime_error with a
// descriptive message for anything it does not understand (not a Phantom dump, untagged or small dump
// missing required arrays, APR levels, types it cannot tell apart, ...).
SPHSnapshot load_phantom_dump(const std::string& path, const PhantomUnits& units);

} // namespace kiwi
