#!/usr/bin/env python3
"""Interactive config generator: inspects a raw SPH snapshot and writes a
ready-to-run kiwi config.json, so a new run doesn't require hand-editing the
domain box or copying an old config as a template.

Usage:
    suggest_config.py <snapshot_path> [options]

Two domain modes (--mode, default "full"):

  full (default): a cubic domain box (equal x/y/z extent -- see below for
    why) sized to cover a percentile-trimmed particle spread and centered on
    that trimmed spread's own geometric center. Trimming guards against a
    handful of far-flung outlier particles blowing up the box (the same
    class of large-h/far-flung-particle issue documented for the
    coverage_grid buffer strategy in point_set.cpp) -- see --trim-percentile
    below. How many particles fall outside the suggested domain as a result
    is always reported, so nothing is silently dropped without being said.

  zoom: a smaller cubic sub-domain around a chosen center point, for
    focusing a run on one region of a larger snapshot. The center can be:
      1. geometric  -- the bounding-box midpoint of every particle
      2. com        -- the mass-weighted center of mass
      3. densest     -- the density-peak location, found via an iterative
                         coarse-to-fine sliding-box search (see
                         --densest-precision-fraction below)
      4. user        -- exact coordinates you supply (--center X Y Z)
    Pick a sub-domain either manually (--center-mode + --box-side, fully
    scriptable, no prompts) or via the interactive assistance: a live
    matplotlib window plots every particle (x-y and x-z projections) with
    all four candidate centers marked (a persistent marker for 'user' is
    only shown once you've entered coordinates), plus the current
    candidate box as a red square on both panels. Each round you type a new
    center choice and/or box side length at the terminal and the plot
    redraws in place; repeat until you accept. Launched automatically in
    zoom mode whenever --center-mode/--box-side aren't both given (or
    always, if --zoom-interactive is passed); requires matplotlib and a
    terminal, so it's unavailable with --non-interactive. If your
    matplotlib backend has no live display, it falls back to writing a PNG
    after each round instead so you still have something to look at -- the
    terminal prompts work the same either way. The final accepted
    selection is always saved as a PNG too, as a record of what was chosen
    (same idea as the buffer/Lloyd diagnostic plots).

What it asks you in both modes (kept minimal -- everything else defaults to
the tool's known-good settings: coverage_grid buffer, 0 Lloyd iterations,
analytic kernel integral):
    - run_name (defaults to the snapshot's own basename; confirm or override)
    - work_folder (no sensible default -- always asked, unless --work-folder
      is passed)
    - whether to use the tabulated kernel integral instead of the analytic
      default, and if so, the density-table path (unless --density-table is
      passed, which implies yes)
Any of these can be supplied as flags instead, in which case that particular
prompt is skipped -- pass all three (--run-name, --work-folder, and either
nothing or --density-table), plus (for --mode zoom) --center-mode and
--box-side, for a fully non-interactive run.

Why the domain is always a cube: POLARIS's own Voronoi grid loader
(GridVoronoi::loadGridFromBinaryFile in the POLARIS source) stores a single
scalar `max_len` in the binary header and treats the domain as a cube
spanning [-0.5*max_len, +0.5*max_len] on every axis -- this is baked into
POLARIS itself, not a convention of ours. write_polaris_grid() (see
writer_polaris.cpp) already re-centers coordinates at the domain's own
midpoint when it writes the POLARIS file, so the *snapshot* never needs to be
physically shifted for POLARIS's sake -- only the *domain* needs equal
extent in x/y/z (true of both --mode full and --mode zoom's sub-domain
alike). (A real non-cubic domain currently makes write_polaris_grid emit a
loud warning rather than a silently wrong grid -- see its "domain is not a
cube" check.)

--center-mode densest (--mode zoom only) finds the density peak with an iterative
coarse-to-fine box search: start with a cubic box spanning every particle, repeatedly
halve the box size, slide it across the current box in fixed steps of size
`dx = --densest-precision-fraction * the starting box's side` (default 1/128), keep
whichever position holds the most particles, and recenter the next (finer) round there --
until the box shrinks to `dx`, at which point its center is the answer.

--shift-snapshot is a separate, NOT-default option: pass it if you want a
physically CoM-shifted copy of the raw snapshot written to disk (positions
only; vx/vy/vz/m/h/rho columns are untouched) and the generated config
pointed at that copy instead, with the domain re-expressed centered at the
origin. Works the same in either mode (shifts by whatever the final chosen
domain's own center turned out to be). Useful for other tooling/analysis
that expects origin-centered coordinates, but not required for a correct
POLARIS grid.

Requires numpy always; matplotlib only if the interactive zoom assistance
actually runs (numpy is the same dependency the earlier cube-extraction test
scripts in this project already used).

Snapshot format: two raw SPH snapshot formats are recognized, matching
load_sph_snapshot() in snapshot.cpp exactly -- "plain_ascii" (a
headerless text file, 9 columns per particle: x y z vx vy vz m h rho) and
"splash_ascii" (output of `splash to ascii`, from the SPLASH repo:
github.com/danieljprice/splash -- a commented header/preamble followed by
15-column rows: x y z m h rho vx vy vz pressure alpha div_v poten dt itype;
only itype == 1 gas particles are kept, others are dropped with a note
printed).
Auto-detected by default (sniffs the file for the SPLASH signature comment);
override with --snapshot-format if detection ever guesses wrong.
"""
import argparse
import json
import os
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def detect_snapshot_format(path):
    """Mirrors detect_snapshot_format() in snapshot.cpp exactly: sniffs the
    file's first line for the `splash to ascii` signature comment; anything
    else is assumed to be the plain_ascii format (no header at
    all). Returns "plain_ascii" or "splash_ascii" -- never "auto"."""
    with open(path) as f:
        first_line = f.readline()
    if first_line.startswith("#") and "SPLASH" in first_line:
        return "splash_ascii"
    return "plain_ascii"


def _splash_ascii_header_line_count(path):
    """Number of lines to skip before the per-particle data block starts, in
    a `splash to ascii` file. Uses the same rule as
    find_splash_ascii_data_start() in snapshot.cpp: scan forward skipping
    '#' comment lines; whenever an uncommented line is found, peek at the
    line after it -- if that one is ALSO uncommented, the first line is the
    start of the data block (every following row is uncommented too);
    otherwise the first line was a single metadata values line (time/npart/
    units) sandwiched between comments, so skip it and keep scanning. This
    doesn't hardcode how many metadata blocks precede the data, or how many
    columns any of them has (npart's column count in particular varies with
    how many particle types SPLASH knows about)."""

    def is_comment_or_blank(s):
        s = s.strip()
        return s == "" or s.startswith("#")

    with open(path) as f:
        lines = f.readlines()

    i = 0
    n = len(lines)
    while i < n:
        if is_comment_or_blank(lines[i]):
            i += 1
            continue
        if i + 1 < n and not is_comment_or_blank(lines[i + 1]):
            return i
        i += 1
    raise ValueError(f"{path}: could not locate the start of particle data in the "
                      f"splash-ascii header preamble")


def _normalize_splash_label(raw):
    """Strips a bracketed unit suffix (e.g. " [pc]", " [M_{Sun}]") and
    collapses whitespace/underscores/backslashes, so "particle mass
    [M_{Sun}]", "v_x [km/s]", and "\\alpha [1/s]" all normalize down to
    "particlemass", "vx", and "alpha" -- matching normalize_splash_label()
    in snapshot.cpp exactly."""
    out = []
    in_brackets = False
    for c in raw:
        if c == "[":
            in_brackets = True
            continue
        if c == "]":
            in_brackets = False
            continue
        if in_brackets:
            continue
        if c in "\\_ \t":
            continue
        out.append(c.lower())
    return "".join(out)


# Only the columns this tool actually needs -- or needs to recognize in
# order to treat specially, like itype -- are listed here. Any other
# column present in a given file (pressure, alpha, div_v, poten, dt, or a
# physics-module column this codebase has never heard of) is simply
# ignored rather than rejected. Matches splash_ascii_column_aliases() in
# snapshot.cpp exactly.
_SPLASH_ASCII_COLUMN_ALIASES = {
    "x": "x", "y": "y", "z": "z",
    "particlemass": "m", "mass": "m", "m": "m",
    "h": "h", "smoothinglength": "h",
    "density": "rho", "rho": "rho",
    "vx": "vx", "vy": "vy", "vz": "vz",
    "itype": "itype",
}


def _extract_label_unit(raw):
    """Unit string inside the LAST '[...]' bracket of a raw (un-normalized)
    column label, e.g. "density [g/cm^3]" -> "g/cm^3". "" if there's no
    bracket. Matches extract_label_unit() in snapshot.cpp."""
    close = raw.rfind("]")
    if close == -1:
        return ""
    open_ = raw.rfind("[", 0, close)
    if open_ == -1:
        return ""
    return raw[open_ + 1:close]


def _normalize_unit(raw):
    return "".join(c.lower() for c in raw if c not in "{}_ \t")


def _is_msun_unit(raw):
    return _normalize_unit(raw) == "msun"


def _is_pc_unit(raw):
    return _normalize_unit(raw) == "pc"


def _is_gcm3_unit(raw):
    return _normalize_unit(raw) in ("g/cm^3", "g/cm3", "gcm^-3", "gcm-3")


def _msun_pc3_to_g_cm3():
    """g/cm^3 per Msun/pc^3 -- same constants as this project's own
    SPLASH-comparison scripts (compare_splash.py, investigate_discrepancy.py,
    verify_h_hypothesis.py) and snapshot.cpp's msun_pc3_to_g_cm3(), kept in
    one place per language so they can never quietly drift apart."""
    msun_g = 1.9891e33
    pc_cm = 3.08567758e18
    return msun_g / pc_cm ** 3


def _split_splash_ascii_labels(header_line, n_data_cols, path):
    """Splits a splash_ascii column-header comment line ("# x [pc]  y [pc]
    ... itype") into one label string per data column, in file order.
    Mirrors split_splash_ascii_labels() in snapshot.cpp exactly -- see that
    function's comment for why this needs fixed-width field parsing rather
    than a plain whitespace split (multi-word labels like "particle mass
    [M_{Sun}]" or "div v [erg]" can't otherwise be told apart from their
    neighbors), and why a trailing `itype` column, when present, is a
    special unpadded case rather than a normal fixed-width field.

    `n_data_cols` is the column count already known unambiguously from a
    whitespace split of an actual data row (numbers never contain embedded
    spaces, so that split is never ambiguous)."""
    line = header_line.rstrip("\r\n")
    if len(line) < 2 or line[0] != "#":
        raise ValueError(
            f"splash-ascii snapshot: expected the line immediately before the particle data "
            f"to be the '#'-prefixed column-header line, but it wasn't: {path}")
    body = line[1:]
    if body.startswith(" "):
        body = body[1:]  # the format's `1x` right after `#`

    # SPLASH always writes the literal lowercase token `itype` with no
    # padding after it, iff there's more than one particle type -- see the
    # comment in split_splash_ascii_labels() (snapshot.cpp) for why this
    # suffix check is safe even if trailing whitespace elsewhere on the
    # line was stripped by some intervening tool.
    has_itype_tail = body.endswith("itype")

    n_fixed = n_data_cols - 1 if has_itype_tail else n_data_cols
    if n_fixed <= 0:
        raise ValueError(
            f"splash-ascii snapshot: column-header line doesn't have enough fields for the "
            f"{n_data_cols} column(s) seen in the particle data: {path}")

    fixed_part = body[:-5] if has_itype_tail else body

    # Field width is one of exactly two values fixed by SPLASH's own source
    # (23 or 15 content characters + a 1-space separator, depending on
    # whether SPLASH itself was built in double or single precision) --
    # not re-derived by dividing this line's length, which would silently
    # misalign every field if trailing padding on the very last one had
    # been stripped. See split_splash_ascii_labels() (snapshot.cpp).
    field_w = None
    for candidate in (24, 16):
        if has_itype_tail:
            if len(fixed_part) == n_fixed * candidate:
                field_w = candidate
                break
        else:
            min_len = (n_fixed - 1) * candidate + 1
            max_len = n_fixed * candidate
            if min_len <= len(fixed_part) <= max_len:
                field_w = candidate
                break
    if field_w is None:
        raise ValueError(
            f"splash-ascii snapshot: could not reliably split the column-header line into "
            f"{n_fixed} field(s) -- its length doesn't match either of SPLASH's known label "
            f"field widths (23 or 15 characters): {path}")

    labels = []
    for i in range(n_fixed):
        if not has_itype_tail and i == n_fixed - 1:
            field = fixed_part[i * field_w:]  # rest of the string -- tolerant of stripped
                                               # trailing padding
        else:
            field = fixed_part[i * field_w:(i + 1) * field_w]
        labels.append(field.strip())
    if has_itype_tail:
        labels.append("itype")
    return labels


def _load_plain_ascii(path):
    """Read the 9-column, no-header ASCII format (x y z vx vy vz m h rho) --
    the exact format load_sph_snapshot()'s plain_ascii branch (snapshot.cpp)
    reads. Returns a dict of numpy arrays, one per column."""
    import numpy as np

    data = np.loadtxt(path)
    if data.ndim != 2 or data.shape[1] != 9:
        raise ValueError(
            f"expected 9 whitespace-separated columns (x y z vx vy vz m h rho), "
            f"got shape {data.shape} from {path}"
        )
    cols = ["x", "y", "z", "vx", "vy", "vz", "m", "h", "rho"]
    return {name: data[:, i] for i, name in enumerate(cols)}


def _load_splash_ascii(path):
    """Read a `splash to ascii` file -- mirrors load_splash_ascii() in
    snapshot.cpp (see that function and its helpers for full detail).
    Column identity is read from the '#'-prefixed column-header line
    itself (fixed-width fields; see _split_splash_ascii_labels()) rather
    than assumed by position, since `splash to ascii`'s own column count
    and order both genuinely vary between files -- they depend on which
    physics modules a given Phantom run had switched on (see
    write_sphdump(), SPLASH's own write_sphdata.f90), and even whether an
    itype column is written at all depends on whether the dump has more
    than one particle type. `pressure` (present in some snapshots)
    is recognized as a known label but never required or kept,
    since nothing downstream of this script uses it.

    If either `h` or `density` is missing, the other is derived via
    h = 1.2*(m/rho)^(1/3) (Price et al. 2018, eq. 6 -- note the *positive*
    1/3 exponent: (m/rho)^(-1/3) would make h grow with density, which is
    backwards -- denser gas is more highly resolved, i.e. smaller h),
    converting through Msun/pc^3 first since `splash to ascii` writes mass
    in M_Sun, h in pc, and density in g/cm^3 -- three different unit
    systems that h = hfact*(m/rho)^(1/3) can't be applied to directly. This
    refuses with a clear error, rather than silently misapplying that
    conversion, if the header's own unit brackets say anything other than
    those three units.

    Returns the same dict shape as _load_plain_ascii (x/y/z/vx/vy/vz/m/h/
    rho arrays).

    NOTE: this used to hardcode a fixed 15-column, fixed-order layout (x y
    z m h rho vx vy vz pressure alpha div_v poten dt itype); an earlier
    version of that was even missing the `pressure` column entirely, an
    off-by-one shared independently with this format's C++ reader -- see
    the "Nineteenth delivery" section of the project progress doc for that
    story. Both were real bugs; this rewrite (reading column identity from
    the header rather than assuming it) fixes the whole class, not just
    that one instance."""
    import numpy as np

    with open(path) as f:
        lines = f.readlines()
    skip = _splash_ascii_header_line_count(path)
    if skip == 0:
        raise ValueError(
            f"splash-ascii snapshot: no header preamble found before the particle data -- "
            f"expected at least the column-header comment line: {path}")
    header_line = lines[skip - 1]

    data = np.loadtxt(path, skiprows=skip)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.ndim != 2 or data.shape[0] == 0:
        raise ValueError(f"splash-ascii snapshot: no particle data found: {path}")
    n_data_cols = data.shape[1]
    if n_data_cols < 8:  # sanity floor: fewer than x,y,z,m,vx,vy,vz + (h or rho) can't be real
        raise ValueError(
            f"splash-ascii snapshot: particle data has only {n_data_cols} column(s), too few "
            f"to be a real splash_ascii file: {path}")

    labels = _split_splash_ascii_labels(header_line, n_data_cols, path)

    col_index = {}  # semantic name -> column index
    for i, label in enumerate(labels):
        semantic = _SPLASH_ASCII_COLUMN_ALIASES.get(_normalize_splash_label(label))
        if semantic is None:
            continue  # an unrecognized/unneeded column -- ignore it
        if semantic in col_index:
            raise ValueError(
                f"splash-ascii snapshot: the column-header line appears to name "
                f"'{semantic}' more than once (columns {col_index[semantic]} and {i}) -- "
                f"the header likely failed to parse correctly: {path}")
        col_index[semantic] = i

    for req in ("x", "y", "z", "m", "vx", "vy", "vz"):
        if req not in col_index:
            raise ValueError(
                f"splash-ascii snapshot: could not find a '{req}' column in the "
                f"column-header line: {path}")
    have_h = "h" in col_index
    have_rho = "rho" in col_index
    if not have_h and not have_rho:
        raise ValueError(
            f"splash-ascii snapshot: could not find either an 'h' (smoothing length) or a "
            f"'density' column in the column-header line -- need at least one, to derive the "
            f"other from it: {path}")
    have_itype = "itype" in col_index

    if not have_itype:
        eprint(f"  note: no 'itype' column in the header -- this file has only one particle "
               f"type (SPLASH omits itype entirely in that case); treating every row as gas: "
               f"{path}")

    if have_h != have_rho:
        m_unit = _extract_label_unit(labels[col_index["m"]])
        present_unit = _extract_label_unit(labels[col_index["h" if have_h else "rho"]])
        units_ok = _is_msun_unit(m_unit) and (
            _is_pc_unit(present_unit) if have_h else _is_gcm3_unit(present_unit))
        if not units_ok:
            raise ValueError(
                f"splash-ascii snapshot: only a '{'h' if have_h else 'density'}' column was "
                f"found, and this tool only knows how to derive the other from it when mass "
                f"is in M_Sun, h is in pc, and density is in g/cm^3 (the h = 1.2*(m/rho)^(1/3) "
                f"relation needs m, h, and rho in a mutually consistent unit system) -- got "
                f"mass unit '{m_unit}' and {'h' if have_h else 'density'} unit "
                f"'{present_unit}' instead. Extend this derivation before using it on a file "
                f"in different units: {path}")
        eprint(f"  note: only a '{'h' if have_h else 'density'}' column was found in the "
               f"header -- deriving {'density from h' if have_h else 'h from density'} via "
               f"h = 1.2*(m/rho)^(1/3) (Price et al. 2018, eq. 6), converting through "
               f"Msun/pc^3 since density is in g/cm^3 but mass and h aren't: {path}")

    itype = data[:, col_index["itype"]] if have_itype else np.ones(data.shape[0])
    gas = itype == 1
    n_dropped = int((~gas).sum())
    if n_dropped:
        eprint(f"  note: {n_dropped:,} of {data.shape[0]:,} particle row(s) with itype != 1 "
               f"(non-gas -- sink/dust/star/etc.) in {path} were excluded -- only gas "
               f"particles are used for the domain/config suggestion.")
    data = data[gas]

    m = data[:, col_index["m"]]
    if have_h and have_rho:
        h = data[:, col_index["h"]]
        rho = data[:, col_index["rho"]]
    elif have_h:
        # m is in Msun and h is in pc (checked above), but rho must come out in g/cm^3 to
        # match every other splash_ascii density column -- so derive in Msun/pc^3 first (the
        # unit system the hfact relation actually needs) and convert at the end.
        h = data[:, col_index["h"]]
        ratio = h / 1.2
        safe_ratio = np.where(ratio > 0.0, ratio, 1.0)
        rho_msun_pc3 = np.where(ratio > 0.0, m / safe_ratio ** 3, 0.0)
        rho = rho_msun_pc3 * _msun_pc3_to_g_cm3()
    else:
        # Inverse: density column is in g/cm^3, so convert to Msun/pc^3 before applying
        # h = hfact*(m/rho)^(1/3), so the resulting h comes out in pc as expected.
        rho = data[:, col_index["rho"]]
        rho_msun_pc3 = rho / _msun_pc3_to_g_cm3()
        valid = (m > 0.0) & (rho_msun_pc3 > 0.0)
        safe_rho_msun_pc3 = np.where(rho_msun_pc3 > 0.0, rho_msun_pc3, 1.0)
        h = np.where(valid, 1.2 * np.cbrt(m / safe_rho_msun_pc3), 0.0)

    cols = {
        "x": data[:, col_index["x"]], "y": data[:, col_index["y"]], "z": data[:, col_index["z"]],
        "vx": data[:, col_index["vx"]], "vy": data[:, col_index["vy"]],
        "vz": data[:, col_index["vz"]],
        "m": m, "h": h, "rho": rho,
    }
    return cols


def load_snapshot(path, fmt="auto"):
    """Read a raw SPH snapshot in either recognized format (see the module
    docstring). fmt is "auto" (default -- resolved via
    detect_snapshot_format), "plain_ascii", or "splash_ascii". Returns a
    dict of numpy arrays, one per column (x/y/z/vx/vy/vz/m/h/rho), regardless
    of which format it came from."""
    if fmt == "auto":
        fmt = detect_snapshot_format(path)
    if fmt == "plain_ascii":
        return _load_plain_ascii(path)
    if fmt == "splash_ascii":
        return _load_splash_ascii(path)
    raise ValueError(f"unknown snapshot format: {fmt!r} (expected auto, plain_ascii, or "
                      f"splash_ascii)")


def suggest_domain(snap, trim_percentile, pad_fraction):
    """Percentile-trimmed, cubic domain suggestion (the "full" mode).
    Trimming is per-axis and marginal (each axis trimmed independently, not
    a joint/radial trim) -- simple and matches how the earlier true-density
    test-cube extraction in this project already worked with per-axis
    bounding boxes.

    Returns (domain_dict, info_dict) where info_dict carries the numbers
    worth reporting to the user: raw vs. trimmed extents per axis, the
    center used, and how many particles end up outside the final box.
    """
    import numpy as np

    x, y, z = snap["x"], snap["y"], snap["z"]
    n_total = x.size

    raw_bounds = {
        "x": (float(x.min()), float(x.max())),
        "y": (float(y.min()), float(y.max())),
        "z": (float(z.min()), float(z.max())),
    }

    lo_pct, hi_pct = trim_percentile, 100.0 - trim_percentile
    trimmed_bounds = {}
    for name, arr in (("x", x), ("y", y), ("z", z)):
        lo, hi = np.percentile(arr, [lo_pct, hi_pct])
        trimmed_bounds[name] = (float(lo), float(hi))

    center = {name: (lo + hi) / 2.0 for name, (lo, hi) in trimmed_bounds.items()}
    extents = {name: hi - lo for name, (lo, hi) in trimmed_bounds.items()}
    side = max(extents.values()) * (1.0 + pad_fraction)
    half = side / 2.0

    domain = {
        "x_min": center["x"] - half, "x_max": center["x"] + half,
        "y_min": center["y"] - half, "y_max": center["y"] + half,
        "z_min": center["z"] - half, "z_max": center["z"] + half,
    }

    inside = (
        (x >= domain["x_min"]) & (x <= domain["x_max"]) &
        (y >= domain["y_min"]) & (y <= domain["y_max"]) &
        (z >= domain["z_min"]) & (z <= domain["z_max"])
    )
    n_inside = int(inside.sum())

    info = {
        "n_total": n_total,
        "n_inside": n_inside,
        "n_outside": n_total - n_inside,
        "raw_bounds": raw_bounds,
        "trimmed_bounds": trimmed_bounds,
        "center": center,
        "side": side,
        "trim_percentile": trim_percentile,
        "pad_fraction": pad_fraction,
    }
    return domain, info


def print_domain_summary(info):
    print(f"  {info['n_total']:,} particles read")
    for axis in ("x", "y", "z"):
        rlo, rhi = info["raw_bounds"][axis]
        tlo, thi = info["trimmed_bounds"][axis]
        print(f"  {axis}: raw [{rlo:.6g}, {rhi:.6g}]  "
              f"trimmed (p{info['trim_percentile']:g}-p{100 - info['trim_percentile']:g}) "
              f"[{tlo:.6g}, {thi:.6g}]")
    print(f"  cube side: {info['side']:.6g}  "
          f"(padding {info['pad_fraction'] * 100:.1f}%)")
    print(f"  center: ({info['center']['x']:.6g}, {info['center']['y']:.6g}, "
          f"{info['center']['z']:.6g})")
    if info["n_outside"] > 0:
        pct = 100.0 * info["n_outside"] / info["n_total"]
        print(f"  NOTE: {info['n_outside']:,} particle(s) ({pct:.3f}%) fall outside "
              f"this domain and will be excluded from the run -- these are the "
              f"trimmed outliers. Re-run with a smaller --trim-percentile (or 0) "
              f"if you want the domain to cover every particle instead.")
    else:
        print("  all particles fall within this domain.")


# ---------------------------------------------------------------------------
# --mode zoom: sub-domain center candidates and helpers
# ---------------------------------------------------------------------------

def compute_geometric_center(snap):
    """Midpoint of each axis's raw [min, max] -- the bounding-box center of
    every particle in the snapshot (no trimming). One of the four candidate
    sub-domain centers offered in --mode zoom. (Distinct from --mode full's
    own "center", which is the midpoint of the *trimmed* spread -- there is
    no trimming here, zoom's whole point is picking a deliberate region.)"""
    return {
        "x": (float(snap["x"].min()) + float(snap["x"].max())) / 2.0,
        "y": (float(snap["y"].min()) + float(snap["y"].max())) / 2.0,
        "z": (float(snap["z"].min()) + float(snap["z"].max())) / 2.0,
    }


def compute_center_of_mass(snap):
    """Mass-weighted center of mass, using the snapshot's own 'm' column."""
    import numpy as np

    m = snap["m"]
    total_m = float(np.sum(m))
    return {
        "x": float(np.sum(snap["x"] * m)) / total_m,
        "y": float(np.sum(snap["y"] * m)) / total_m,
        "z": float(np.sum(snap["z"] * m)) / total_m,
    }


def compute_densest_region(snap, precision_fraction=1.0 / 128.0):
    """Iterative sliding-box density-peak search: start with a cubic
    box spanning the whole particle set; repeatedly halve the box size,
    slide it across the *current* box in steps of a fixed precision `dx`
    (= precision_fraction * the starting box's side, held constant across
    every round), find the position with the most particles inside, and
    recenter the search there for the next, finer round -- until the box
    side reaches `dx`. The final box's center is the answer. Coarse-to-fine
    rather than one single fine-grained pass, so each round's search is
    cheap even though the final precision is fine.

    Implemented with a 3D histogram + integral image (summed-area table)
    per round rather than literally sliding a window particle-by-particle:
    equivalent to the brute-force description (same box sizes, same dx
    step, same "most particles wins" rule -- verified in testing against an
    actual brute-force reference on synthetic data), just fast enough to
    run on real snapshots. Each round bins only the particles still inside
    the current search box (so later, much smaller rounds are cheap too),
    builds a summed-area table over that histogram, then reads every
    candidate window's particle count off it in one vectorized array
    expression instead of a python loop over window positions.
    """
    import numpy as np

    x, y, z = snap["x"], snap["y"], snap["z"]
    pts = np.column_stack([x, y, z])
    lo = pts.min(axis=0)
    hi = pts.max(axis=0)

    center = (lo + hi) / 2.0
    side = float(np.max(hi - lo))
    if side <= 0:
        return {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])}

    dx = side * precision_fraction
    if dx <= 0:
        dx = side  # degenerate precision_fraction -- fall back to one round only

    max_rounds = 200  # safety cap; real usage needs ~log2(1/precision_fraction) rounds
    rounds = 0
    while side > dx * (1.0 + 1e-9) and rounds < max_rounds:
        rounds += 1
        box_side = side / 2.0
        half = side / 2.0
        box_lo = center - half
        box_hi = center + half

        mask = np.all((pts >= box_lo) & (pts <= box_hi), axis=1)
        active = pts[mask]
        if active.shape[0] == 0:
            break  # nothing left to refine against -- keep the current center

        n_bins = max(2, int(round(side / dx)))
        box_bins = max(1, min(n_bins - 1, int(round(box_side / dx))))

        edges = [np.linspace(box_lo[d], box_hi[d], n_bins + 1) for d in range(3)]
        hist, _ = np.histogramdd(active, bins=edges)

        # 3D summed-area table (integral image): P[i,j,k] = sum of hist over
        # every bin with index < i,< j,< k. Zero-padded on the low side so the
        # inclusion-exclusion window formula below needs no boundary special-case.
        padded = np.pad(hist, ((1, 0), (1, 0), (1, 0)))
        sat = padded.cumsum(axis=0).cumsum(axis=1).cumsum(axis=2)

        b = box_bins
        window_counts = (
            sat[b:, b:, b:] - sat[:-b, b:, b:] - sat[b:, :-b, b:] - sat[b:, b:, :-b]
            + sat[:-b, :-b, b:] + sat[:-b, b:, :-b] + sat[b:, :-b, :-b]
            - sat[:-b, :-b, :-b]
        )

        i0, j0, k0 = np.unravel_index(int(np.argmax(window_counts)), window_counts.shape)
        best_lo = np.array([edges[0][i0], edges[1][j0], edges[2][k0]])
        center = best_lo + box_side / 2.0
        side = box_side

    return {"x": float(center[0]), "y": float(center[1]), "z": float(center[2])}


def domain_from_center_side(center, side):
    half = side / 2.0
    return {
        "x_min": center["x"] - half, "x_max": center["x"] + half,
        "y_min": center["y"] - half, "y_max": center["y"] + half,
        "z_min": center["z"] - half, "z_max": center["z"] + half,
    }


def count_inside_cube(snap, center, side):
    half = side / 2.0
    x, y, z = snap["x"], snap["y"], snap["z"]
    inside = (
        (x >= center["x"] - half) & (x <= center["x"] + half) &
        (y >= center["y"] - half) & (y <= center["y"] + half) &
        (z >= center["z"] - half) & (z <= center["z"] + half)
    )
    return int(inside.sum()), int(x.size)


def print_zoom_summary(center, center_mode, side, n_inside, n_total):
    print(f"  sub-domain: center=({center['x']:.6g}, {center['y']:.6g}, {center['z']:.6g})  "
          f"[{center_mode}]  side={side:.6g}")
    n_outside = n_total - n_inside
    if n_outside > 0:
        pct = 100.0 * n_outside / n_total
        print(f"  {n_inside:,}/{n_total:,} particles fall inside this sub-domain "
              f"({n_outside:,}, {pct:.3f}%, fall outside and will be excluded from the run).")
    else:
        print(f"  all {n_total:,} particles fall inside this sub-domain.")


def resolve_manual_zoom_center(center_mode, center_arg, geometric_center, com_center,
                                densest_center):
    """Resolves --center-mode (+ --center, for 'user') to a concrete center
    dict. Raises ValueError with a message meant to be printed directly to
    the user on anything unresolvable (missing --center for 'user',
    'densest' before it's implemented, an unknown mode)."""
    if center_mode == "geometric":
        return geometric_center
    if center_mode == "com":
        return com_center
    if center_mode == "densest":
        if densest_center is None:
            raise ValueError(
                "--center-mode densest is not implemented yet (placeholder) -- use "
                "geometric, com, or user (--center X Y Z), or omit --center-mode/--box-side "
                "to use the interactive assistance instead."
            )
        return densest_center
    if center_mode == "user":
        if center_arg is None:
            raise ValueError("--center-mode user requires --center X Y Z")
        return {"x": center_arg[0], "y": center_arg[1], "z": center_arg[2]}
    raise ValueError(f"unknown --center-mode: {center_mode!r} (expected geometric, com, "
                      f"densest, or user)")


def _interactive_backend_available():
    """Best-effort check for whether matplotlib's current backend can
    actually show a live window. Non-interactive backends (Agg above all --
    the common default on a headless machine) can still build figures and
    call savefig() just fine; they just can't show() anything, so
    run_interactive_zoom_assistance() falls back to writing a PNG each
    round when this is False rather than silently showing nothing."""
    import matplotlib

    backend = matplotlib.get_backend().lower()
    return backend not in ("agg", "pdf", "ps", "svg", "svgz", "template", "pgf", "cairo")


def run_interactive_zoom_assistance(snap, geometric_center, com_center, densest_center,
                                     initial_center, initial_center_mode, initial_side,
                                     frame_save_prefix):
    """Live matplotlib assistance for picking a --mode zoom sub-domain: a
    two-panel (x-y, x-z) scatter of every particle, the three currently-
    available center candidates marked, and the current candidate box drawn
    as a red square on both panels. Each round, prompts at the terminal for
    a new center choice and/or box side; Enter keeps the current value.
    Returns (center, center_mode, side) once the user accepts, or exits the
    process (code 1) if they abort with 'q'.

    Falls back to writing a PNG after every round (see
    _interactive_backend_available()) when the matplotlib backend has no
    live display -- the prompt loop itself doesn't care either way, only
    what happens after redraw() runs.
    """
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    interactive_backend = _interactive_backend_available()
    if not interactive_backend:
        eprint(f"  note: matplotlib backend '{matplotlib.get_backend()}' has no live display -- "
               f"the plot window won't actually be visible. Writing a PNG after every round "
               f"instead (open it yourself to see the current box); the terminal prompts below "
               f"still work the same way.")

    x, y, z = snap["x"], snap["y"], snap["z"]
    n_total = x.size

    center = dict(initial_center)
    center_mode = initial_center_mode
    side = initial_side
    # Tracks the last coordinates typed at the 'u' prompt, if any -- shown as its own
    # persistent marker (distinct from the current box) so it stays visible even after
    # switching to compare against geometric/com/densest.
    user_center = dict(initial_center) if initial_center_mode == "user" else None

    fig, (ax_xy, ax_xz) = plt.subplots(1, 2, figsize=(13, 6.5))
    frame_n = 0

    def redraw():
        nonlocal frame_n
        for ax in (ax_xy, ax_xz):
            ax.clear()
        ax_xy.scatter(x, y, s=1, c="0.4", label=f"particles ({n_total:,})")
        ax_xz.scatter(x, z, s=1, c="0.4")

        def mark(ax, cx, cy, marker, color, label):
            ax.plot([cx], [cy], marker=marker, ms=13, mew=2, mfc="none", mec=color,
                    ls="none", label=label)

        mark(ax_xy, geometric_center["x"], geometric_center["y"], "+", "tab:blue", "geometric center")
        mark(ax_xz, geometric_center["x"], geometric_center["z"], "+", "tab:blue", None)
        mark(ax_xy, com_center["x"], com_center["y"], "x", "tab:orange", "center of mass")
        mark(ax_xz, com_center["x"], com_center["z"], "x", "tab:orange", None)
        if densest_center is not None:
            mark(ax_xy, densest_center["x"], densest_center["y"], "*", "tab:green", "densest region")
            mark(ax_xz, densest_center["x"], densest_center["z"], "*", "tab:green", None)
        if user_center is not None:
            mark(ax_xy, user_center["x"], user_center["y"], "D", "tab:purple", "user-specified center")
            mark(ax_xz, user_center["x"], user_center["z"], "D", "tab:purple", None)

        half = side / 2.0
        ax_xy.add_patch(Rectangle((center["x"] - half, center["y"] - half), side, side,
                                   fill=False, edgecolor="crimson", lw=2,
                                   label=f"current box ({center_mode})"))
        ax_xz.add_patch(Rectangle((center["x"] - half, center["z"] - half), side, side,
                                   fill=False, edgecolor="crimson", lw=2))

        ax_xy.set_xlabel("x"); ax_xy.set_ylabel("y"); ax_xy.set_title("x-y")
        ax_xz.set_xlabel("x"); ax_xz.set_ylabel("z"); ax_xz.set_title("x-z")
        for ax in (ax_xy, ax_xz):
            ax.set_aspect("equal", adjustable="datalim")
        ax_xy.legend(loc="upper right", fontsize=8, markerscale=1.5)

        n_inside, _ = count_inside_cube(snap, center, side)
        fig.suptitle(f"center=({center['x']:.6g}, {center['y']:.6g}, {center['z']:.6g}) "
                     f"[{center_mode}]  side={side:.6g}  --  "
                     f"{n_inside:,}/{n_total:,} particles inside")
        fig.tight_layout()
        try:
            plt.show(block=False)
            fig.canvas.draw()
            plt.pause(0.001)
        except Exception:  # noqa: BLE001 -- best-effort display, never fatal
            pass
        if not interactive_backend:
            frame_n += 1
            frame_path = f"{frame_save_prefix}_iter{frame_n}.png"
            fig.savefig(frame_path, dpi=130)
            print(f"  wrote {frame_path} -- open it to see the current box")

    print("\nInteractive sub-domain assistance -- x-y and x-z projections of every particle, "
          "with the geometric center (+), center of mass (x), densest region (*, if available), "
          "the last user-specified center (D, once you've entered one), and the current "
          "candidate box (red square, both panels) marked.")
    redraw()

    while True:
        print(f"\n  current: center=({center['x']:.6g}, {center['y']:.6g}, {center['z']:.6g}) "
              f"[{center_mode}]  side={side:.6g}")
        choice = prompt("  center: [g]eometric / [c]enter-of-mass / [d]ensest / [u]ser x,y,z / "
                         "Enter to keep current", "")
        choice = (choice or "").strip().lower()
        if choice.startswith("g"):
            center, center_mode = dict(geometric_center), "geometric"
        elif choice.startswith("c"):
            center, center_mode = dict(com_center), "com"
        elif choice.startswith("d"):
            if densest_center is None:
                eprint("  densest-region is not available, pick another option.")
                continue
            center, center_mode = dict(densest_center), "densest"
        elif choice.startswith("u"):
            coords = prompt("    x y z", "")
            try:
                cx, cy, cz = (float(v) for v in (coords or "").split())
            except ValueError:
                eprint("    could not parse three numbers, try again.")
                continue
            center, center_mode = {"x": cx, "y": cy, "z": cz}, "user"
            user_center = dict(center)
        elif choice:
            eprint("  not understood, try again.")
            continue

        prev_side = side
        side_in = prompt("  box side (cube side length)", f"{side:.6g}")
        try:
            new_side = float(side_in)
            if new_side <= 0:
                eprint("  box side must be positive, keeping the previous side.")
                side = prev_side
            else:
                side = new_side
        except (TypeError, ValueError):
            eprint("  could not parse a number, keeping the previous side.")
            side = prev_side

        redraw()

        done = prompt("  [a]ccept this box / Enter to keep adjusting / [q]uit without selecting", "")
        done = (done or "").strip().lower()
        if done.startswith("a"):
            break
        if done.startswith("q"):
            plt.close(fig)
            eprint("suggest_config.py: aborted during interactive sub-domain selection, "
                   "no config written.")
            sys.exit(1)

    final_path = f"{frame_save_prefix}_selection.png"
    fig.savefig(final_path, dpi=150)
    print(f"  wrote {final_path} (record of the selected sub-domain)")
    plt.close(fig)
    return center, center_mode, side


def _shift_plain_ascii(snapshot_path, center, output_path):
    """Write a copy of a plain_ascii raw snapshot with (x, y, z) re-centered
    on `center`; vx/vy/vz/m/h/rho columns pass through unchanged. Preserves
    the original file's column formatting precision (repr-round-trippable)."""
    with open(snapshot_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            parts = line.split()
            if len(parts) != 9:
                continue
            x, y, z = (float(parts[0]) - center["x"],
                       float(parts[1]) - center["y"],
                       float(parts[2]) - center["z"])
            rest = parts[3:]
            fout.write(f"{x!r} {y!r} {z!r} " + " ".join(rest) + "\n")


def _shift_splash_ascii(snapshot_path, center, output_path):
    """Same idea as _shift_plain_ascii, but for the splash_ascii format:
    the header/preamble (time/npart/units/column-header lines) is copied
    through byte-for-byte -- only x/y/z ever move, so units, particle
    counts, and the column layout all stay valid -- and every data row has
    just its x/y/z columns shifted, with every other column (m, h,
    density, velocities, pressure, alpha, div_v, poten, dt, itype, or any
    other physics-module column) passed through completely untouched,
    INCLUDING any non-gas (itype != 1) rows: those are passed through
    unshifted-in-identity-but-still-shifted-in-position rows rather than
    dropped here -- load_sph_snapshot()/load_snapshot() do the itype
    filtering, not this rewrite, so a shifted copy remains a complete,
    format-valid snapshot in its own right.

    Which token is x, which is y, which is z is read from the header line
    itself (same as _load_splash_ascii) rather than assumed to be the
    first three columns, since column order isn't guaranteed -- see
    _load_splash_ascii's docstring for why."""
    n_header = _splash_ascii_header_line_count(snapshot_path)
    with open(snapshot_path) as f:
        lines = f.readlines()
    header_line = lines[n_header - 1]
    n_data_cols = len(lines[n_header].split())
    labels = _split_splash_ascii_labels(header_line, n_data_cols, snapshot_path)
    col_index = {}
    for i, label in enumerate(labels):
        semantic = _SPLASH_ASCII_COLUMN_ALIASES.get(_normalize_splash_label(label))
        if semantic in ("x", "y", "z") and semantic not in col_index:
            col_index[semantic] = i
    for req in ("x", "y", "z"):
        if req not in col_index:
            raise ValueError(
                f"splash-ascii snapshot: could not find a '{req}' column in the "
                f"column-header line, needed to shift particle positions: {snapshot_path}")
    ix, iy, iz = col_index["x"], col_index["y"], col_index["z"]

    with open(snapshot_path) as fin, open(output_path, "w") as fout:
        for i, line in enumerate(fin):
            if i < n_header:
                fout.write(line)
                continue
            parts = line.split()
            if len(parts) != n_data_cols:
                continue  # defensive: skip any stray blank/short trailing line
            parts[ix] = repr(float(parts[ix]) - center["x"])
            parts[iy] = repr(float(parts[iy]) - center["y"])
            parts[iz] = repr(float(parts[iz]) - center["z"])
            fout.write(" ".join(parts) + "\n")


def shift_snapshot_file(snapshot_path, center, output_path, fmt):
    """Write a CoM-shifted copy of the raw snapshot (positions only) in
    whichever format `fmt` is -- "plain_ascii" or "splash_ascii" (never
    "auto"; resolve it before calling)."""
    if fmt == "plain_ascii":
        _shift_plain_ascii(snapshot_path, center, output_path)
    elif fmt == "splash_ascii":
        _shift_splash_ascii(snapshot_path, center, output_path)
    else:
        raise ValueError(f"unknown snapshot format: {fmt!r} (expected plain_ascii or "
                          f"splash_ascii)")


def prompt(question, default=None):
    suffix = f" [{default}]" if default is not None else ""
    answer = input(f"{question}{suffix}: ").strip()
    return answer if answer else default


def main():
    ap = argparse.ArgumentParser(
        description="Suggest a domain and write a kiwi config.json from a raw SPH snapshot.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("snapshot", help="path to the raw SPH snapshot (plain_ascii or "
                                      "splash_ascii format -- see --snapshot-format)")
    ap.add_argument("--snapshot-format", choices=["auto", "plain_ascii", "splash_ascii"],
                     default="auto",
                     help="raw SPH snapshot format (default: auto -- detected by sniffing the "
                          "file for the `splash to ascii` signature comment). Override if "
                          "detection ever guesses wrong.")
    ap.add_argument("--run-name", help="run_name for the config (default: snapshot basename)")
    ap.add_argument("--work-folder", help="work_folder for the config (no default -- always required)")
    ap.add_argument("--density-table", help="path to a tabulated kernel-integral file; "
                                             "passing this implies use_tabulated_integral=true")
    ap.add_argument("--mode", choices=["full", "zoom"], default="full",
                     help="'full' (default): the percentile-trimmed whole-domain suggestion. "
                          "'zoom': a smaller cubic sub-domain around a chosen center -- see "
                          "--center-mode/--center/--box-side and --zoom-interactive.")
    ap.add_argument("--trim-percentile", type=float, default=0.05,
                     help="[--mode full] per-tail percentile trimmed from each axis before "
                          "sizing the domain (default: 0.05, i.e. keep the [0.05, 99.95] "
                          "percentile range per axis). 0 disables trimming (domain covers every "
                          "particle). Also used as the interactive zoom assistance's starting "
                          "box-side guess when --box-side isn't given.")
    ap.add_argument("--pad-fraction", type=float, default=0.02,
                     help="[--mode full] extra fractional padding added to the cube side beyond "
                          "the trimmed extent, so no particle sits exactly on the domain wall "
                          "(default: 0.02)")
    ap.add_argument("--center-mode", choices=["geometric", "com", "densest", "user"],
                     help="[--mode zoom] how to choose the sub-domain center: 'geometric' "
                          "(bounding-box midpoint of all particles), 'com' (mass-weighted center "
                          "of mass), 'densest' (density-peak location, see "
                          "--densest-precision-fraction), or 'user' (exact coordinates via "
                          "--center). Given together with --box-side, fully specifies the "
                          "sub-domain and skips the interactive assistance (unless "
                          "--zoom-interactive is also passed). Leave both unset to launch the "
                          "interactive assistance instead.")
    ap.add_argument("--densest-precision-fraction", type=float, default=1.0 / 128.0,
                     help="[--mode zoom] target precision of the --center-mode densest search, "
                          "as a fraction of the whole particle set's bounding-cube side (default: "
                          "1/128). Smaller = finer final precision, more search rounds "
                          "(rounds ~= log2(1/this)); the search itself stays fast at any "
                          "reasonable value since each round only rescans the current box.")
    ap.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"),
                     help="[--mode zoom] exact sub-domain center coordinates; required with "
                          "--center-mode user, and used as the interactive assistance's starting "
                          "center (as 'user') if given without --center-mode")
    ap.add_argument("--box-side", type=float,
                     help="[--mode zoom] cube side length of the sub-domain; required (together "
                          "with --center-mode) to skip the interactive assistance")
    ap.add_argument("--zoom-interactive", action="store_true",
                     help="[--mode zoom] launch the interactive assistance even if --center-mode "
                          "and --box-side are both given (they're used as the starting point, so "
                          "you can preview/adjust before accepting)")
    ap.add_argument("--block-scale-factor", type=float, default=3.0,
                     help="grid_blocks.block_scale_factor (default: 3.0, the tool's own default)")
    ap.add_argument("--shift-snapshot", action="store_true",
                     help="NOT the default -- also write a CoM-shifted copy of the snapshot and "
                          "point the config at it instead (see module docstring)")
    ap.add_argument("--output", help="path to write the config JSON to "
                                      "(default: <run_name>.json in the current directory)")
    ap.add_argument("--non-interactive", action="store_true",
                     help="never prompt; fail instead if run_name/work_folder aren't supplied as "
                          "flags, or (in --mode zoom) if --center-mode and --box-side aren't both "
                          "supplied (the interactive assistance needs a terminal)")
    args = ap.parse_args()

    if not os.path.isfile(args.snapshot):
        eprint(f"suggest_config.py: no such file: {args.snapshot}")
        return 1

    try:
        resolved_format = args.snapshot_format
        if resolved_format == "auto":
            resolved_format = detect_snapshot_format(args.snapshot)
        snap = load_snapshot(args.snapshot, resolved_format)
    except Exception as e:  # noqa: BLE001 -- surface any parse error plainly, this is a CLI tool
        eprint(f"suggest_config.py: failed to read snapshot: {e}")
        return 1
    print(f"  snapshot format: {resolved_format}"
          + ("" if args.snapshot_format != "auto" else " (auto-detected)"))

    if args.mode == "zoom":
        geometric_center = compute_geometric_center(snap)
        com_center = compute_center_of_mass(snap)
        print(f"Inspecting {args.snapshot} ({snap['x'].size:,} particles) for sub-domain selection:")
        print(f"  geometric center: ({geometric_center['x']:.6g}, {geometric_center['y']:.6g}, "
              f"{geometric_center['z']:.6g})")
        print(f"  center of mass:   ({com_center['x']:.6g}, {com_center['y']:.6g}, "
              f"{com_center['z']:.6g})")
        densest_center = compute_densest_region(snap, args.densest_precision_fraction)
        print(f"  densest region:   ({densest_center['x']:.6g}, {densest_center['y']:.6g}, "
              f"{densest_center['z']:.6g})")
        print()

        manual_ready = args.center_mode is not None and args.box_side is not None

        if manual_ready and not args.zoom_interactive:
            try:
                center = resolve_manual_zoom_center(args.center_mode, args.center,
                                                     geometric_center, com_center, densest_center)
            except ValueError as e:
                eprint(f"suggest_config.py: {e}")
                return 1
            side = args.box_side
            center_mode = args.center_mode
        else:
            if args.non_interactive:
                if manual_ready:
                    eprint("suggest_config.py: --zoom-interactive requires a terminal, can't "
                           "combine with --non-interactive")
                else:
                    eprint("suggest_config.py: --mode zoom requires --center-mode and --box-side "
                           "when --non-interactive is set (the interactive assistance needs a "
                           "terminal and a display)")
                return 1

            # Starting point for the interactive assistance: whatever manual flags were
            # given, else the geometric center and the full-mode suggestion's own side as a
            # first guess.
            if args.center_mode is not None:
                try:
                    initial_center = resolve_manual_zoom_center(
                        args.center_mode, args.center, geometric_center, com_center,
                        densest_center)
                    initial_center_mode = args.center_mode
                except ValueError as e:
                    eprint(f"suggest_config.py: {e}")
                    return 1
            elif args.center is not None:
                initial_center = {"x": args.center[0], "y": args.center[1], "z": args.center[2]}
                initial_center_mode = "user"
            else:
                initial_center = geometric_center
                initial_center_mode = "geometric"

            if args.box_side is not None:
                initial_side = args.box_side
            else:
                _, fallback_info = suggest_domain(snap, args.trim_percentile, args.pad_fraction)
                initial_side = fallback_info["side"]

            frame_save_prefix = os.path.splitext(os.path.basename(args.snapshot))[0] + "_zoom"
            center, center_mode, side = run_interactive_zoom_assistance(
                snap, geometric_center, com_center, densest_center,
                initial_center, initial_center_mode, initial_side, frame_save_prefix,
            )

        domain = domain_from_center_side(center, side)
        n_inside, n_total = count_inside_cube(snap, center, side)
        print_zoom_summary(center, center_mode, side, n_inside, n_total)
        print()
    else:
        domain, info = suggest_domain(snap, args.trim_percentile, args.pad_fraction)
        print(f"Inspecting {args.snapshot}:")
        print_domain_summary(info)
        print()

    domain_center = {
        "x": (domain["x_min"] + domain["x_max"]) / 2.0,
        "y": (domain["y_min"] + domain["y_max"]) / 2.0,
        "z": (domain["z_min"] + domain["z_max"]) / 2.0,
    }

    default_run_name = os.path.splitext(os.path.basename(args.snapshot))[0]
    run_name = args.run_name
    if run_name is None:
        if args.non_interactive:
            eprint("suggest_config.py: --run-name is required with --non-interactive")
            return 1
        run_name = prompt("run_name", default_run_name)

    work_folder = args.work_folder
    if work_folder is None:
        if args.non_interactive:
            eprint("suggest_config.py: --work-folder is required with --non-interactive")
            return 1
        work_folder = prompt("work_folder (base output directory, no default)")
        while not work_folder:
            eprint("  work_folder is required.")
            work_folder = prompt("work_folder (base output directory, no default)")

    density_table = args.density_table
    use_tabulated = density_table is not None
    if density_table is None and not args.non_interactive:
        want_tabulated = prompt(
            "use the tabulated kernel integral instead of the analytic default? (y/N)", "n"
        )
        if want_tabulated.strip().lower().startswith("y"):
            density_table = prompt("density_table path")
            use_tabulated = bool(density_table)

    input_snapshot = args.snapshot
    if args.shift_snapshot:
        shifted_path = os.path.join(
            os.path.dirname(os.path.abspath(args.snapshot)) or ".",
            os.path.splitext(os.path.basename(args.snapshot))[0] + "_comshift.out",
        )
        print(f"\n--shift-snapshot: writing CoM-shifted copy to {shifted_path} ...")
        shift_snapshot_file(args.snapshot, domain_center, shifted_path, resolved_format)
        print("  done.")
        input_snapshot = shifted_path
        half_x = (domain["x_max"] - domain["x_min"]) / 2.0
        half_y = (domain["y_max"] - domain["y_min"]) / 2.0
        half_z = (domain["z_max"] - domain["z_min"]) / 2.0
        domain = {
            "x_min": -half_x, "x_max": half_x,
            "y_min": -half_y, "y_max": half_y,
            "z_min": -half_z, "z_max": half_z,
        }

    config = {
        "run_name": run_name,
        "domain": domain,
        "grid_blocks": {"block_scale_factor": args.block_scale_factor},
        "input_snapshot": input_snapshot,
        "input_snapshot_format": resolved_format,
        "work_folder": work_folder,
        "buffer": {"enabled": True, "strategy": "coverage_grid"},
        "lloyd": {"iterations": 0},
        "cartesian_map": {"resolution": 128},
    }
    if use_tabulated:
        config["mapping"] = {"use_tabulated_integral": True, "density_table": density_table}

    output_path = args.output or f"{run_name}.json"
    with open(output_path, "w") as f:
        json.dump(config, f, indent=2)
        f.write("\n")

    print(f"\nWrote {output_path}:")
    print(json.dumps(config, indent=2))
    print(f"\nDefaults used (not prompted): buffer.strategy=coverage_grid, "
          f"lloyd.iterations=0, mapping.use_tabulated_integral={str(use_tabulated).lower()}, "
          f"cartesian_map.resolution=128. "
          f"Edit {output_path} directly to change any of these.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
