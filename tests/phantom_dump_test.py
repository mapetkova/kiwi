#!/usr/bin/env python3
"""Synthetic Phantom-binary-dump test for KIWI's direct dump reader (standard library only).

Writes small Phantom dumps in several layouts (byte order, integer/real widths, extra particle types,
dead particles, sink blocks, MPI-style repeated blocks, units), writes the same particles as a plain
`plain_ascii` text file, runs `kiwi sph-to-grid` on both and requires the two grids to agree. Also
checks that malformed dumps are rejected with an understandable message.

The writer follows the format description at
https://phantomsph.readthedocs.io/en/latest/user-guide/dumpfile.html (Fortran unformatted records,
capture pattern, tagged header of eight typed sections, per-block tagged arrays).

Usage: phantom_dump_test.py <path-to-kiwi-executable> [work_dir]
"""
import json
import math
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile

HALF = 6.0
RES = 24
SPLASH_UNITS = {"length_cm": 3.086e18, "mass_g": 1.989e33, "velocity_cms": 1.0e5}  # KIWI's defaults


# ------------------------------------------------------------------ dump writer

def f32(x):
    return struct.unpack("<f", struct.pack("<f", x))[0]


class DumpWriter:
    """Writes a Phantom dump. `e` is '<' or '>', `ib`/`rb` the default int/real widths in bytes."""

    def __init__(self, e="<", ib=4, rb=8):
        self.e, self.ib, self.rb = e, ib, rb
        self.fmt = [("i" if ib == 4 else "q"), "b", "h", "i", "q", ("d" if rb == 8 else "f"), "f", "d"]

    def rec(self, payload):
        n = struct.pack(self.e + "i", len(payload))
        return n + payload + n

    def capture(self):
        i, r = self.fmt[0], self.fmt[5]
        return self.rec(struct.pack(self.e + i + r + i + i + i, 60769, 60878.0, 60878, 1, 690706))

    def header(self, sections):
        """sections: {kind: [(tag, value), ...]} for kinds 0..7."""
        out = b""
        for k in range(8):
            items = sections.get(k, [])
            out += self.rec(struct.pack(self.e + "i", len(items)))
            if items:
                out += self.rec(b"".join(t.ljust(16).encode() for t, _ in items))
                out += self.rec(struct.pack(self.e + str(len(items)) + self.fmt[k], *[v for _, v in items]))
        return out

    def blocks(self, blocks):
        """blocks: [(n, {kind: [(name, values), ...]}), ...]."""
        out = self.rec(struct.pack(self.e + "i", len(blocks)))
        for n, arrays in blocks:
            nums = [len(arrays.get(k, [])) for k in range(8)]
            out += self.rec(struct.pack(self.e + "q8i", n, *nums))
        for n, arrays in blocks:
            for k in range(8):
                for name, vals in arrays.get(k, []):
                    out += self.rec(name.ljust(16).encode())
                    out += self.rec(struct.pack(self.e + str(len(vals)) + self.fmt[k], *vals))
        return out


def make_particles(rng, n):
    rows = []
    while len(rows) < n:
        p = [rng.gauss(0, 1.0) for _ in range(3)]
        if sum(c * c for c in p) > 6.0:
            continue
        v = [rng.uniform(-3, 3) for _ in range(3)]
        h = rng.uniform(0.25, 0.6)
        rows.append(p + v + [h])
    return rows


def write_case(path, rows, *, e="<", ib=4, rb=8, units=(1.0, 1.0, 1.0), mass=0.02, hfact=1.2, dust_every=0,
               dead=3, mpi=1, ident="FT:Phantom:test (hydro): 01/01/2026", drop=(), real4_positions=False,
               with_itype=None):
    """Writes the dump `path` and returns the particle rows a correct reader must keep (code units already
    rounded to the precision stored in the file). `dust_every`>0 makes every n-th particle type 2."""
    w = DumpWriter(e, ib, rb)
    posf = f32 if (rb == 4 or real4_positions) else (lambda x: x)
    realk = 6 if real4_positions else 5            # kind holding x,y,z,vx,vy,vz
    types = [2 if (dust_every and (i + 1) % dust_every == 0) else 1 for i in range(len(rows))]
    use_itype = (dust_every > 0) if with_itype is None else with_itype
    keep = []
    allp = []                                      # (row, type, alive) in file order
    for r, t in zip(rows, types):
        allp.append((r, t, True))
    for k in range(dead):                          # dead particles: h < 0, positions inside the box
        allp.append(([0.1 * k, 0.2, 0.3, 1.0, 1.0, 1.0, -1.0], 1, False))
    rng = random.Random(7)
    rng.shuffle(allp)                              # dead particles are interleaved with live ones
    nper = math.ceil(len(allp) / mpi)
    hdr_ints = [("ntypes", 8), ("nblocks", mpi)]
    npt = [sum(1 for r, t, a in allp if t == 1)] + [sum(1 for r, t, a in allp if t == 2)] + [0] * 6
    stored = f32 if rb == 4 else (lambda x: x)      # massoftype and hfact live in the default real
    mass, hfact = stored(mass), stored(hfact)
    massoftype = [mass, stored(mass * 3.0)] + [0.0] * 6
    sections = {0: hdr_ints, 4: [("npartoftype", v) for v in npt],
                5: [("massoftype", v) for v in massoftype] + [("hfact", hfact)],
                7: [("udist", units[0]), ("umass", units[1]), ("utime", units[2])]}
    blocks = []
    for m in range(mpi):
        chunk = allp[m * nper:(m + 1) * nper]
        cols = {c: [] for c in ("x", "y", "z", "vx", "vy", "vz", "h", "itype")}
        for r, t, alive in chunk:
            for c, v in zip(("x", "y", "z", "vx", "vy", "vz"), r[:6]):
                cols[c].append(posf(v))
            cols["h"].append(f32(r[6]))
            cols["itype"].append(t)
            if alive and t == 1:
                keep.append([posf(v) for v in r[:6]] + [f32(r[6])])
        arrays = {realk: [(c, cols[c]) for c in ("x", "y", "z", "vx", "vy", "vz") if c not in drop]}
        arrays[6] = arrays.get(6, []) + [("h", cols["h"])]
        if use_itype:
            arrays[1] = [("itype", cols["itype"])]
        blocks.append((len(chunk), arrays))
        sink = {5: [("x", [50.0, 60.0]), ("y", [0.0, 0.0]), ("z", [0.0, 0.0]), ("m", [9.0, 9.0]), ("h", [1.0, 1.0])]}
        if rb == 4:
            sink = {5: [(n_, [f32(v) for v in vals]) for n_, vals in sink[5]]}
        blocks.append((2, sink))
    data = w.capture() + w.rec(ident.ljust(100).encode()) + w.header(sections) + w.blocks(blocks)
    with open(path, "wb") as f:
        f.write(data)
    return keep, massoftype[0], hfact


# ------------------------------------------------------------------ running KIWI

def run(kiwi, stage, cfg_path, expect_fail=None):
    r = subprocess.run([kiwi, stage, cfg_path], capture_output=True, text=True)
    if expect_fail is not None:
        ok = r.returncode != 0 and expect_fail in r.stderr
        return ok, r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ""
    if r.returncode != 0:
        print(f"FAIL: kiwi {stage} exited with {r.returncode}\n{r.stderr}")
        sys.exit(1)
    return r


def load_grid(folder):
    p = os.path.join(folder, os.listdir(folder)[0], "cartesian_grid.dat")
    with open(p, "rb") as f:
        assert f.read(8) == b"VGTCGRD2"
        dims = struct.unpack("<3i", f.read(12))
        f.read(48)
        n = dims[0] * dims[1] * dims[2]
        return [struct.unpack(f"<{n}d", f.read(8 * n)) for _ in range(4)]


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_phdump_")
    os.makedirs(work, exist_ok=True)
    rows = make_particles(random.Random(3), 150)
    domain = {f"{a}_{b}": s * HALF for a in "xyz" for b, s in (("min", -1), ("max", 1))}

    def cfg(tag, snapshot, fmt, extra=None):
        c = {"run_name": tag, "domain": domain, "input_snapshot": snapshot, "input_snapshot_format": fmt,
             "work_folder": os.path.join(work, "runs_" + tag), "cartesian_grid": {"resolution": RES}}
        c.update(extra or {})
        shutil.rmtree(c["work_folder"], ignore_errors=True)
        p = os.path.join(work, tag + ".json")
        with open(p, "w") as f:
            json.dump(c, f)
        return c, p

    def compare(tag, dump_kwargs, kiwi_units, ascii_conv, kiwi_extra=None):
        """Dump vs equivalent ascii. ascii_conv converts code-unit values to the units KIWI works in."""
        dump = os.path.join(work, tag + ".dump")
        keep, m_code, hfact = write_case(dump, rows, **dump_kwargs)
        asc = os.path.join(work, tag + ".txt")
        with open(asc, "w") as f:
            for x, y, z, vx, vy, vz, h in keep:
                X, Y, Z, VX, VY, VZ, M, H = ascii_conv(x, y, z, vx, vy, vz, m_code, h)
                f.write(f"{X!r} {Y!r} {Z!r} {VX!r} {VY!r} {VZ!r} {M!r} {H!r} 1.0\n")
        extra = dict(kiwi_extra or {})
        if kiwi_units:
            extra["phantom_units"] = kiwi_units
        c1, p1 = cfg(tag + "_dump", dump, "auto", extra)
        c2, p2 = cfg(tag + "_txt", asc, "plain_ascii")
        run(kiwi, "sph-to-grid", p1)
        run(kiwi, "sph-to-grid", p2)
        a, b = load_grid(c2["work_folder"]), load_grid(c1["work_folder"])
        mx = max(a[0])
        ok = check(f"{tag}: kept {len(keep)} particles, density grid == ascii grid",
                   max(abs(p - q) for p, q in zip(a[0], b[0])) <= 1e-9 * mx,
                   f"max diff/max = {max(abs(p - q) for p, q in zip(a[0], b[0])) / mx:.1e}")
        sig = [i for i, r_ in enumerate(a[0]) if r_ > 1e-6 * mx]
        dv = max(max(abs(a[k][i] - b[k][i]) for i in sig) for k in (1, 2, 3))
        ok &= check(f"{tag}: velocity grid == ascii grid", dv <= 1e-9, f"max |dv| = {dv:.1e}")
        return ok

    ident = lambda x, y, z, vx, vy, vz, m, h: (x, y, z, vx, vy, vz, m, h)
    ok = True
    print("dumps vs the same particles as plain ascii")
    ok &= compare("le_int32_real64", {}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("big_endian", {"e": ">"}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("int64_real32", {"ib": 8, "rb": 4}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("big_endian_int64_real64", {"e": ">", "ib": 8, "rb": 8}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("real4_positions", {"real4_positions": True}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("dust_type_dropped", {"dust_every": 4}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("mpi_blocks", {"mpi": 3}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    ok &= compare("itype_all_gas", {"with_itype": True}, {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, ident)
    # a small dump (real*4 everywhere, no velocity arrays) is usable for density when explicitly allowed
    zero_v = lambda x, y, z, vx, vy, vz, m, h: (x, y, z, 0.0, 0.0, 0.0, m, h)
    ok &= compare("small_dump_zero_velocity", {"rb": 4, "drop": ("vx", "vy", "vz"), "ident": "ST:Phantom:test (hydro)"},
                  {"length_cm": 1.0, "mass_g": 1.0, "velocity_cms": 1.0}, zero_v, {"allow_missing_velocity": True})
    # unit conversion with KIWI's default (SPLASH) constants: code units -> pc, M_Sun, km/s
    ud, um, ut = 3.086e18, 1.9891e38, 1.488114062362035e12
    def conv(x, y, z, vx, vy, vz, m, h):
        lf, vf = ud / SPLASH_UNITS["length_cm"], (ud / ut) / SPLASH_UNITS["velocity_cms"]
        return x * lf, y * lf, z * lf, vx * vf, vy * vf, vz * vf, m * um / SPLASH_UNITS["mass_g"], h * lf
    ok &= compare("default_units", {"units": (ud, um, ut)}, None, conv)

    print("malformed dumps are rejected clearly")
    bad = {
        "untagged": ({"ident": "FU:Phantom:old"}, "tagged"),
        "small_dump_without_velocities": ({"drop": ("vx", "vy", "vz")}, "allow_missing_velocity"),
        "missing_position_array": ({"drop": ("z",)}, "'z' not found"),
    }
    for tag, (kw, msg) in bad.items():
        dump = os.path.join(work, tag + ".dump")
        write_case(dump, rows, **kw)
        c, p = cfg(tag, dump, "auto")
        good, last = run(kiwi, "sph-to-grid", p, expect_fail=msg)
        ok &= check(f"{tag}: rejected with '{msg}'", good, last[:110])
    # several particle types listed in the header but no itype array
    dump = os.path.join(work, "no_itype.dump")
    write_case(dump, rows, dust_every=4, with_itype=False)
    c, p = cfg("no_itype", dump, "auto")
    good, last = run(kiwi, "sph-to-grid", p, expect_fail="itype")
    ok &= check("several types without an itype array: rejected", good, last[:110])
    junk = os.path.join(work, "junk.dump")
    with open(junk, "wb") as f:
        f.write(os.urandom(4096))
    c, p = cfg("junk", junk, "phantom_binary")
    good, last = run(kiwi, "sph-to-grid", p, expect_fail="capture")
    ok &= check("not a Phantom dump: rejected", good, last[:110])

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
