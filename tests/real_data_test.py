#!/usr/bin/env python3
"""Real-data test for KIWI: a Phantom/SPLASH snapshot from an actual simulation (standard library only).

The snapshot lives in tests/data/ (not tracked by git: it is hundreds of MB). It was run with the cubic
spline kernel; the non-cubic kernels are run on the same data purely as a stress test -- the physics is
"wrong" for them but mass and momentum conservation do not depend on that, so they still must hold.

If the Phantom binary dump the ascii was exported from sits next to it (same name without ".ascii"),
the cubic-spline run is repeated from the binary and must give the same grid as the ascii run.

If a Phantom small dump of the preceding step sits next to it (brick_03926: positions and h only, no
velocities), it is run with allow_missing_velocity and must give all-zero velocities and a density grid
close to the full dump's (the two steps are ~130 yr apart).

Checks, all on the full snapshot with sph-to-grid:
  * cubic spline (default fast path): total mass and total momentum on the grid equal the particle sums
    (the domain is padded so that no kernel is cut off), density >= 0, nothing NaN;
  * cubic spline forced through the generic integral (kernel.generic_path) == the default path;
  * tabulated cubic integral: close to the analytic one (mass within a loose tolerance);
  * quintic spline and Wendland C6 (the widest and the highest-order kernel): mass and momentum conserved.

Usage: real_data_test.py <path-to-kiwi-executable> [snapshot] [work_dir]
Exit code 0 on success, 1 on failure, 77 (ctest "skipped") if the snapshot is not there.
"""
import json
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_SNAPSHOT = os.path.join(HERE, "data", "brick_03927.ascii")
CELL = 1.0  # grid cell size in the snapshot's length unit (pc)


# Column order of the `splash to ascii` header in this snapshot (checked against the header line below).
COLUMNS = ["x", "y", "z", "particle mass", "h", "density", "v_x", "v_y", "v_z"]
ITYPE = 14  # last column


def read_snapshot(path):
    """Gas rows (itype == 1) of a `splash to ascii` file: particle count, mass, momentum, sum m|v|,
    extents and max h (all in the file's own units, which is all the conservation checks need)."""
    n = 0
    msum = 0.0
    mom = [0.0, 0.0, 0.0]
    mabsv = 0.0
    lo, hi = [math.inf] * 3, [-math.inf] * 3
    hmax = 0.0
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                if line.startswith("# x ["):  # the column-header line
                    pos = [line.index(name) for name in COLUMNS]
                    assert pos == sorted(pos) and line.rstrip().endswith("itype"), "unexpected column layout"
                continue
            t = line.split()
            if int(float(t[ITYPE])) != 1:
                continue
            m, h = float(t[3]), float(t[4])
            v = [float(t[6]), float(t[7]), float(t[8])]
            n += 1
            msum += m
            for d in range(3):
                x = float(t[d])
                mom[d] += m * v[d]
                lo[d] = min(lo[d], x)
                hi[d] = max(hi[d], x)
            mabsv += m * math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)
            hmax = max(hmax, h)
    return dict(n=n, mass=msum, momentum=mom, mabsv=mabsv, lo=lo, hi=hi, hmax=hmax)


def load_grid(path):
    with open(path, "rb") as f:
        assert f.read(8) == b"VGTCGRD2", "bad magic"
        nx, ny, nz = struct.unpack("<3i", f.read(12))
        b = struct.unpack("<6d", f.read(48))
        n = nx * ny * nz
        arrays = [struct.unpack(f"<{n}d", f.read(8 * n)) for _ in range(4)]
    return (nx, ny, nz), b, arrays


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    snap_path = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_SNAPSHOT
    if not os.path.exists(snap_path):
        print(f"SKIP: snapshot not found: {snap_path}")
        return 77
    work = sys.argv[3] if len(sys.argv) > 3 else tempfile.mkdtemp(prefix="kiwi_real_")
    os.makedirs(work, exist_ok=True)

    t0 = time.time()
    s = read_snapshot(snap_path)
    print(f"snapshot: {s['n']} gas particles, mass {s['mass']:.6g}, max h {s['hmax']:.3g}  "
          f"({time.time() - t0:.1f} s to read)")
    pad = 3.0 * s["hmax"] * 1.02  # widest kernel supported (quintic: 3 h)
    dom = {}
    res = {}
    for d, ax in enumerate("xyz"):
        a, b = s["lo"][d] - pad, s["hi"][d] + pad
        n = max(1, int(math.ceil((b - a) / CELL)))
        b = a + n * CELL  # whole number of cells
        dom[ax + "_min"], dom[ax + "_max"], res["resolution_" + ax] = a, b, n
    print(f"domain padded by {pad:.3g}: {res}")

    def run(tag, kernel=None, tabulated=False, snapshot=None, extra=None):
        cfg = {"run_name": tag, "domain": dom, "input_snapshot": snapshot or snap_path,
               "work_folder": os.path.join(work, tag),
               "cartesian_grid": dict(res, use_tabulated_integral=tabulated)}
        if kernel:
            cfg["kernel"] = kernel
        cfg.update(extra or {})
        shutil.rmtree(cfg["work_folder"], ignore_errors=True)
        path = os.path.join(work, tag + ".json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=1)
        t = time.time()
        r = subprocess.run([kiwi, "sph-to-grid", path], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: kiwi sph-to-grid ({tag}) exited with {r.returncode}\n{r.stderr[-2000:]}")
            sys.exit(1)
        out = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "cartesian_grid.dat")
        print(f"  ran {tag} in {time.time() - t:.1f} s")
        return load_grid(out)

    def conservation(label, grid, mass_tol, mom_tol):
        dims, b, (rho, vx, vy, vz) = grid
        vol = ((b[1] - b[0]) / dims[0]) * ((b[3] - b[2]) / dims[1]) * ((b[5] - b[4]) / dims[2])
        mass = sum(rho) * vol
        mom = [sum(r * v for r, v in zip(rho, vv)) * vol for vv in (vx, vy, vz)]
        ok = check(f"{label}: no NaN/inf, density >= 0",
                   all(math.isfinite(x) for x in rho + vx + vy + vz) and min(rho) >= 0.0)
        rel = abs(mass - s["mass"]) / s["mass"]
        ok &= check(f"{label}: total mass", rel < mass_tol, f"rel.err={rel:.2e} tol={mass_tol:.0e}")
        dmom = max(abs(a - c) for a, c in zip(mom, s["momentum"])) / s["mabsv"]
        ok &= check(f"{label}: total momentum", dmom < mom_tol, f"err/sum(m|v|)={dmom:.2e} tol={mom_tol:.0e}")
        return ok

    ok = True
    print("cubic spline, default (analytic) path")
    cubic = run("cubic")
    ok &= conservation("cubic", cubic, 1e-9, 1e-9)

    binary = snap_path[:-len(".ascii")] if snap_path.endswith(".ascii") else None
    if binary and os.path.exists(binary):
        print("cubic spline from the Phantom binary dump the ascii was exported from")
        bgrid = run("cubic_binary", snapshot=binary)
        a, b = cubic[2], bgrid[2]
        sig = [i for i, r in enumerate(a[0]) if r > 1e-6 * max(a[0])]
        dr = max(abs(a[0][i] - b[0][i]) for i in sig) / max(a[0])
        dv = max(max(abs(a[k][i] - b[k][i]) for i in sig) for k in (1, 2, 3))
        ok &= check("binary == ascii: density", dr < 1e-9, f"max diff/max = {dr:.1e}")
        ok &= check("binary == ascii: velocity (cells with rho > 1e-6 max)", dv < 1e-6, f"max |dv| = {dv:.1e} km/s")

    small = binary[:-1] + "6" if binary else None
    if small and os.path.exists(small) and os.path.exists(binary):
        print("Phantom small dump of the previous step (no velocities), density only")
        sgrid = run("cubic_small", snapshot=small, extra={"allow_missing_velocity": True})
        a, s_ = cubic[2], sgrid[2]
        vol = ((cubic[1][1] - cubic[1][0]) / cubic[0][0]) * ((cubic[1][3] - cubic[1][2]) / cubic[0][1]) * \
              ((cubic[1][5] - cubic[1][4]) / cubic[0][2])
        ok &= check("small dump: velocities are all zero", all(v == 0.0 for k in (1, 2, 3) for v in s_[k]))
        mass_ratio = sum(s_[0]) / sum(a[0])
        ok &= check("small dump: total mass close to the full dump's", abs(mass_ratio - 1.0) < 2e-3,
                    f"ratio {mass_ratio:.6f}")
        l1 = sum(abs(p - q) for p, q in zip(a[0], s_[0])) / sum(a[0])
        ok &= check("small dump: density grid close to the full dump's", l1 < 0.10, f"L1 difference {l1:.2%}")

    print("cubic spline through the generic integral")
    gen = run("cubic_generic", {"generic_path": True})
    rho_c, rho_g = cubic[2][0], gen[2][0]
    worst = max(abs(a - b) for a, b in zip(rho_c, rho_g)) / max(rho_c)
    ok &= check("generic == fast path", worst < 1e-8, f"max |d density|/max density = {worst:.2e}")

    print("cubic spline, tabulated integral")
    tab = run("cubic_tab", tabulated=True)
    dims, b, (rho_t, *_rest) = tab
    vol = ((b[1] - b[0]) / dims[0]) * ((b[3] - b[2]) / dims[1]) * ((b[5] - b[4]) / dims[2])
    rel = abs(sum(rho_t) * vol - s["mass"]) / s["mass"]
    ok &= check("tabulated: total mass (approximate by design)", rel < 1e-3, f"rel.err={rel:.2e} tol=1e-3")

    for name in ("quintic_spline", "wendland_c6"):
        print(f"{name} on the same data (conservation stress test)")
        ok &= conservation(name, run(name, {"name": name}), 1e-9, 1e-9)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
