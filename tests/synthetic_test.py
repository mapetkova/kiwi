#!/usr/bin/env python3
"""Synthetic end-to-end test for KIWI (standard library only, no numpy needed).

Builds a small, fully known SPH "snapshot" (a cloud of equal-mass particles moving with a
uniform bulk velocity), runs KIWI on it, and checks physical invariants that must hold:

  sph-to-grid   total mass on the grid == total particle mass (kernel fully inside the domain),
                every cell with density has the bulk velocity, no NaN / negative density.
  cartesian-map the same checks through the full Voronoi pipeline
                (point-set -> tessellate -> map -> cartesian-map), with a looser mass tolerance
                because that path samples the Voronoi cells onto the Cartesian grid.
  kernels       the same two checks for every other kernel (quartic, quintic, Wendland C2/C4/C6),
                with each snapshot's h scaled so the kernel support radius is the same as the cubic
                spline's; plus the cubic spline forced through the generic integral
                (kernel.generic_path) must reproduce the default cubic result.

Usage: synthetic_test.py <path-to-kiwi-executable> [work_dir]
Exit code 0 on success, 1 on any failure (this is what ctest looks at).
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

N_PART = 400
PART_MASS = 0.01
BULK_V = (2.0, -1.0, 0.5)
HALF_DOMAIN = 6.0  # domain is [-6, 6]^3, a cube
RES = 48


# name -> support radius in units of h (SPLASH/Phantom convention, the config default)
KERNELS = {"quartic_spline": 2.5, "quintic_spline": 3.0, "wendland_c2": 2.0, "wendland_c4": 2.0,
           "wendland_c6": 2.0}
# the (slower) Voronoi pipeline is run for the kernels that stress it most: the widest support and the
# highest polynomial order
VORONOI_KERNELS = ("quintic_spline", "wendland_c6")


def write_snapshot(path, h_scale=1.0):
    """dale_phantom format: `x y z vx vy vz m h rho`, no header. Particles sit in a Gaussian blob
    (sigma 1, truncated at 2.5 sigma) well inside the domain; h is set from a local-density
    estimate that is only roughly right, which is fine -- mass conservation doesn't need an
    accurate h, just a positive one."""
    rng = random.Random(12345)
    rows = []
    while len(rows) < N_PART:
        p = [rng.gauss(0.0, 1.0) for _ in range(3)]
        if math.sqrt(sum(c * c for c in p)) > 2.5:
            continue
        rows.append(p)
    with open(path, "w") as f:
        for x, y, z in rows:
            r2 = x * x + y * y + z * z
            rho = N_PART * PART_MASS / (2 * math.pi) ** 1.5 * math.exp(-0.5 * r2)
            rho = max(rho, 1e-3)
            h = h_scale * 1.2 * (PART_MASS / rho) ** (1.0 / 3.0)
            f.write(f"{x:.8e} {y:.8e} {z:.8e} {BULK_V[0]} {BULK_V[1]} {BULK_V[2]} "
                    f"{PART_MASS:.8e} {h:.8e} {rho:.8e}\n")
    return N_PART * PART_MASS


def run(kiwi, stage, cfg_path):
    r = subprocess.run([kiwi, stage, cfg_path], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: `kiwi {stage}` exited with {r.returncode}\n{r.stderr}")
        sys.exit(1)


def read_doubles(f, n):
    data = f.read(8 * n)
    if len(data) != 8 * n:
        raise ValueError("truncated file")
    return struct.unpack(f"<{n}d", data)


def load_grid(path, magic, res_fields):
    with open(path, "rb") as f:
        if f.read(8) != magic:
            raise ValueError(f"{path}: bad magic")
        dims = struct.unpack(f"<{res_fields}i", f.read(4 * res_fields))
        bounds = read_doubles(f, 6)
        n = 1
        for d in dims:
            n *= d
        density, vx, vy, vz = (read_doubles(f, n) for _ in range(4))
    return dims, bounds, density, vx, vy, vz


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def verify(label, grid, true_mass, mass_tol):
    dims, b, density, vx, vy, vz = grid
    vol = ((b[1] - b[0]) / dims[0]) * ((b[3] - b[2]) / dims[1]) * ((b[5] - b[4]) / dims[2])
    total = sum(density) * vol
    ok = True
    ok &= check(f"{label}: no NaN/inf", all(math.isfinite(x) for x in density + vx + vy + vz))
    ok &= check(f"{label}: density >= 0", min(density) >= 0.0, f"min={min(density):.3g}")
    rel = abs(total - true_mass) / true_mass
    ok &= check(f"{label}: total mass", rel < mass_tol,
                f"grid={total:.6g} true={true_mass:.6g} rel.err={rel:.2e} tol={mass_tol:.0e}")
    thresh = 1e-6 * max(density)
    worst = 0.0
    for d, a, b_, c in zip(density, vx, vy, vz):
        if d > thresh:
            worst = max(worst, abs(a - BULK_V[0]), abs(b_ - BULK_V[1]), abs(c - BULK_V[2]))
    ok &= check(f"{label}: bulk velocity recovered", worst < 1e-6, f"max dev={worst:.2e}")
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_synth_")
    os.makedirs(work, exist_ok=True)
    snap = os.path.join(work, "synthetic_snapshot.out")
    true_mass = write_snapshot(snap)

    def config(name, extra):
        cfg = {
            "run_name": name,
            "domain": {"x_min": -HALF_DOMAIN, "x_max": HALF_DOMAIN, "y_min": -HALF_DOMAIN,
                       "y_max": HALF_DOMAIN, "z_min": -HALF_DOMAIN, "z_max": HALF_DOMAIN},
            "input_snapshot": snap,
            "input_snapshot_format": "dale_phantom",
            "work_folder": os.path.join(work, "runs_" + name),
            "buffer": {"enabled": False},
            "lloyd": {"iterations": 0, "save_each_iteration": False},
        }
        cfg.update(extra)
        # KIWI caches every stage's output, so start from a clean run folder or a stale result
        # could make this test pass without actually running anything.
        shutil.rmtree(cfg["work_folder"], ignore_errors=True)
        path = os.path.join(work, name + ".json")
        with open(path, "w") as f:
            json.dump(cfg, f, indent=2)
        return cfg, path

    ok = True

    print("sph-to-grid (direct particles -> Cartesian grid)")
    cfg, path = config("direct", {"cartesian_grid": {"resolution": RES}})
    run(kiwi, "sph-to-grid", path)
    out = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "cartesian_grid.dat")
    out_direct = out
    dims, *rest = load_grid(out, b"VGTCGRD2", 3)  # magic kept for file-format compatibility
    ok &= check("direct: grid shape", dims == (RES, RES, RES), str(dims))
    ok &= verify("direct", (dims, *rest), true_mass, 1e-6)

    print("Voronoi pipeline (point-set -> tessellate -> map -> cartesian-map)")
    cfg, path = config("voronoi", {"cartesian_map": {"resolution": RES}})
    run(kiwi, "cartesian-map", path)
    out = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "cartesian_map.dat")
    with open(out, "rb") as f:
        assert f.read(8) == b"VGTCMAP1"
        res = struct.unpack("<i", f.read(4))[0]
        b = read_doubles(f, 6)
        n = res ** 3
        arrays = [read_doubles(f, n) for _ in range(4)]
    ok &= check("voronoi: grid shape", res == RES, str(res))
    ok &= verify("voronoi", ((res, res, res), b, *arrays), true_mass, 2e-2)

    print("cubic spline through the generic integral (kernel.generic_path) vs the default fast path")
    cfg_g, path_g = config("direct_generic", {"cartesian_grid": {"resolution": RES}, "kernel": {"generic_path": True}})
    run(kiwi, "sph-to-grid", path_g)
    out_g = os.path.join(cfg_g["work_folder"], os.listdir(cfg_g["work_folder"])[0], "cartesian_grid.dat")
    ref, gen = load_grid(out_direct, b"VGTCGRD2", 3), load_grid(out_g, b"VGTCGRD2", 3)
    dmax = max(ref[2])
    worst = max(abs(a - b) for a, b in zip(ref[2], gen[2])) / dmax
    ok &= check("generic cubic == fast-path cubic", worst < 1e-8, f"max |d density|/max density = {worst:.2e}")

    for name, radkern in KERNELS.items():
        print(f"kernel {name} (support {radkern} h)")
        ksnap = os.path.join(work, f"synthetic_snapshot_{name}.out")
        kmass = write_snapshot(ksnap, h_scale=2.0 / radkern)

        def kconfig(tag, extra):
            c, p = config(tag, dict(extra, kernel={"name": name}))
            c["input_snapshot"] = ksnap
            with open(p, "w") as f:
                json.dump(c, f, indent=2)
            return c, p

        cfg, path = kconfig(name + "_direct", {"cartesian_grid": {"resolution": RES}})
        run(kiwi, "sph-to-grid", path)
        gout = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "cartesian_grid.dat")
        ok &= verify(name + " direct", load_grid(gout, b"VGTCGRD2", 3), kmass, 1e-6)

        if name not in VORONOI_KERNELS:
            continue
        cfg, path = kconfig(name + "_voronoi", {"cartesian_map": {"resolution": RES}})
        run(kiwi, "cartesian-map", path)
        vout = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "cartesian_map.dat")
        with open(vout, "rb") as f:
            assert f.read(8) == b"VGTCMAP1"
            res = struct.unpack("<i", f.read(4))[0]
            b = read_doubles(f, 6)
            arrays = [read_doubles(f, res ** 3) for _ in range(4)]
        ok &= verify(name + " voronoi", ((res, res, res), b, *arrays), kmass, 2e-2)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
