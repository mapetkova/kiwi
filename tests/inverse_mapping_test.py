#!/usr/bin/env python3
"""`kiwi inverse-mapping-test` (experimental; a first step towards mapping a grid back onto particles).

On a synthetic cloud it maps the particles onto the Voronoi grid, then maps the cell masses back onto the
same particles, and checks the invariants the stage guarantees (see mapping.hpp for what it cannot detect):

  * domain larger than the cloud: every particle's kernel is fully inside, so every particle recovers
    its own mass (coverage_fraction == 1);
  * domain that cuts the cloud: particles outside the domain recover nothing (they have no cell), particles
    near the walls recover less than their own mass (coverage < 1), and none recovers more;
  * in both cases the total recovered mass equals the total mass on the grid (exactly, by construction),
    and the original total is the mass of all the particles.

Standard library only.

Usage: inverse_mapping_test.py <path-to-kiwi-executable> [work_dir]
Exit code 0 on success, 1 on any failure (this is what ctest looks at).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from synthetic_test import write_snapshot  # noqa: E402


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_inverse_")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    snap = os.path.join(work, "snapshot.out")
    total_mass = write_snapshot(snap)
    inside_2 = []  # particles inside the tight domain
    for line in open(snap):
        v = [float(t) for t in line.split()]
        if all(abs(c) < 2.0 for c in v[:3]):
            inside_2.append(v[6])

    def run(tag, half):
        cfg = {"run_name": tag, "domain": {a + b: s * half for a in "xyz" for b, s in (("_min", -1), ("_max", 1))},
               "input_snapshot": snap, "input_snapshot_format": "plain_ascii",
               "work_folder": os.path.join(work, tag), "buffer": {"enabled": False}, "lloyd": {"iterations": 0}}
        path = os.path.join(work, tag + ".json")
        json.dump(cfg, open(path, "w"))
        r = subprocess.run([kiwi, "inverse-mapping-test", path], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"FAIL: kiwi inverse-mapping-test exited with {r.returncode}\n{r.stderr[-1500:]}")
            sys.exit(1)
        out = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0], "inverse-mapping-test.dat")
        lines = [l.split() for l in open(out)]
        totals = [float(t) for t in lines[1]]  # original, grid, recovered
        rows = [(float(t[1]), float(t[2]), float(t[3])) for t in lines[2:]]
        return totals, rows

    ok = True
    print("domain larger than the cloud")
    (t_orig, t_grid, t_rec), rows = run("wide", 6.0)
    cov = [r[2] for r in rows]
    ok &= check("every particle is listed", len(rows) == 400, f"{len(rows)} rows")
    ok &= check("every particle recovers its own mass", all(abs(c - 1.0) < 1e-8 for c in cov),
                f"coverage in [{min(cov):.10f}, {max(cov):.10f}]")
    ok &= check("recovered total == grid total", abs(t_rec - t_grid) < 1e-9 * t_grid)
    ok &= check("grid total == particle mass", abs(t_grid - total_mass) < 1e-9 * total_mass,
                f"{t_grid:.12g} vs {total_mass:.12g}")

    print("domain that cuts the cloud")
    (t_orig, t_grid, t_rec), rows = run("tight", 2.0)
    cov = [r[2] for r in rows]
    n_outside = 400 - len(inside_2)
    n_zero = sum(1 for c in cov if c == 0.0)
    partial = [c for c in cov if 0.0 < c < 0.99]
    ok &= check("particles outside the domain recover nothing", n_outside > 0 and n_zero == n_outside,
                f"{n_zero} with coverage 0, {n_outside} outside")
    ok &= check("particles near the walls recover only part of their mass", len(partial) > 0,
                f"{len(partial)} with 0 < coverage < 0.99")
    ok &= check("no particle recovers more than its mass", max(cov) < 1.0 + 1e-8, f"max {max(cov):.10f}")
    ok &= check("recovered total == grid total", abs(t_rec - t_grid) < 1e-9 * t_grid)
    ok &= check("original total is the mass of all particles", abs(t_orig - total_mass) < 1e-9 * total_mass)
    ok &= check("less mass is recovered than the particles hold", t_rec < t_orig - 0.01 * total_mass,
                f"{t_rec:.4g} of {t_orig:.4g}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
