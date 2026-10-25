#!/usr/bin/env python3
"""`kiwi convert-snapshot` test (standard library only): synthetic Phantom dumps in several layouts are
converted to plain_ascii and must give exactly the particles the dump holds (gas only, dead particles and
other particle types dropped), in the units the reader promises (pc, Msun, km/s; rho in g/cm^3), with rho = m (hfact/h)^3 in code units.

Usage: convert_snapshot_test.py <path-to-kiwi-executable> [work_dir]
Exit code 0 on success, 1 on any failure (this is what ctest looks at).
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from phantom_dump_test import make_particles, write_case  # noqa: E402


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_convert_")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    rows = make_particles(random.Random(3), 150)
    # code units chosen so that the reader's default units (pc, Msun, km/s) are the identity
    pc_units = (3.086e18, 1.989e33, 3.086e18 / 1.0e5)

    def convert(tag, dump_kwargs, extra=(), units=(1.0, 1.0, 1.0)):
        dump = os.path.join(work, tag + ".dump")
        keep, m_code, hfact = write_case(dump, rows, **dump_kwargs)
        u = dump_kwargs.get("units", (1.0, 1.0, 1.0))
        rho_cgs = u[1] / u[0] ** 3  # code density -> g/cm^3
        out = os.path.join(work, tag + ".txt")
        r = subprocess.run([kiwi, "convert-snapshot", dump, out, *extra], capture_output=True, text=True)
        if r.returncode != 0:
            return check(f"{tag}: convert-snapshot exits 0", False, r.stderr.strip()[-200:])
        got = [[float(v) for v in line.split()] for line in open(out)]
        ok = check(f"{tag}: {len(keep)} gas particles written", len(got) == len(keep), f"got {len(got)}")
        if len(got) != len(keep):
            return False
        worst = 0.0
        for g, k in zip(got, keep):
            h = k[6]
            want = k[:6] + [m_code, h, m_code * (hfact / h) ** 3 * rho_cgs]
            if "drop" in dump_kwargs:
                want[3:6] = [0.0, 0.0, 0.0]
            worst = max(worst, max(abs(a - b) / max(abs(b), 1e-300) if b else abs(a) for a, b in zip(g, want)))
        ok &= check(f"{tag}: columns x y z vx vy vz m h rho match the dump", worst < 1e-9, f"max rel.diff={worst:.1e}")
        return ok

    ok = True
    ident = ("--length-cm", "1", "--mass-g", "1", "--velocity-cms", "1")
    ok &= convert("units_override", {}, ident)
    ok &= convert("default_units_pc", {"units": pc_units})
    ok &= convert("big_endian_dust_dropped", {"e": ">", "dust_every": 4, "units": pc_units})
    ok &= convert("int64_real32_mpi", {"ib": 8, "rb": 4, "mpi": 3, "units": pc_units})
    ok &= convert("small_dump_no_velocity",
                  {"rb": 4, "drop": ("vx", "vy", "vz"), "ident": "ST:Phantom:test (hydro)", "units": pc_units},
                  ("--allow-missing-velocity",))

    r = subprocess.run([kiwi, "convert-snapshot", os.path.join(work, "small_dump_no_velocity.dump"),
                        os.path.join(work, "x.txt")], capture_output=True, text=True)
    ok &= check("a small dump without --allow-missing-velocity is refused", r.returncode == 1 and "[kiwi] error" in r.stderr)
    r = subprocess.run([kiwi, "convert-snapshot", "a"], capture_output=True, text=True)
    ok &= check("missing arguments: usage, exit 1", r.returncode == 1 and "usage" in r.stderr)
    r = subprocess.run([kiwi, "convert-snapshot", "a", "b", "--nope"], capture_output=True, text=True)
    ok &= check("unknown option: exit 1", r.returncode == 1 and "unknown option" in r.stderr)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
