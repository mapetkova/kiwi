#!/usr/bin/env python3
"""POLARIS writer test for KIWI (standard library only).

Runs `kiwi write-polaris` on the synthetic cloud from synthetic_test.py (cubic domain, no buffer) and
decodes grid-polaris.dat with decode_polaris_grid.py, then checks it against the same run's mapping.dat:

  header      sim id, the fixed quantity list, cell count, box side length in metres;
  geometry    cell volumes (metres^3) add up to the box volume, centres lie inside the box, neighbour
              lists are symmetric and refer to existing rows (or to a domain-wall marker < 0);
  payload     per cell: uniform gas temperature, velocity and density equal mapping.dat's, scaled by the
              polaris_writer settings (float32 precision), rows in ascending cell id.

(validate_writer_polaris.py is the one-off cross-check against the legacy grid-file.cc output; it needs
that legacy file and is not run here.)

Usage: polaris_writer_test.py <path-to-kiwi-executable> [work_dir]
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
from decode_polaris_grid import decode  # noqa: E402
from synthetic_test import HALF_DOMAIN, write_snapshot  # noqa: E402

SCALES = {"length_scale_m": 2.5e16, "velocity_scale": 1500.0, "density_scale": 3.0e-18,
          "gas_temperature_k": 42.0, "sim_id": 7}


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_polaris_")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    snap = os.path.join(work, "snapshot.out")
    write_snapshot(snap)
    cfg = {"run_name": "polaris",
           "domain": {a + s: v * HALF_DOMAIN for a in "xyz" for s, v in (("_min", -1), ("_max", 1))},
           "input_snapshot": snap, "input_snapshot_format": "plain_ascii",
           "work_folder": os.path.join(work, "runs"),
           "buffer": {"enabled": False}, "lloyd": {"iterations": 0},
           "polaris_writer": SCALES}
    cfg_path = os.path.join(work, "polaris.json")
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=1)
    r = subprocess.run([kiwi, "write-polaris", cfg_path], capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: kiwi write-polaris exited with {r.returncode}\n{r.stderr[-2000:]}")
        return 1
    run_dir = os.path.join(cfg["work_folder"], os.listdir(cfg["work_folder"])[0])
    g = decode(os.path.join(run_dir, "grid-polaris.dat"))
    cells, neighbors = g["cells"], g["neighbors"]

    mapping = {}
    with open(os.path.join(run_dir, "mapping.dat")) as f:  # id x y z volume density vx vy vz
        for line in f:
            t = line.split()
            mapping[int(t[0])] = [float(x) for x in t[4:9]]
    ids = sorted(mapping)

    ok = True
    L = SCALES["length_scale_m"]
    side = 2 * HALF_DOMAIN * L
    print("header")
    ok &= check("sim id and quantity list", g["sim"] == SCALES["sim_id"] and g["params"] == 5
                and list(g["plist"]) == [3, 7, 8, 9, 28], f"sim={g['sim']} plist={list(g['plist'])}")
    ok &= check("cell count matches mapping.dat", g["cells_d"] == len(cells) == len(ids) > 0, f"{len(cells)} cells")
    ok &= check("box side length (m)", abs(g["length_scale"] - side) < 1e-9 * side, f"{g['length_scale']:.6e}")

    print("geometry")
    vol = sum(c[3] for c in cells)
    ok &= check("cell volumes fill the box", abs(vol - side ** 3) < 1e-6 * side ** 3,
                f"rel.diff={abs(vol - side ** 3) / side ** 3:.1e}")
    inside = all(abs(c[k]) <= 0.5 * side * (1 + 1e-6) for c in cells for k in range(3))
    ok &= check("cell centres lie inside the box", inside)
    n = len(cells)
    in_range = all(-n < j < n for nb in neighbors for j in nb)
    ok &= check("neighbour ids are rows or wall markers", in_range)
    sets = [set(j for j in nb if j >= 0) for nb in neighbors]
    asym = sum(1 for i, s in enumerate(sets) for j in s if i not in sets[j])
    ok &= check("neighbour lists are symmetric", asym == 0, f"{asym} one-sided pairs")
    ok &= check("no cell lists itself or a neighbour twice",
                all(i not in s for i, s in enumerate(sets)) and all(len(set(nb)) == len(nb) for nb in neighbors))

    print("payload (vs mapping.dat, scaled)")
    worst_t = worst_d = worst_v = 0.0
    for cid, c in zip(ids, cells):
        _vol, dens, vx, vy, vz = mapping[cid]
        temp, wvx, wvy, wvz, wdens = c[4]
        worst_t = max(worst_t, abs(temp - SCALES["gas_temperature_k"]))
        want = dens * SCALES["density_scale"]
        worst_d = max(worst_d, abs(wdens - want) / want if want else abs(wdens))
        for w, v in ((wvx, vx), (wvy, vy), (wvz, vz)):
            worst_v = max(worst_v, abs(w - v * SCALES["velocity_scale"]) / SCALES["velocity_scale"])
    ok &= check("gas temperature", worst_t < 1e-4, f"max |dT|={worst_t:.1e} K")
    ok &= check("density (float32)", worst_d < 1e-6, f"max rel.err={worst_d:.1e}")
    ok &= check("velocity (float32)", worst_v < 1e-6, f"max err/scale={worst_v:.1e}")
    vol_m = max(abs(c[3] - mapping[cid][0] * L ** 3) / c[3] for cid, c in zip(ids, cells))
    ok &= check("cell volume vs mapping.dat (m^3)", vol_m < 1e-9, f"max rel.err={vol_m:.1e}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
