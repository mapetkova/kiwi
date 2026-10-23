#!/usr/bin/env python3
"""Bad-input tests for KIWI (standard library only).

Every malformed config or snapshot below must make `kiwi` exit with status 1 and a one-line
"[kiwi] error: ..." message that names the problem -- never a crash (signal), a hang, or a silent
exit 0 on data that was only partly read.

Usage: bad_input_test.py <path-to-kiwi-executable> [work_dir]
Exit code 0 on success, 1 on any failure (this is what ctest looks at).
"""
import copy
import json
import os
import shutil
import subprocess
import sys
import tempfile

GOOD_ROWS = ["0.1 0.1 0.1 0 0 0 1 0.3 1", "0.2 -0.3 0.1 0 0 0 1 0.3 1", "-0.5 0.4 0.3 0 0 0 1 0.3 1"]
BASE = {
    "run_name": "bad",
    "domain": {"x_min": -1, "x_max": 1, "y_min": -1, "y_max": 1, "z_min": -1, "z_max": 1},
    "input_snapshot": "good.out",
    "input_snapshot_format": "dale_phantom",
    "work_folder": "work",
}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    kiwi = os.path.abspath(sys.argv[1])
    root = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_bad_")
    shutil.rmtree(root, ignore_errors=True)
    os.makedirs(root)
    os.chdir(root)  # relative paths in the configs below resolve here

    def write(name, text):
        with open(name, "w") as f:
            f.write(text)

    write("good.out", "\n".join(GOOD_ROWS) + "\n")
    write("empty.out", "")
    write("garbage.out", "a b c\n1 2\n")
    write("nan_row.out", "\n".join([GOOD_ROWS[0], "0.1 0.1 nan 0 0 0 1 0.3 1"] + GOOD_ROWS[1:]) + "\n")
    write("text_row.out", "\n".join([GOOD_ROWS[0], "0.1 0.1 abc 0 0 0 1 0.3 1"] + GOOD_ROWS[1:]) + "\n")
    write("short_row.out", "\n".join([GOOD_ROWS[0], "0.1 0.1"] + GOOD_ROWS[1:]) + "\n")
    write("long_row.out", "\n".join([GOOD_ROWS[0] + " 7"] + GOOD_ROWS[1:]) + "\n")
    write("no_header_cols.out", "# x y z m vx vy vz\n0.1 0.1 0.1 1 0 0 0\n")
    with open("random.bin", "wb") as f:
        f.write(os.urandom(300))

    def cfg(**changes):
        c = copy.deepcopy(BASE)
        for k, v in changes.items():
            if isinstance(v, dict) and isinstance(c.get(k), dict):
                c[k].update(v)
            else:
                c[k] = v
        return c

    # (name, stage, config dict or raw text, expected substring of the error message)
    cases = [
        ("missing snapshot", "sph-to-grid", cfg(input_snapshot="nope.out"), "could not open SPH snapshot"),
        ("empty snapshot", "sph-to-grid", cfg(input_snapshot="empty.out"), "no rows"),
        ("text instead of data", "sph-to-grid", cfg(input_snapshot="garbage.out"), "line 1"),
        ("NaN in a row", "sph-to-grid", cfg(input_snapshot="nan_row.out"), "line 2"),
        ("non-number in a row", "sph-to-grid", cfg(input_snapshot="text_row.out"), "line 2"),
        ("row with too few columns", "sph-to-grid", cfg(input_snapshot="short_row.out"), "line 2"),
        ("row with too many columns", "sph-to-grid", cfg(input_snapshot="long_row.out"), "more than 9 columns"),
        ("splash file without data", "sph-to-grid",
         cfg(input_snapshot="no_header_cols.out", input_snapshot_format="splash_ascii"), "header preamble"),
        ("random bytes as Phantom dump", "sph-to-grid",
         cfg(input_snapshot="random.bin", input_snapshot_format="phantom_binary"), "capture pattern"),
        ("unknown snapshot format", "sph-to-grid", cfg(input_snapshot_format="fits"), "input_snapshot_format"),
        ("unknown kernel", "sph-to-grid", cfg(kernel={"name": "gaussian"}), "unknown kernel"),
        ("zero resolution", "sph-to-grid", cfg(cartesian_grid={"resolution": 0}), "resolution must be > 0"),
        ("negative resolution", "sph-to-grid", cfg(cartesian_grid={"resolution": -4}), "resolution must be > 0"),
        ("negative Lloyd iterations", "point-set", cfg(lloyd={"iterations": -1}), "iterations must be >= 0"),
        ("inverted domain", "sph-to-grid",
         cfg(domain={"x_min": 1, "x_max": -1}), "x_min must be < domain.x_max"),
        ("domain value is a string", "sph-to-grid", cfg(domain={"x_min": "a"}), "type must be number"),
        ("missing domain key", "sph-to-grid", {"run_name": "x", "domain": {"x_min": -1}}, "missing required field"),
        ("invalid JSON", "sph-to-grid", "{not json", "not valid JSON"),
        ("missing config file", "sph-to-grid", None, "could not open config file"),
    ]

    ok = True

    def check(name, good, detail=""):
        nonlocal ok
        print(f"  [{'PASS' if good else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
        ok &= good

    def run(args):
        shutil.rmtree("work", ignore_errors=True)  # cached stage output must not hide an error
        try:
            return subprocess.run([kiwi] + args, capture_output=True, text=True, timeout=120)
        except subprocess.TimeoutExpired:
            return None

    for i, (name, stage, c, expect) in enumerate(cases):
        path = f"case{i}.json"
        if c is not None:
            write(path, c if isinstance(c, str) else json.dumps(c))
        r = run([stage, path])
        if r is None:
            check(name, False, "timed out")
            continue
        err = [l for l in r.stderr.splitlines() if l.startswith("[kiwi] error")]
        check(name, r.returncode == 1 and len(err) == 1 and expect in err[0],
              f"exit={r.returncode}, " + (err[0][:110] if err else "no [kiwi] error line"))

    print("command line")
    r = run([])
    check("no arguments: usage, exit 1", r.returncode == 1 and "usage" in r.stderr)
    write("good.json", json.dumps(cfg()))
    r = run(["no-such-stage", "good.json"])
    check("unknown stage: exit 1", r.returncode == 1 and "unknown stage" in r.stderr)
    r = run(["sph-to-grid", "good.json"])
    check("the unmodified good config runs", r.returncode == 0, f"exit={r.returncode}")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
