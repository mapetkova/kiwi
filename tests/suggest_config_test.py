#!/usr/bin/env python3
"""tools/suggest_config.py on a Phantom binary dump (needs numpy; skipped with exit 77 without it).

The script reads the dump through `kiwi convert-snapshot`. Checks: the config it writes points at the
original dump with format phantom_binary, the suggested domain equals the one from the same particles as
plain ascii, that config runs (`kiwi point-set`), --shift-snapshot gives a plain_ascii copy centred on the
origin, a small dump needs --allow-missing-velocity (which lands in the config), and a missing kiwi
executable is a clear error.

Usage: suggest_config_test.py <path-to-kiwi-executable> [work_dir]
"""
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SCRIPT = os.path.join(HERE, "..", "tools", "suggest_config.py")


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def main():
    try:
        import numpy  # noqa: F401
    except ImportError:
        print("SKIP: numpy is not installed")
        return 77
    from phantom_dump_test import make_particles, write_case

    kiwi = os.path.abspath(sys.argv[1])
    work = sys.argv[2] if len(sys.argv) > 2 else tempfile.mkdtemp(prefix="kiwi_suggest_test_")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    rows = make_particles(random.Random(3), 150)
    pc_units = (3.086e18, 1.989e33, 3.086e18 / 1.0e5)
    dump = os.path.join(work, "snap.dump")
    write_case(dump, rows, units=pc_units)
    small = os.path.join(work, "small.dump")
    write_case(small, rows, rb=4, drop=("vx", "vy", "vz"), ident="ST:Phantom:test (hydro)", units=pc_units)
    txt = os.path.join(work, "snap.txt")
    subprocess.run([kiwi, "convert-snapshot", dump, txt], check=True, capture_output=True)

    def suggest(snapshot, tag, *extra, kiwi_path=kiwi):
        out = os.path.join(work, tag + ".json")
        r = subprocess.run([sys.executable, SCRIPT, snapshot, "--non-interactive", "--run-name", tag,
                            "--work-folder", os.path.join(work, "runs_" + tag), "--output", out,
                            "--kiwi", kiwi_path, *extra], capture_output=True, text=True, cwd=work)
        return r, (json.load(open(out)) if r.returncode == 0 and os.path.exists(out) else None)

    ok = True
    r, c_dump = suggest(dump, "dump")
    ok &= check("suggest_config.py on a Phantom dump succeeds", r.returncode == 0, (r.stderr.strip()[-200:] if r.returncode else ""))
    r2, c_txt = suggest(txt, "txt")
    ok &= check("... and on the same particles as plain ascii", r2.returncode == 0, (r2.stderr.strip()[-200:] if r2.returncode else ""))
    if c_dump and c_txt:
        ok &= check("config keeps the original dump and format phantom_binary",
                    c_dump["input_snapshot"] == dump and c_dump["input_snapshot_format"] == "phantom_binary")
        same = all(abs(c_dump["domain"][k] - c_txt["domain"][k]) < 1e-9 for k in c_dump["domain"])
        ok &= check("suggested domain equals the plain-ascii one", same)
        r = subprocess.run([kiwi, "point-set", os.path.join(work, "dump.json")], capture_output=True, text=True)
        ok &= check("the suggested config runs (kiwi point-set)", r.returncode == 0, (r.stderr.strip()[-200:] if r.returncode else ""))

    r, c = suggest(dump, "shift", "--shift-snapshot")
    ok &= check("--shift-snapshot succeeds", r.returncode == 0, (r.stderr.strip()[-200:] if r.returncode else ""))
    if c:
        shifted = c["input_snapshot"]
        ok &= check("... writes a plain_ascii copy and points the config at it",
                    c["input_snapshot_format"] == "plain_ascii" and os.path.isfile(shifted) and shifted != dump)
        mid = [(c["domain"][a + "_min"] + c["domain"][a + "_max"]) / 2 for a in "xyz"]
        ok &= check("... with a domain centred on the origin", all(abs(m) < 1e-9 for m in mid))

    r, _ = suggest(small, "small")
    ok &= check("a small dump without --allow-missing-velocity fails with kiwi's message",
                r.returncode == 1 and "failed to read snapshot" in r.stderr)
    r, c = suggest(small, "small_ok", "--allow-missing-velocity")
    ok &= check("with --allow-missing-velocity it works and the config says so",
                r.returncode == 0 and c is not None and c.get("allow_missing_velocity") is True)

    r, _ = suggest(dump, "nokiwi", kiwi_path=os.path.join(work, "no-such-kiwi"))
    ok &= check("a missing kiwi executable is a clear error",
                r.returncode == 1 and "kiwi executable was not found" in r.stderr)

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
