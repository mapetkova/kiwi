#!/usr/bin/env python3
"""`kiwi --version` must equal the `version:` in CITATION.cff (the two are bumped by hand together).

Usage: version_consistent_test.py <path-to-kiwi-executable> <CITATION.cff>
"""
import re
import subprocess
import sys


def main():
    kiwi, cff = sys.argv[1], sys.argv[2]
    out = subprocess.run([kiwi, "--version"], capture_output=True, text=True).stdout.strip()
    m = re.search(r"^version:\s*['\"]?([^\s'\"]+)", open(cff).read(), re.M)
    if not m:
        print(f"FAIL: no top-level `version:` in {cff}")
        return 1
    ok = out == f"kiwi {m.group(1)}"
    print(f"  [{'PASS' if ok else 'FAIL'}] kiwi --version ({out!r}) vs CITATION.cff version ({m.group(1)!r})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
