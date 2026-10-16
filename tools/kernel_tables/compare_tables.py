"""Compare the thesis appendix-A coefficient tables with tables derived from the kernel definitions.

Usage: python3 compare_tables.py            (needs thesis_tables.json from parse_thesis.py)
Prints every disagreement (kernel, table, entry, thesis value, derived value); exit status 1 if any.
"""
import json, os, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import KERNELS
from derive import derive_2d, derive_3d

def diff(where, thesis, derived, out):
    for key in sorted(set(thesis) | set(derived), key=lambda s: (str(s))):
        t = thesis.get(key, Fr(0)); d = derived.get(key, Fr(0))
        if t != d:
            out.append((where, key, t, d))

def main():
    here = os.path.dirname(os.path.abspath(__file__))
    T = json.load(open(os.path.join(here, "thesis_tables.json")))
    P = lambda dct: {int(k): Fr(v) for k, v in dct.items()}
    bad, nchecked = [], 0
    for key, K in KERNELS.items():
        for dim in (2, 3):
            derived_A = K.norm(dim)
            if Fr(T["A"][key][str(dim)]) != derived_A:
                bad.append((f"{key} A_{dim}D", "", Fr(T["A"][key][str(dim)]), derived_A))
        for k, d in derive_2d(K).items():
            diff(f"{key} 2D region {k}", P(T["2d"][key].get(str(k), {})), d, bad); nchecked += 1
        a, B = derive_3d(K)
        for k, d in a.items():
            diff(f"{key} 3D a_i region {k}", P(T["3d_a"][key].get(str(k), {})), d, bad); nchecked += 1
        for (j, k), d in B.items():
            tb = T["3d_b"][key].get(f"{j},{k}")
            if tb is None:
                bad.append((f"{key} 3D B_{k} (q0 region {j})", "row missing", "", d)); continue
            diff(f"{key} 3D B_{k} (q0 region {j})", P(tb), d, bad); nchecked += 1
        for jk in T["3d_b"][key]:
            if tuple(map(int, jk.split(","))) not in B:
                bad.append((f"{key} 3D B row {jk}", "not derivable", "", ""))
    print(f"compared {nchecked} coefficient rows for {len(KERNELS)} kernels")
    for where, key, t, d in bad:
        print(f"MISMATCH  {where:34s} entry {str(key):>4s}: thesis {str(t):>14s}   derived {str(d)}")
    print("all tables agree" if not bad else f"{len(bad)} mismatch(es)")
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main())
