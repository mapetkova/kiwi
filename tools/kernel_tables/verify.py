"""Numerically verify coefficient tables against the brute-force reference (reference.py).

For random (r0, R0, phi) the table-based I_P (integral over phi' of  A * sum_i a_i mu^i q0^(1-i) + A*B,
with the region switching where the vertex radius crosses a breakpoint) is compared with the
reference integral, for (1) the tables derived here and (2) the thesis tables as printed.

Usage: python3 verify.py [n_points_per_kernel]
"""
import json, math, os, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kernels import KERNELS
from derive import derive_2d, derive_3d
from reference import integrate, pyramid_2d, pyramid_3d

def region(c, x):
    """1-based region of radius x: first k with x <= c_k, n+1 if x > 1."""
    for k, ck in enumerate(c, 1):
        if x <= ck:
            return k
    return len(c) + 1

def table_3d(a, B, c, A, r0, R0, phi, n=20):
    q0 = r0
    j = region(c, r0)
    cuts = [0.0] + [math.acos(R0 / math.sqrt(x * x - r0 * r0)) for x in c
                    if x * x - r0 * r0 > R0 * R0] + [phi]
    cuts = sorted(set(x for x in cuts if x <= phi))
    tot = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        k = region(c, math.sqrt(r0 * r0 + (R0 / math.cos(0.5 * (lo + hi))) ** 2))
        def g(ph):
            mu = r0 / math.sqrt(r0 * r0 + (R0 / math.cos(ph)) ** 2)
            s = sum(float(v) * mu ** i * q0 ** (1 - i) for i, v in a[k].items())
            return s + sum(float(v) * q0 ** p for p, v in B[(j, k)].items())
        tot += integrate(g, lo, hi, n)
    return A * tot

def table_2d(a, c, A, r0, phi, n=24):
    cuts = [0.0] + [math.acos(r0 / x) for x in c if x > r0] + [phi]
    cuts = sorted(set(x for x in cuts if x <= phi))
    tot = 0.0
    for lo, hi in zip(cuts, cuts[1:]):
        k = region(c, r0 / math.cos(0.5 * (lo + hi)))
        g = lambda ph: sum(float(v) * (r0 / math.cos(ph)) ** p for p, v in a[k].items())
        tot += integrate(g, lo, hi, n)
    return A * tot

TOL = 1e-9   # relative; the table form loses a few digits to cancellation at small r0

def main(npts):
    here = os.path.dirname(os.path.abspath(__file__))
    T = json.load(open(os.path.join(here, "thesis_tables.json")))
    P = lambda d: {int(k): Fr(v) for k, v in d.items()}
    rng = random.Random(1)
    worst = {}
    for key, K in KERNELS.items():
        c = [float(x) for x in K.breaks]
        A2, A3 = float(K.norm(2)) / math.pi, float(K.norm(3)) / math.pi
        d2 = derive_2d(K); a3, B3 = derive_3d(K)
        th2 = {int(k): P(v) for k, v in T["2d"][key].items()}
        th3a = {int(k): P(v) for k, v in T["3d_a"][key].items()}
        th3b = {tuple(map(int, k.split(","))): P(v) for k, v in T["3d_b"][key].items()}
        for _ in range(npts):
            r0 = rng.uniform(0.03, 1.2); R0 = rng.uniform(0.05, 1.6)
            phi = rng.uniform(0.2, 1.5)
            r2 = pyramid_2d(K, A2, r0, phi)
            r3 = pyramid_3d(K, A3, r0, R0, phi)
            res = {("2D", "derived"): table_2d(d2, c, A2, r0, phi),
                   ("2D", "thesis"): table_2d(th2, c, A2, r0, phi),
                   ("3D", "derived"): table_3d(a3, B3, c, A3, r0, R0, phi)}
            try:
                res[("3D", "thesis")] = table_3d(th3a, th3b, c, A3, r0, R0, phi)
            except KeyError:
                res[("3D", "thesis")] = float("nan")
            for (dim, src), v in res.items():
                ref = r2 if dim == "2D" else r3
                err = abs(v - ref) / max(abs(ref), 1e-300)
                w = worst.setdefault((key, dim, src), 0.0)
                worst[(key, dim, src)] = max(w, err) if err == err else float("inf")
    print(f"max relative error of table-based I_P vs brute-force quadrature ({npts} random points per kernel)")
    print(f"{'kernel':6s} {'2D derived':>12s} {'2D thesis':>12s} {'3D derived':>12s} {'3D thesis':>12s}")
    for key in KERNELS:
        print(f"{key:6s} {worst[(key,'2D','derived')]:12.2e} {worst[(key,'2D','thesis')]:12.2e} "
              f"{worst[(key,'3D','derived')]:12.2e} {worst[(key,'3D','thesis')]:12.2e}")
    bad = [k for k, v in worst.items() if k[2] == "derived" and not v < TOL]
    print("derived tables: " + ("OK (all within %.0e)" % TOL if not bad else "FAILED for %s" % bad))
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 8))
