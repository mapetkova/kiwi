"""Exact definitions of the piecewise-polynomial SPH kernels from thesis appendix A.

Every kernel is W(r) = A_d / H^d * w(q), q = r / H, with support radius H (w = 0 for q >= 1).
Polynomials are dicts {power: Fraction}; negative powers are allowed (Laurent polynomials).
Standard library only, all arithmetic exact.
"""
from fractions import Fraction as Fr
from math import pi

# ---------------------------------------------------------------- polynomial helpers

def padd(a, b, sb=1):
    out = dict(a)
    for p, c in b.items():
        out[p] = out.get(p, Fr(0)) + sb * c
    return {p: c for p, c in out.items() if c != 0}

def pmul(a, b):
    out = {}
    for p, c in a.items():
        for q, d in b.items():
            out[p + q] = out.get(p + q, Fr(0)) + c * d
    return {p: c for p, c in out.items() if c != 0}

def pscale(a, s):
    return {p: c * s for p, c in a.items() if c * s != 0}

def ppow(a, n):
    out = {0: Fr(1)}
    for _ in range(n):
        out = pmul(out, a)
    return out

def pint(a):
    """Antiderivative without constant term. A 1/x term would give a log: not supported."""
    out = {}
    for p, c in a.items():
        if p == -1:
            raise ValueError("log term in integral")
        out[p + 1] = c / (p + 1)
    return out

def peval(a, x):
    return sum((c * x ** p for p, c in a.items()), Fr(0) if isinstance(x, (int, Fr)) else 0.0)

# ---------------------------------------------------------------- kernel definitions

class Kernel:
    def __init__(self, key, name, family, terms=None, wend=None):
        self.key, self.name, self.family = key, name, family
        self.terms = terms          # splines: sum coef * (shift - q)_+^n
        self.wend = wend            # Wendland: (n, poly): (1-q)_+^n * poly(q)
        shifts = sorted({s for _, s, _ in terms}) if terms else [Fr(1)]
        self.breaks = shifts        # c_1 < ... < c_n = 1

    def w_direct(self, q):
        """Shape function straight from the definition (floats), independent of pieces()."""
        if q >= 1.0:
            return 0.0
        if self.terms:
            return sum(c * max(0.0, float(s) - q) ** n for c, s, n in self.terms)
        n, poly = self.wend
        return (1.0 - q) ** n * sum(float(c) * q ** p for p, c in poly.items())

    def pieces(self):
        """[(lo, hi, poly)] for the n pieces on (c_{k-1}, c_k)."""
        out, lo = [], Fr(0)
        for hi in self.breaks:
            if self.terms:
                poly = {}
                for c, s, n in self.terms:
                    if s >= hi:
                        poly = padd(poly, pscale(ppow({0: s, 1: Fr(-1)}, n), c))
            else:
                n, p = self.wend
                poly = pmul(ppow({0: Fr(1), 1: Fr(-1)}, n), p)
            out.append((lo, hi, poly))
            lo = hi
        return out

    def norm(self, dim):
        """Exact a_d such that A_d = a_d / pi (so that W integrates to 1)."""
        tot = Fr(0)
        for lo, hi, poly in self.pieces():
            P = pint(pmul(poly, {dim - 1: Fr(1)}))
            tot += peval(P, hi) - peval(P, lo)
        return 1 / ({2: 2, 3: 4}[dim] * tot)

def _w(*coefs):
    return {i: Fr(c) for i, c in enumerate(coefs)}

KERNELS = {k.key: k for k in [
    Kernel("M4", "Cubic spline", "bspline", terms=[(Fr(1), Fr(1), 3), (Fr(-4), Fr(1, 2), 3)]),
    Kernel("M5", "Quartic spline", "bspline",
           terms=[(Fr(1), Fr(1), 4), (Fr(-5), Fr(3, 5), 4), (Fr(10), Fr(1, 5), 4)]),
    Kernel("M6", "Quintic spline", "bspline",
           terms=[(Fr(1), Fr(1), 5), (Fr(-6), Fr(2, 3), 5), (Fr(15), Fr(1, 3), 5)]),
    Kernel("WC2", "Wendland C2", "wendland", wend=(4, _w(1, 4))),
    Kernel("WC4", "Wendland C4", "wendland", wend=(6, _w(1, 6, Fr(35, 3)))),
    Kernel("WC6", "Wendland C6", "wendland", wend=(8, _w(1, 8, 25, 32))),
]}
