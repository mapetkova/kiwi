"""Derive the I_P coefficient tables of thesis appendix A from the kernel definition alone.

Setup (H = 1, kernel shape w(q); the overall factor A_d is divided out, as in the thesis tables).
A "vertex pyramid" has its apex at the particle and a triangular base in a wall plane at distance
r0, spanned by the foot of the perpendicular, the foot on the wall edge (distance R0 from it) and
the vertex (azimuth phi).  Its kernel integral is I_P(phi), and dI_P/dphi is a Laurent polynomial
in mu = r0 / r (r = distance of the vertex from the particle).

 3D:  dI_P/dphi / A = sum_i a_i mu^i q0^(1-i)  +  B(q0)          (i <= 1, i != 0)
      a_i  : per vertex region k (q = r in (c_{k-1}, c_k]);
      B    : per (region j of q0, vertex region k >= j), a polynomial in q0 = r0.
 2D:  dI_P/dphi / A = sum_p a_p q0^p sec(phi)^p                   (p >= 0), per vertex region k.

Regions are numbered 1..n for the pieces and n+1 for q > 1 (where the kernel vanishes).
"""
from fractions import Fraction as Fr
from kernels import padd, pmul, pscale, pint, peval

def _mhat(kernel, dim):
    """Enclosed-mass shape function m_k(rho) = int_0^rho w s^(dim-1) ds, per region (poly, const-incl.)."""
    pieces = kernel.pieces()
    ms, prev_val, prev_c = [], Fr(0), Fr(0)
    for lo, hi, poly in pieces:
        P = pint(pmul(poly, {dim - 1: Fr(1)}))
        const = prev_val - peval(P, lo)
        m = dict(P)
        if const != 0:
            m[0] = m.get(0, Fr(0)) + const
        ms.append(m)
        prev_val = peval(m, hi)
    ms.append({0: prev_val} if prev_val != 0 else {})     # outside the support
    return ms

def derive_2d(kernel):
    """{k: {p: a_p}}, k = 1..n+1."""
    return {k + 1: dict(m) for k, m in enumerate(_mhat(kernel, 2))}

def derive_3d(kernel):
    """Returns (a, B): a = {k: {i: a_i}} (i != 0), B = {(j, k): {power: b}} with j <= k."""
    ms = _mhat(kernel, 3)
    nreg = len(ms)                                   # n + 1
    cs = [Fr(0)] + list(kernel.breaks)               # c_0 = 0, c_1 .. c_n = 1
    F = [pint(pmul(m, {-2: Fr(1)})) for m in ms]     # F_k(rho) = int m_k / rho^2
    a = {}
    for k in range(1, nreg + 1):
        a[k] = {-p: c for p, c in F[k - 1].items()}          # r0 * rho^p  ->  a_{-p}
    def Fv(k, x):                                     # F_k evaluated at x (exact)
        return peval(F[k - 1], x)
    B = {}
    for j in range(1, nreg + 1):
        for k in range(j, nreg + 1):
            q0 = {1: Fr(1)}
            # constant part of r0 * (F_j(c_j) - F_j(r0) + middle regions - F_k(c_{k-1})), k > j
            if k == j:
                b = pscale(pmul(q0, F[j - 1]), Fr(-1))
            else:
                const = Fv(j, cs[j]) if j < nreg else Fr(0)
                for l in range(j + 1, k):
                    const += Fv(l, cs[l]) - Fv(l, cs[l - 1])
                const -= Fv(k, cs[k - 1])
                b = padd(pscale(pmul(q0, F[j - 1]), Fr(-1)), {1: const})
            B[(j, k)] = {p: c for p, c in b.items() if c != 0}
    return a, B
