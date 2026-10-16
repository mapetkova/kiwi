"""Brute-force numerical vertex-pyramid kernel integrals (the independent reference).

I_P(phi) = integral of W over the pyramid with apex at the particle and a triangular base in a
wall plane at distance r0: foot of the perpendicular O', foot P on the wall edge (|O'P| = R0),
vertex Q at azimuth phi (2D: the segment from the foot to the vertex).  Computed with Gauss-
Legendre quadrature of the kernel as *defined* (Kernel.w_direct), not from any derived formula.
"""
import math

_GL = {}
def gl(n):
    """Gauss-Legendre nodes/weights on [-1, 1] (Newton iteration)."""
    if n in _GL:
        return _GL[n]
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        xs.append(x); ws.append(2 / ((1 - x * x) * dp * dp))
    _GL[n] = (xs, ws)
    return _GL[n]

def integrate(f, a, b, n):
    if b <= a:
        return 0.0
    xs, ws = gl(n)
    h, m = 0.5 * (b - a), 0.5 * (b + a)
    return h * sum(w * f(m + h * x) for x, w in zip(xs, ws))

def enclosed(kernel, A, dim, rho):
    """m(rho) = int_0^rho W(r) r^(dim-1) dr, W = A * w(r) (H = 1)."""
    cuts = [0.0] + [float(c) for c in kernel.breaks if float(c) < rho] + [rho]
    cuts = sorted(set(min(c, 1.0) for c in cuts))
    return A * sum(integrate(lambda r: kernel.w_direct(r) * r ** (dim - 1), cuts[i], cuts[i + 1], 12)
                   for i in range(len(cuts) - 1))

def pyramid_3d(kernel, A, r0, R0, phi, n=20):
    """A = a_3 / pi (float).  Integral over the vertex pyramid; the wedge at azimuth phi' has the
    kernel mass  int_0^{Rmax} m(rho) r0 R / rho^3 dR,  rho = sqrt(r0^2+R^2)."""
    cuts = [0.0] + [math.acos(R0 / math.sqrt(c * c - r0 * r0)) for c in map(float, kernel.breaks)
                    if c * c - r0 * r0 > R0 * R0] + [phi]
    cuts = sorted(set(c for c in cuts if c <= phi))
    def wedge(ph):
        Rmax = R0 / math.cos(ph)
        rc = [0.0] + [math.sqrt(c * c - r0 * r0) for c in map(float, kernel.breaks)
                      if c * c > r0 * r0 and math.sqrt(c * c - r0 * r0) < Rmax] + [Rmax]
        return sum(integrate(lambda R: enclosed(kernel, A, 3, math.sqrt(r0 * r0 + R * R))
                             * r0 * R / (r0 * r0 + R * R) ** 1.5, rc[i], rc[i + 1], 16)
                   for i in range(len(rc) - 1))
    return sum(integrate(wedge, cuts[i], cuts[i + 1], n) for i in range(len(cuts) - 1))

def pyramid_2d(kernel, A, r0, phi, n=24):
    """A = a_2 / pi.  Integral over the triangle apex-foot-vertex: int_0^phi m(r0/cos phi') dphi'."""
    cuts = [0.0] + [math.acos(r0 / c) for c in map(float, kernel.breaks) if c > r0] + [phi]
    cuts = sorted(set(c for c in cuts if c <= phi))
    return sum(integrate(lambda p: enclosed(kernel, A, 2, r0 / math.cos(p)), cuts[i], cuts[i + 1], n)
               for i in range(len(cuts) - 1))
