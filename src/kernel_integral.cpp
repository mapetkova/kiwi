#include "kiwi/kernel_integral.hpp"
#include "kiwi/kernel_integral_generic.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace kiwi {

double DensityLookupTable::at(int i, int j, int k) const {
    return values[(static_cast<size_t>(i) * (nR_0 + 1) + j) * (nphi + 1) + k];
}
double& DensityLookupTable::at(int i, int j, int k) {
    return values[(static_cast<size_t>(i) * (nR_0 + 1) + j) * (nphi + 1) + k];
}

DensityLookupTable load_density_table(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open density lookup table: " + path);

    DensityLookupTable t;
    f >> t.nr01 >> t.nR_01 >> t.nphi1;
    f >> t.rl >> t.mul >> t.cphil;
    f >> t.nr02 >> t.nR_02 >> t.nphi2;
    if (!f) throw std::runtime_error("density lookup table header is malformed: " + path);

    t.nr0 = t.nr01 + t.nr02;
    t.nR_0 = t.nR_01 + t.nR_02;
    t.nphi = t.nphi1 + t.nphi2;

    t.values.assign(static_cast<size_t>(t.nr0 + 2) * (t.nR_0 + 1) * (t.nphi + 1), 0.0);

    for (int i = 1; i <= t.nr0; ++i) {
        for (int j = 0; j < t.nR_0; ++j) {
            for (int k = 0; k < t.nphi; ++k) {
                if (!(f >> t.at(i, j, k)))
                    throw std::runtime_error("density lookup table is truncated (short read): " + path);
                t.at(0, j, k) = 0.0;
            }
        }
    }

    // Pad one extra shell in r0 by repeating the outermost table row -- the
    // interpolation below reads at i+1 of the last valid i.
    for (int j = 0; j < t.nR_0; ++j) {
        for (int k = 0; k < t.nphi; ++k) {
            t.at(t.nr0 + 1, j, k) = t.at(t.nr0, j, k);
        }
    }

    return t;
}

DensityLookupTable generate_density_table() {
    // Canonical grid shape, matching density-map-grid-250-300-300.dat's own
    // header exactly (50 150 150 / 0.1 0.98 0.98 / 199 150 150) -- see the
    // header comment on the declaration for provenance.
    DensityLookupTable t;
    t.nr01 = 50;
    t.nr02 = 199;
    t.nR_01 = 150;
    t.nR_02 = 150;
    t.nphi1 = 150;
    t.nphi2 = 150;
    t.rl = 0.1;
    t.mul = 0.98;
    t.cphil = 0.98;
    t.nr0 = t.nr01 + t.nr02;
    t.nR_0 = t.nR_01 + t.nR_02;
    t.nphi = t.nphi1 + t.nphi2;

    const double h = 1.0; // the table is built dimensionless, in units of h -- see gridded_integral()

    t.values.assign(static_cast<size_t>(t.nr0 + 2) * (t.nR_0 + 1) * (t.nphi + 1), 0.0);
    // Row 0 (i=0) is left at its assigned 0.0 default, same as
    // load_density_table()'s explicit t.at(0, j, k) = 0.0 -- gridding()
    // itself never writes this row either, only i=1..nr0.

    // n here is gridding()'s single shared sub-grid size (nR_01 == nR_02 ==
    // nphi1 == nphi2 == 150 in the canonical grid) -- kept as one local to
    // mirror the source this is ported from, not implying it could differ
    // per axis; if it ever needs to, the four t.nR_0*/t.nphi* fields above
    // already carry the real (currently-equal) values independently.
    const int n = t.nR_01;

    for (int i = 1; i <= t.nr0; ++i) {
        const double r0 = (i <= t.nr01) ? (t.rl / t.nr01) * i * h
                                         : t.rl + ((2.0 - t.rl) / t.nr02) * (i - t.nr01) * h;

        for (int j = 0; j < t.nR_0; ++j) {
            const double mu0 = (j < n) ? (t.mul / (n - 1)) * j
                                        : t.mul + ((1.0 - t.mul) / n) * (j - n + 1);
            // j == 0 -> mu0 == 0.0 -> R_0 == +inf (IEEE-754 double
            // division) here, deliberately not special-cased -- see the
            // header comment on generate_density_table() for why this
            // still produces the correct finite result.
            const double R_0 = r0 * std::sqrt(1.0 - mu0 * mu0) / mu0;

            for (int k = 0; k < t.nphi; ++k) {
                const double phi = (k < n) ? std::acos((t.cphil / (n - 1)) * k)
                                            : std::acos(t.cphil + ((1.0 - t.cphil) / n) * (k - n + 1));
                t.at(i, j, k) = full_integral_new(h, phi, r0, R_0);
            }
        }
    }

    // Pad one extra shell in r0 by repeating the outermost table row, same
    // as load_density_table() does after reading a file from disk -- the
    // interpolation in gridded_integral() reads at i+1 of the last valid i.
    for (int j = 0; j < t.nR_0; ++j) {
        for (int k = 0; k < t.nphi; ++k) {
            t.at(t.nr0 + 1, j, k) = t.at(t.nr0, j, k);
        }
    }

    return t;
}

namespace {
// Clamp an interpolation index into the table's valid range instead of
// reading out of bounds -- the legacy code relied on a fixed-size C array
// generously over-allocated beyond the nominal table dimensions to survive
// an occasional out-of-range i/j/k without crashing (it printed a warning
// but didn't clamp). We clamp explicitly and keep the warning.
int clamp_index(int v, int lo, int hi, const char* name) {
    if (v < lo || v > hi) {
        std::cerr << "[kiwi] warning: gridded_integral index " << name << "=" << v
                  << " out of range [" << lo << ", " << hi << "], clamping\n";
        if (v < lo) v = lo;
        if (v > hi) v = hi;
    }
    return v;
}
} // namespace

double gridded_integral(const DensityLookupTable& t, double h, double phi, double r0_old, double R_0_old) {
    if (r0_old == 0.0) return 0.0;
    if (R_0_old == 0.0) return 0.0;
    if (phi == 0.0) return 0.0;

    double r0 = r0_old / h;
    double R_0 = R_0_old / h;
    double cphi = std::cos(phi);

    if (r0 < 0.0) r0 = 0.0;
    if (R_0 < 0.0) R_0 = 0.0;

    double mu0 = r0 / std::sqrt(r0 * r0 + R_0 * R_0);
    if (r0 > 2.0) r0 = 2.0;

    const int nr01 = t.nr01, nr02 = t.nr02, nR_01 = t.nR_01, nR_02 = t.nR_02, nphi1 = t.nphi1, nphi2 = t.nphi2;
    const double rl = t.rl, mul = t.mul, cphil = t.cphil;

    int i, j, k;
    double frac;

    if (r0 < rl) i = static_cast<int>(r0 / (rl / nr01));
    else i = nr01 + static_cast<int>((r0 - rl) / ((2.0 - rl) / nr02));

    if (mu0 < mul) j = static_cast<int>(mu0 / (mul / (nR_01 - 1)));
    else j = nR_01 - 1 + static_cast<int>((mu0 - mul) / ((1.0 - mul) / nR_02));

    if (cphi < cphil) k = static_cast<int>(cphi / (cphil / (nphi1 - 1)));
    else k = nphi1 - 1 + static_cast<int>((cphi - cphil) / ((1.0 - cphil) / nphi2));

    i = clamp_index(i, 0, t.nr0, "i");
    j = clamp_index(j, 0, t.nR_0 - 1, "j");
    k = clamp_index(k, 0, t.nphi - 1, "k");

    if (r0 < rl) frac = (r0 - (rl / nr01) * i) / (rl / nr01);
    else frac = (r0 - rl - ((2.0 - rl) / nr02) * (i - nr01)) / ((2.0 - rl) / nr02);

    double fx1 = frac * t.at(i + 1, j, k) + (1 - frac) * t.at(i, j, k);
    double fx2 = frac * t.at(i + 1, j + 1, k) + (1 - frac) * t.at(i, j + 1, k);
    double fx3 = frac * t.at(i + 1, j, k + 1) + (1 - frac) * t.at(i, j, k + 1);
    double fx4 = frac * t.at(i + 1, j + 1, k + 1) + (1 - frac) * t.at(i, j + 1, k + 1);

    if (mu0 < mul) frac = (mu0 - mul / (nR_01 - 1) * j) / (mul / (nR_01 - 1));
    else frac = (mu0 - mul - (1.0 - mul) / nR_02 * (j - nR_01 + 1)) / ((1.0 - mul) / nR_02);
    double fy1 = frac * fx2 + (1 - frac) * fx1;
    double fy2 = frac * fx4 + (1 - frac) * fx3;

    if (cphi < cphil) frac = (cphi - cphil / (nphi1 - 1) * k) / (cphil / (nphi1 - 1));
    else frac = (cphi - cphil - (1.0 - cphil) / nphi2 * (k - nphi1 + 1)) / ((1.0 - cphil) / nphi2);
    double fz = frac * fy2 + (1 - frac) * fy1;

    if ((j == (nR_01 + nR_02 - 1)) || (k == (nphi1 + nphi2 - 1))) fz = 0.0;

    return fz;
}

namespace {
// sin(x) and cos(x) of the same angle, computed together in one shared
// range-reduction pass where the platform's libm exposes that as a single
// entry point, instead of two separate std::sin()/std::cos() calls each
// paying for their own range reduction -- both are needed together at
// every phi1/phi2/phi evaluation below. Falls back to the portable
// two-call form on any other platform: functionally identical either way,
// this only affects speed, never the result.
inline void sincos_local(double x, double& s, double& c) {
#if defined(__APPLE__)
    // Darwin's libm extension, declared in <math.h> in C++ mode.
    __sincos(x, &s, &c);
#elif defined(__GLIBC__)
    // glibc's own sincos() extension.
    ::sincos(x, &s, &c);
#else
    s = std::sin(x);
    c = std::cos(x);
#endif
}
} // namespace

// Tier 1 (once per wall/face): the B1/B2/B3 boundary constants, functions
// of r0 and h alone (thesis eq. 4.64-4.66) -- see the "Per the thesis
// derivation" comment on the header declaration for why these are safe to
// hoist above the per-edge, per-vertex loops that call this.
//
// Reproduces full_integral_new()'s r0==0 -> 0 guard via `valid`: everything
// below r0h_2/r0h_3 divides by r0, so we skip that arithmetic entirely
// rather than compute-then-discard a NaN.
WallIntegralConstants compute_wall_constants(double h, double r0) {
    WallIntegralConstants w;
    w.h = h;
    w.r0 = r0;
    if (r0 == 0.0) return w;
    w.valid = true;

    w.h2 = h * h;
    w.r03 = r0 * r0 * r0;
    w.r0h2 = r0 / h * r0 / h;
    w.r0h3 = w.r0h2 * r0 / h;
    w.r0h_2 = h / r0 * h / r0;
    w.r0h_3 = w.r0h_2 * h / r0;

    if (r0 >= 2.0 * h) {
        w.B3 = w.h2 * h / 4.;
    } else if (r0 > h) {
        w.B3 = w.r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * w.r0h2 + 1. / 30. * w.r0h3 - 1. / 15. * w.r0h_3 + 8. / 5. * w.r0h_2);
        w.B2 = w.r03 / 4. * (-4. / 3. + (r0 / h) - 0.3 * w.r0h2 + 1. / 30. * w.r0h3 - 1. / 15. * w.r0h_3);
    } else {
        w.B3 = w.r03 / 4. * (-2. / 3. + 0.3 * w.r0h2 - 0.1 * w.r0h3 + 7. / 5. * w.r0h_2);
        w.B2 = w.r03 / 4. * (-2. / 3. + 0.3 * w.r0h2 - 0.1 * w.r0h3 - 1. / 5. * w.r0h_2);
        w.B1 = w.r03 / 4. * (-2. / 3. + 0.3 * w.r0h2 - 0.1 * w.r0h3);
    }
    return w;
}

// Tier 2 (once per edge): the D2/D3 region-matching constants -- the
// "phi1 business"/"phi2 business" blocks, functions of r0, R_0 and h only
// (never the caller's per-vertex phi). Reproduces the original's R_0==0 ->
// 0 guard via `valid`.
EdgeIntegralConstants compute_edge_constants(const WallIntegralConstants& w, double R_0) {
    EdgeIntegralConstants e;
    e.R_0 = R_0;
    if (!w.valid || R_0 == 0.0) return e;
    e.valid = true;

    const double h = w.h, r0 = w.r0, h2 = w.h2;
    const double r03 = w.r03, r0h2 = w.r0h2, r0h3 = w.r0h3, r0h_3 = w.r0h_3;
    const double B1 = w.B1, B2 = w.B2, B3 = w.B3;

    double mu, logs, u, sinp, cosp, cosp2, tanp, phi1, phi2;
    double I0, I1, I_1, I_2, I_3, I_4, I_5;
    double D2 = 0.0, D3 = 0.0;

    double a = R_0 / r0;
    double a2 = a * a;
    e.a = a;
    e.a2 = a2;

    double linedist2 = r0 * r0 + R_0 * R_0;

    if (linedist2 <= h2) {
        ////// phi1 business /////
        phi1 = acos(R_0 / sqrt(h * h - r0 * r0));

        sincos_local(phi1, sinp, cosp);
        cosp2 = cosp * cosp;
        mu = cosp / a / sqrt(1. + cosp2 / a2);

        tanp = sinp / cosp;

        I0 = phi1;
        I_2 = phi1 + a2 * tanp;
        I_4 = phi1 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

        u = sinp * sqrt(1 - mu * mu);
        logs = log(1 + u) - log(1 - u);
        I1 = atan(u / a);

        I_1 = a / 2. * logs + I1;
        I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
        I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. * ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

        D2 = -1. / 6. * I_2 + 0.25 * (r0 / h) * I_3 - 0.15 * r0h2 * I_4 + 1. / 30. * r0h3 * I_5 - 1. / 60. * r0h_3 * I1 + (B1 - B2) / r03 * I0;

        ////// phi2 business /////
        phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

        sincos_local(phi2, sinp, cosp);
        cosp2 = cosp * cosp;
        mu = cosp / a / sqrt(1. + cosp2 / a2);

        tanp = sinp / cosp;

        I0 = phi2;
        I_2 = phi2 + a2 * tanp;
        I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

        u = sinp * sqrt(1 - mu * mu);
        logs = log(1 + u) - log(1 - u);
        I1 = atan(u / a);

        I_1 = a / 2. * logs + I1;
        I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
        I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. * ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

        D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 - 1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 + D2;
    } else if (linedist2 <= 4.0 * h2) {
        ////// phi2 business /////
        phi2 = acos(R_0 / sqrt(4.0 * h * h - r0 * r0));

        sincos_local(phi2, sinp, cosp);
        cosp2 = cosp * cosp;
        mu = cosp / a / sqrt(1. + cosp2 / a2);

        tanp = sinp / cosp;

        I0 = phi2;
        I_2 = phi2 + a2 * tanp;
        I_4 = phi2 + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

        u = sinp * sqrt(1 - mu * mu);
        logs = log(1 + u) - log(1 - u);
        I1 = atan(u / a);

        I_1 = a / 2. * logs + I1;
        I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
        I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. * ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

        // D2 stays 0.0 here, same as the original: this branch never sets it.
        D3 = 1. / 3. * I_2 - 0.25 * (r0 / h) * I_3 + 3. / 40. * r0h2 * I_4 - 1. / 120. * r0h3 * I_5 + 4. / 15. * r0h_3 * I1 + (B2 - B3) / r03 * I0 + D2;
    }
    // else linedist2 > 4h2: D2 = D3 = 0.0, same as the original's defaults.

    e.D2 = D2;
    e.D3 = D3;
    return e;
}

// Tier 3 (once per vertex -- twice per edge, at phi1 and phi2): the only
// part of the computation that actually depends on the caller's phi.
// Reproduces the original's phi==0 -> 0 guard, and returns 0 if the wall or
// edge were invalid (r0==0 or R_0==0), matching full_integral_new() exactly.
double vertex_integral(const WallIntegralConstants& w, const EdgeIntegralConstants& e, double phi) {
    if (!w.valid || !e.valid) return 0.0;
    if (phi == 0.0) return 0.0;

    const double pi = M_PI;
    const double h2 = w.h2, r03 = w.r03, r0h2 = w.r0h2, r0h3 = w.r0h3, r0h_3 = w.r0h_3;
    const double B1 = w.B1, B2 = w.B2, B3 = w.B3;
    const double r0 = w.r0, h = w.h;
    const double R_0 = e.R_0, a = e.a, a2 = e.a2, D2 = e.D2, D3 = e.D3;

    double sinp, cosp;
    sincos_local(phi, sinp, cosp);

    double R = R_0 / cosp;
    double r2 = r0 * r0 + R * R;

    double cosp2 = cosp * cosp;
    double mu = cosp / a / sqrt(1. + cosp2 / a2);
    double tanp = sinp / cosp;

    double I0 = phi;
    double I_2 = phi + a2 * tanp;
    double I_4 = phi + 2 * a2 * tanp + 1. / 3. * a2 * a2 * tanp * (2. + 1. / cosp2);

    double u = sinp * sqrt(1 - mu * mu);
    double logs = log(1 + u) - log(1 - u);
    double I1 = atan(u / a);

    double I_1 = a / 2. * logs + I1;
    double I_3 = I_1 + a * (1. + a2) / 4. * (2 * u / (1 - u * u) + logs);
    double I_5 = I_3 + a * (1. + a2) * (1. + a2) / 16. * ((10 * u - 6 * u * u * u) / (1 - u * u) / (1 - u * u) + 3. * logs);

    double full_int;
    if (r2 < h2) {
        full_int = r0h3 / pi * (1. / 6. * I_2 - 3. / 40. * r0h2 * I_4 + 1. / 40. * r0h3 * I_5 + B1 / r03 * I0);
    } else if (r2 < 4.0 * h2) {
        full_int = r0h3 / pi * (0.25 * (4. / 3. * I_2 - (r0 / h) * I_3 + 0.3 * r0h2 * I_4 -
                                         1. / 30. * r0h3 * I_5 + 1. / 15. * r0h_3 * I1) +
                                 B2 / r03 * I0 + D2);
    } else {
        full_int = r0h3 / pi * (-0.25 * r0h_3 * I1 + B3 / r03 * I0 + D3);
    }

    return full_int;
}

// Direct port of the legacy full_integral_new(), the analytic closed-form
// SPH cubic-spline kernel line integral that gridded_integral() above
// approximates via table lookup. Now a thin wrapper over the three tiers
// above -- kept for spot-checking/API compatibility, and always does the
// full recompute (no hoisting), since it has no way to know a caller wants
// to reuse wall/edge constants across several calls. mapping.cpp's
// compute_cell_integral() calls the three tiers directly instead, hoisting
// the wall-level and edge-level tiers above its per-vertex loop.
double full_integral_new(double h, double phi, double r0, double R_0) {
    if (r0 == 0.0) return 0.0;
    if (R_0 == 0.0) return 0.0;
    if (phi == 0.0) return 0.0;

    WallIntegralConstants w = compute_wall_constants(h, r0);
    EdgeIntegralConstants e = compute_edge_constants(w, R_0);
    return vertex_integral(w, e, phi);
}

// Direct extraction of mapping.cpp's compute_cell_integral() per-face loop
// body -- same operations, same order, generalized from voro++'s packed
// f_vert/f_coord indexing to a plain flattened vertex array. See the header
// declaration for the winding requirement and for how mapping.cpp and
// cartesian_grid.cpp each use this.
namespace {

// The per-face / per-edge / per-vertex kernel-integral evaluators that face_integral_impl() below is
// parameterised on, so the face walk (geometry, sign logic, neighbour-discovery distance) exists
// exactly once for both the cubic-spline fast path and the generic-kernel path.
//   set_wall(ar0)   once per face, with the particle's unsigned distance to the face plane
//   set_edge(R_0)   once per edge, with the in-plane distance from the particle's projection to the edge line
//   at(phi)         once per edge vertex
struct CubicFaceEvaluator {
    double h;
    bool tabulated;
    const DensityLookupTable* table;
    WallIntegralConstants wall;
    EdgeIntegralConstants edge;
    double ar0 = 0.0, R_0 = 0.0;

    void set_wall(double ar0_) {
        ar0 = ar0_;
        // Wall-level kernel-integral constants (B1/B2/B3): a function of (h, ar0) alone, computed once per
        // face and reused for every edge. Only meaningful (and only computed) on the analytic path.
        if (!tabulated) wall = compute_wall_constants(h, ar0);
    }
    void set_edge(double R_0_) {
        R_0 = R_0_;
        // Edge-level constants (D2/D3): a function of (h, ar0, R_0) alone, shared between this edge's two
        // vertex evaluations.
        if (!tabulated) edge = compute_edge_constants(wall, R_0);
    }
    double at(double phi) const {
        return tabulated ? gridded_integral(*table, h, phi, ar0, R_0) : vertex_integral(wall, edge, phi);
    }
};

struct GenericFaceEvaluator {
    const KernelSpec* K;
    double H;
    GenericWallConstants wall;
    GenericEdgeConstants edge;

    void set_wall(double ar0) { wall = generic_wall_constants(*K, H, ar0); }
    void set_edge(double R_0) { edge = generic_edge_constants(wall, R_0); }
    double at(double phi) const { return generic_vertex_integral(wall, edge, phi); }
};

// support_radius: the kernel's support (2h for the cubic spline); a face whose closest point to the
// particle is at least this far away cannot lie inside the kernel, which is what out_vert reports.
template <class Evaluator>
double face_integral_impl(double x0, double y0, double z0, double support_radius,
                          const std::vector<double>& vertices, int* out_vert, Evaluator& ev) {
    const int n = static_cast<int>(vertices.size() / 3);
    int out_vert_count = 0;
    double face_Msum = 0.0;

    if (n < 3) {
        if (out_vert) *out_vert = out_vert_count;
        return face_Msum;
    }

    auto vx = [&](int i) { return vertices[3 * i]; };
    auto vy = [&](int i) { return vertices[3 * i + 1]; };
    auto vz = [&](int i) { return vertices[3 * i + 2]; };

    // Face-level state (the plane, via A,B,C,D -> r0/xp,yp,zp; and the first
    // vertex x1,y1,z1) computed once from this face's first three vertices,
    // then reused for every edge of the face -- same as the j==0 branch in
    // mapping.cpp's loop.
    const double x1 = vx(0), y1 = vy(0), z1 = vz(0);
    const double x2s = vx(1), y2s = vy(1), z2s = vz(1);
    const double x3s = vx(2), y3s = vy(2), z3s = vz(2);

    double A = (y2s - y1) * (z3s - z1) - (y3s - y1) * (z2s - z1);
    double B = (z2s - z1) * (x3s - x1) - (z3s - z1) * (x2s - x1);
    double C = (x2s - x1) * (y3s - y1) - (x3s - x1) * (y2s - y1);
    double D = -A * x1 - B * y1 - C * z1;

    double norm = std::sqrt(A * A + B * B + C * C);
    double r0 = (A * x0 + B * y0 + C * z0 + D) / norm;
    double ar0 = std::fabs(r0);

    double xp = x0 - r0 * A / norm;
    double yp = y0 - r0 * B / norm;
    double zp = z0 - r0 * C / norm;

    // s2 (and s1, per-edge below) is a scalar triple product used purely to
    // pin down a winding-consistent SIGN for M -- geometrically,
    // (V1-P0).((V2-P0)x(V3-P0)) is 6x the signed volume of the tetrahedron
    // (P0,V1,V2,V3), whose sign flips iff P0 crosses the plane through
    // V1,V2,V3, so any fixed P0 not on that plane works; only its sign
    // matters, and s1/s2 must share the same P0 for their product's sign to
    // stay meaningful. This used to take P0 as the coordinate-system origin
    // implicitly (dropping the "-P0" terms below), which put s2 at exactly
    // 0 -- an unconditional wrong-sign flip via s1*s2*r0<=0, independent of
    // the true geometry -- whenever the face's plane happened to pass
    // through the global origin exactly (for an axis-aligned box/grid face,
    // whenever one of its 6 defining wall coordinates is exactly 0.0, e.g. a
    // domain boundary deliberately placed at x/y/z=0, which is common:
    // found via the lattice-fixture cross-check, and directly relevant to
    // matching SPLASH's own grid, whose y=0 domain wall triggered it).
    // P0 = the particle's own position (x0,y0,z0) instead: for that to
    // reintroduce the same failure, the particle would have to sit exactly
    // on this face's plane, a measure-zero coincidence for continuous
    // particle positions (unlike 0.0 wall coordinates, which are a common,
    // deliberate domain choice) -- and r0 (the particle's own signed
    // distance to the plane) is already exactly 0 in that case too, making
    // the M flip moot since r0's own factor in s1*s2*r0 already zeroes the
    // comparison's other input.
    const double x1r = x1 - x0, y1r = y1 - y0, z1r = z1 - z0;
    const double x2sr = x2s - x0, y2sr = y2s - y0, z2sr = z2s - z0;
    const double x3sr = x3s - x0, y3sr = y3s - y0, z3sr = z3s - z0;
    double s2 = x1r * (y2sr * z3sr - z2sr * y3sr) + y1r * (z2sr * x3sr - x2sr * z3sr) +
                z1r * (x2sr * y3sr - y2sr * x3sr);

    ev.set_wall(ar0);

    // Neighbor-discovery distance: the true 3D distance from the particle to
    // the closest point of the (bounded, convex) face -- NOT to any single
    // vertex or edge in isolation. Two ingredients, tracked across the edge
    // loop below and combined once after it:
    //
    //  - whether the particle's in-plane projection (xp,yp,zp) lands INSIDE
    //    the polygon: tested via cross(edge, edge-to-point) . face_normal
    //    for every edge -- for a convex polygon walked in a consistent
    //    winding, that sign is the same for every edge iff the point is
    //    inside (or exactly on the boundary). If it's inside, the closest
    //    point on the face is the projection itself, at distance |r0|.
    //  - otherwise, the closest point is on the boundary: the minimum, over
    //    all edges, of the true distance from the particle to that EDGE
    //    SEGMENT (clamped to the segment's endpoints when the perpendicular
    //    foot falls outside it, via the same foot-in-segment test the
    //    integral's sum/difference branch below already needs).
    //
    // Both facts are necessarily face-level, not per-edge -- a single edge
    // being "far" says nothing about the face as a whole (this fixes a bug:
    // the previous per-edge vertex-only distance under this same variable name
    // silently dropped whole neighbor cells whenever a particle's kernel
    // crossed a large face away from every corner, e.g. a small-h particle
    // inside a coarse/lattice Voronoi cell).
    bool inside_polygon = true;
    bool have_sign = false;
    double sign_ref = 0.0;
    double min_edge_dist = std::numeric_limits<double>::infinity();

    for (int j = 0; j < n; ++j) {
        double x2 = vx(j), y2 = vy(j), z2 = vz(j);
        double x3, y3, z3;
        if (j == n - 1) {
            x3 = x1;
            y3 = y1;
            z3 = z1;
        } else {
            x3 = vx(j + 1);
            y3 = vy(j + 1);
            z3 = vz(j + 1);
        }

        double r23 = std::sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2) + (z3 - z2) * (z3 - z2));
        double r12 = std::sqrt((xp - x2) * (xp - x2) + (yp - y2) * (yp - y2) + (zp - z2) * (zp - z2));
        double r13 = std::sqrt((xp - x3) * (xp - x3) + (yp - y3) * (yp - y3) + (zp - z3) * (zp - z3));
        double cosa = ((x3 - x2) * (xp - x2) + (y3 - y2) * (yp - y2) + (z3 - z2) * (zp - z2)) / r12 / r23;

        double R_0;
        if (std::fabs(cosa) < 1.0) {
            R_0 = r12 * std::sqrt(1 - cosa * cosa);
        } else {
            R_0 = 0.0; // |cosa|-1 within tolerance (degenerate/collinear edge)
        }

        // Same P0=(x0,y0,z0) shift as s2 above (see that comment) -- must
        // match s2's reference point for s1*s2's sign to retain its
        // intended geometric meaning.
        const double xpr = xp - x0, ypr = yp - y0, zpr = zp - z0;
        const double x2r = x2 - x0, y2r = y2 - y0, z2r = z2 - z0;
        const double x3r = x3 - x0, y3r = y3 - y0, z3r = z3 - z0;
        double s1 = xpr * (y2r * z3r - z2r * y3r) + ypr * (z2r * x3r - x2r * z3r) + zpr * (x2r * y3r - y2r * x3r);

        double phi1 = (R_0 < r12) ? std::acos(R_0 / r12) : 0.0;
        double phi2 = (R_0 < r13) ? std::acos(R_0 / r13) : 0.0;

        // This edge's contribution to the two face-level neighbor-discovery
        // facts (see the comment above the loop): whether (xp,yp,zp) is on
        // the polygon-interior side of this edge, and the true
        // segment-clamped in-plane distance to this edge specifically.
        double ex = x3 - x2, ey = y3 - y2, ez = z3 - z2;
        double tpx = xp - x2, tpy = yp - y2, tpz = zp - z2;
        double side = (ey * tpz - ez * tpy) * A + (ez * tpx - ex * tpz) * B + (ex * tpy - ey * tpx) * C;
        const double side_eps = 1e-9 * norm * r23 * std::max({std::fabs(tpx), std::fabs(tpy), std::fabs(tpz), 1.0});
        if (std::fabs(side) > side_eps) {
            if (!have_sign) {
                sign_ref = side;
                have_sign = true;
            } else if ((side > 0) != (sign_ref > 0)) {
                inside_polygon = false;
            }
        }

        bool foot_beyond_endpoint = (r12 * std::sin(phi1) >= r23) || (r13 * std::sin(phi2) >= r23);
        double in_plane_dist = foot_beyond_endpoint ? std::min(r12, r13) : R_0;
        min_edge_dist = std::min(min_edge_dist, in_plane_dist);

        double M = (s1 * s2 * r0 <= 0) ? -1.0 : 1.0;

        ev.set_edge(R_0);

        auto integral_at = [&](double phi_) { return ev.at(phi_); };

        if (foot_beyond_endpoint) {
            if (phi1 >= phi2) {
                M = M * (integral_at(phi1) - integral_at(phi2));
            } else {
                M = M * (integral_at(phi2) - integral_at(phi1));
            }
        } else {
            M = M * (integral_at(phi1) + integral_at(phi2));
        }

        face_Msum += M;
    }

    // Combine the two face-level facts gathered above into a single
    // reachability decision (see the comment above the edge loop): if the
    // particle's projection lands inside the polygon, the closest point on
    // the face is the projection itself (distance |r0|); otherwise it's on
    // the boundary, at the smallest of the per-edge segment-clamped
    // distances found while walking the edges.
    double face_dist = inside_polygon ? ar0 : std::sqrt(r0 * r0 + min_edge_dist * min_edge_dist);
    out_vert_count = (face_dist >= support_radius) ? n : 0;

    if (out_vert) *out_vert = out_vert_count;
    return face_Msum;
}

} // namespace

double compute_face_integral(double x0, double y0, double z0, double h,
                              bool use_tabulated_integral, const DensityLookupTable* table,
                              const std::vector<double>& vertices, int* out_vert) {
    CubicFaceEvaluator ev{h, use_tabulated_integral, table, {}, {}};
    return face_integral_impl(x0, y0, z0, 2 * h, vertices, out_vert, ev);
}

double compute_face_integral_generic(double x0, double y0, double z0, const KernelSpec& K, double H,
                                      const std::vector<double>& vertices, int* out_vert) {
    GenericFaceEvaluator ev{&K, H, {}, {}};
    return face_integral_impl(x0, y0, z0, H, vertices, out_vert, ev);
}

} // namespace kiwi
