// Regression test for the generic-kernel analytic vertex integral (kiwi/kernel_integral_generic.hpp).
//
//  A. cubic spline: generic_full_integral() must reproduce the hand-written full_integral_new()
//     (H = 2h) at random (r0, R_0, phi);
//  B. all six kernels: generic_full_integral() must agree with brute-force Gauss-Legendre quadrature
//     of the kernel taken straight from its definition (the C++ twin of tools/kernel_tables/reference.py);
//  D. whole-cell masses from compute_face_integral_generic() (boxes): the cubic spline matches
//     compute_face_integral(); a box enclosing the support has mass exactly 1; masses of sub-boxes add
//     up to the mass of the box they partition (this exercises the face-sign logic, including a
//     particle outside the box).
//  C. each kernel is normalised (4 pi * int W r^2 dr = 1) and I_P is continuous across region boundaries.
//
// Usage: test_generic_kernel [n_points_per_kernel]
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "kiwi/kernel_integral.hpp"
#include "kiwi/kernel_integral_generic.hpp"
#include "kiwi/kernel_spec.hpp"

using namespace kiwi;

namespace {

// Shape function w(q) straight from the definitions in thesis appendix A.
double w_direct(const std::string& key, double q) {
    auto p = [](double x, int n) { return x > 0 ? std::pow(x, n) : 0.0; };
    if (q >= 1.0) return 0.0;
    if (key == "M4") return p(1 - q, 3) - 4 * p(0.5 - q, 3);
    if (key == "M5") return p(1 - q, 4) - 5 * p(0.6 - q, 4) + 10 * p(0.2 - q, 4);
    if (key == "M6") return p(1 - q, 5) - 6 * p(2.0 / 3 - q, 5) + 15 * p(1.0 / 3 - q, 5);
    if (key == "WC2") return p(1 - q, 4) * (1 + 4 * q);
    if (key == "WC4") return p(1 - q, 6) * (1 + 6 * q + 35.0 / 3 * q * q);
    if (key == "WC6") return p(1 - q, 8) * (1 + 8 * q + 25 * q * q + 32 * q * q * q);
    std::fprintf(stderr, "unknown kernel %s\n", key.c_str());
    std::exit(2);
}

struct GL {
    std::vector<double> x, w;
    explicit GL(int n) {
        for (int i = 1; i <= n; ++i) {
            double z = std::cos(M_PI * (i - 0.25) / (n + 0.5)), dp = 1;
            for (int it = 0; it < 100; ++it) {
                double p0 = 1, p1 = z;
                for (int k = 2; k <= n; ++k) { double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
                dp = n * (z * p1 - p0) / (z * z - 1);
                double dz = p1 / dp;
                z -= dz;
                if (std::fabs(dz) < 1e-15) break;
            }
            x.push_back(z);
            w.push_back(2 / ((1 - z * z) * dp * dp));
        }
    }
    template <class F> double integrate(F f, double a, double b) const {
        if (b <= a) return 0.0;
        double h = 0.5 * (b - a), m = 0.5 * (b + a), s = 0;
        for (size_t i = 0; i < x.size(); ++i) s += w[i] * f(m + h * x[i]);
        return h * s;
    }
};

// m(rho) = int_0^rho W r^(d-1) dr  (H = 1), W = A w.
double enclosed(const KernelSpec& K, double A, int d, double rho) {
    static const GL g(12);
    std::vector<double> cuts{0.0};
    for (int k = 0; k < K.nreg - 1; ++k) if (K.c[k] < rho) cuts.push_back(K.c[k]);
    cuts.push_back(std::min(rho, 1.0));
    double s = 0;
    for (size_t i = 0; i + 1 < cuts.size(); ++i)
        s += g.integrate([&](double r) { return w_direct(K.key, r) * std::pow(r, d - 1); }, cuts[i], cuts[i + 1]);
    return A * s;
}

// Brute-force integral of W over the vertex pyramid (H = 1): int_0^phi dphi' int_0^{Rmax} m(rho) r0 R / rho^3 dR.
double pyramid(const KernelSpec& K, double r0, double R0, double phi) {
    static const GL gphi(48), gR(48);
    const double A = K.A3_pi / M_PI;
    std::vector<double> pc{0.0};
    for (int k = 0; k < K.nreg - 1; ++k) {
        double c = K.c[k];
        if (c * c - r0 * r0 > R0 * R0) pc.push_back(std::acos(R0 / std::sqrt(c * c - r0 * r0)));
    }
    pc.push_back(phi);
    std::sort(pc.begin(), pc.end());
    auto wedge = [&](double ph) {
        double Rmax = R0 / std::cos(ph);
        std::vector<double> rc{0.0};
        for (int k = 0; k < K.nreg - 1; ++k)
            if (K.c[k] * K.c[k] > r0 * r0 && std::sqrt(K.c[k] * K.c[k] - r0 * r0) < Rmax) rc.push_back(std::sqrt(K.c[k] * K.c[k] - r0 * r0));
        for (double g = r0 * 0.5; g < Rmax; g *= 2.0) rc.push_back(g); // resolve the peak of the integrand near R ~ r0
        rc.push_back(Rmax);
        std::sort(rc.begin(), rc.end());
        double s = 0;
        for (size_t i = 0; i + 1 < rc.size(); ++i)
            s += gR.integrate([&](double R) {
                double rho2 = r0 * r0 + R * R;
                return enclosed(K, A, 3, std::sqrt(rho2)) * r0 * R / (rho2 * std::sqrt(rho2));
            }, rc[i], rc[i + 1]);
        return s;
    };
    double s = 0;
    for (size_t i = 0; i + 1 < pc.size(); ++i)
        if (pc[i + 1] <= phi) s += gphi.integrate(wedge, pc[i], pc[i + 1]);
    return s;
}

// Mass fractions are summed over faces and vertices, so absolute error is what matters; the combined
// measure |a - b| / (ABS_TOL / REL_TOL + max(|a|, |b|)) is <= REL_TOL exactly when |a - b| <= ABS_TOL
// for tiny values and <= REL_TOL * value for large ones (up to a factor 2).
constexpr double ABS_TOL = 1e-12, REL_TOL = 1e-9;
double err(double a, double b) { return std::fabs(a - b) / (ABS_TOL / REL_TOL + std::max(std::fabs(a), std::fabs(b))); }

// The 6 faces of an axis-aligned box, wound as src/cartesian_grid.cpp's box_cell_faces() does
// (the winding compute_face_integral*() needs).
void box_faces(const double lo[3], const double hi[3], std::vector<double> faces[6]) {
    const double xm = lo[0], xp = hi[0], ym = lo[1], yp = hi[1], zm = lo[2], zp = hi[2];
    faces[0] = {xp, ym, zp, xp, yp, zp, xp, yp, zm, xp, ym, zm};
    faces[1] = {xm, yp, zm, xm, yp, zp, xm, ym, zp, xm, ym, zm};
    faces[2] = {xp, yp, zm, xp, yp, zp, xm, yp, zp, xm, yp, zm};
    faces[3] = {xm, ym, zp, xp, ym, zp, xp, ym, zm, xm, ym, zm};
    faces[4] = {xm, yp, zp, xp, yp, zp, xp, ym, zp, xm, ym, zp};
    faces[5] = {xp, ym, zm, xp, yp, zm, xm, yp, zm, xm, ym, zm};
}

double box_mass_generic(const KernelSpec& K, double H, const double p[3], const double lo[3], const double hi[3],
                        int* out_vert_sum = nullptr) {
    std::vector<double> faces[6];
    box_faces(lo, hi, faces);
    double m = 0;
    int ovs = 0;
    for (auto& f : faces) {
        int ov = 0;
        m += compute_face_integral_generic(p[0], p[1], p[2], K, H, f, &ov);
        ovs = ovs * 2 + (ov > 0);
    }
    if (out_vert_sum) *out_vert_sum = ovs;
    return m;
}

double box_mass_cubic(double h, const double p[3], const double lo[3], const double hi[3], int* out_vert_sum = nullptr) {
    std::vector<double> faces[6];
    box_faces(lo, hi, faces);
    double m = 0;
    int ovs = 0;
    for (auto& f : faces) {
        int ov = 0;
        m += compute_face_integral(p[0], p[1], p[2], h, false, nullptr, f, &ov);
        ovs = ovs * 2 + (ov > 0);
    }
    if (out_vert_sum) *out_vert_sum = ovs;
    return m;
}

} // namespace

int main(int argc, char** argv) {
    const int npts = argc > 1 ? std::atoi(argv[1]) : 40;
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    int failures = 0;

    // ---- C1. normalisation
    for (const KernelSpec* K : all_kernels()) {
        double norm = 4 * M_PI * enclosed(*K, K->A3_pi / M_PI, 3, 1.0);
        double e = std::fabs(norm - 1.0);
        std::printf("%-4s normalisation: 4pi int W r^2 dr - 1 = %9.2e\n", K->key, norm - 1.0);
        if (e > 1e-12) { std::printf("  FAIL\n"); ++failures; }
    }

    // ---- A. cubic spline vs the hand-written implementation
    {
        const KernelSpec& K = get_kernel("cubic_spline");
        double worst = 0;
        for (int n = 0; n < 20 * npts; ++n) {
            double H = 0.5 + 6.0 * U(rng);
            double r0 = H * (0.01 + 1.6 * U(rng)), R0 = H * (0.01 + 2.0 * U(rng));
            double phi = 1e-3 + (M_PI / 2 - 2e-3) * U(rng);
            double ref = full_integral_new(H / 2, phi, r0, R0);
            double gen = generic_full_integral(K, H, phi, r0, R0);
            worst = std::max(worst, err(gen, ref));
        }
        std::printf("A. cubic spline, generic vs full_integral_new: max scaled difference %9.2e\n", worst);
        if (worst > REL_TOL) { std::printf("  FAIL\n"); ++failures; }
    }

    // ---- B. all kernels vs brute-force quadrature
    for (const KernelSpec* K : all_kernels()) {
        double worst = 0;
        for (int n = 0; n < npts; ++n) {
            double r0 = 0.02 + 1.4 * U(rng), R0 = 0.02 + 1.6 * U(rng), phi = 0.05 + 1.5 * U(rng);
            double ref = pyramid(*K, r0, R0, phi);
            double gen = generic_full_integral(*K, 1.0, phi, r0, R0);
            worst = std::max(worst, err(gen, ref));
        }
        std::printf("B. %-4s generic vs brute-force quadrature (%d pts): max scaled difference %9.2e\n", K->key, npts, worst);
        if (worst > REL_TOL) { std::printf("  FAIL\n"); ++failures; }
    }

    // ---- B2. stress geometry: tiny r0 (large R_0 / r0), r0 beyond the support, vertices far away
    {
        std::uniform_real_distribution<double> L(std::log(1e-3), std::log(1.0));
        for (const KernelSpec* K : all_kernels()) {
            double worst = 0, wr0 = 0, wR0 = 0, wphi = 0;
            for (int n = 0; n < npts; ++n) {
                double r0 = std::exp(L(rng)) * (U(rng) < 0.15 ? 2.0 : 1.0), R0 = std::exp(L(rng)) * 2.0;
                double phi = 0.05 + 1.5 * U(rng);
                double e = err(generic_full_integral(*K, 1.0, phi, r0, R0), pyramid(*K, r0, R0, phi));
                if (e > worst) { worst = e; wr0 = r0; wR0 = R0; wphi = phi; }
            }
            std::printf("B2. %-4s stress geometry (%d pts): max scaled difference %9.2e (r0=%.3g R0=%.3g phi=%.3g)\n",
                        K->key, npts, worst, wr0, wR0, wphi);
            if (worst > REL_TOL) { std::printf("  FAIL\n"); ++failures; }
        }
    }

    // ---- C2. continuity of I_P across region boundaries
    for (const KernelSpec* K : all_kernels()) {
        double worst = 0;
        for (int n = 0; n < npts; ++n) {
            double r0 = 0.05 + 0.9 * U(rng), R0 = 0.05 + 0.5 * U(rng);
            GenericWallConstants w = generic_wall_constants(*K, 1.0, r0);
            GenericEdgeConstants e = generic_edge_constants(w, R0);
            for (int s = 0; s + 1 < e.nseg; ++s) {
                double pb = e.phi_end[s], d = 1e-9;
                double lo = generic_vertex_integral(w, e, pb - d), hi = generic_vertex_integral(w, e, pb + d);
                worst = std::max(worst, std::fabs(hi - lo));
            }
        }
        std::printf("C. %-4s jump of I_P across region boundaries: max %9.2e\n", K->key, worst);
        if (worst > 1e-7) { std::printf("  FAIL\n"); ++failures; }
    }

    // ---- D. cell masses from boxes
    {
        double worst_cubic = 0;
        int ov_mismatch = 0;
        for (int n = 0; n < 10 * npts; ++n) {
            double H = 0.5 + 3.0 * U(rng), p[3] = {0, 0, 0}, lo[3], hi[3];
            for (int d = 0; d < 3; ++d) {
                double c = H * (U(rng) - 0.5) * 3.0, half = H * (0.05 + 1.2 * U(rng));
                lo[d] = c - half; hi[d] = c + half;
                p[d] = H * (U(rng) - 0.5) * 0.8;
            }
            int ovg = 0, ovc = 0;
            double mg = box_mass_generic(get_kernel("cubic_spline"), H, p, lo, hi, &ovg);
            double mc = box_mass_cubic(H / 2, p, lo, hi, &ovc);
            worst_cubic = std::max(worst_cubic, err(mg, mc));
            if (ovg != ovc) ++ov_mismatch;
        }
        std::printf("D1. cubic box masses, generic vs compute_face_integral: max scaled difference %9.2e, out_vert mismatches %d\n",
                    worst_cubic, ov_mismatch);
        if (worst_cubic > REL_TOL || ov_mismatch) { std::printf("  FAIL\n"); ++failures; }

        for (const KernelSpec* K : all_kernels()) {
            double worst_one = 0, worst_part = 0, worst_grid = 0;
            for (int n = 0; n < 4 * npts; ++n) {
                double H = 0.5 + 3.0 * U(rng), p[3], lo[3], hi[3];
                for (int d = 0; d < 3; ++d) p[d] = H * (U(rng) - 0.5) * 2.0;
                // a box enclosing the whole support has mass 1
                for (int d = 0; d < 3; ++d) { lo[d] = p[d] - H * (1.0 + 2.0 * U(rng)); hi[d] = p[d] + H * (1.0 + 2.0 * U(rng)); }
                worst_one = std::max(worst_one, std::fabs(box_mass_generic(*K, H, p, lo, hi) - 1.0));
                // partition of a random box (particle anywhere, often outside) along a random axis
                double blo[3], bhi[3];
                for (int d = 0; d < 3; ++d) { double c = H * (U(rng) - 0.5) * 2.0, half = H * (0.1 + 0.9 * U(rng)); blo[d] = c - half; bhi[d] = c + half; }
                int ax = n % 3;
                double cut = blo[ax] + (bhi[ax] - blo[ax]) * (0.2 + 0.6 * U(rng));
                double lo2[3] = {blo[0], blo[1], blo[2]}, hi1[3] = {bhi[0], bhi[1], bhi[2]};
                hi1[ax] = cut; lo2[ax] = cut;
                double whole = box_mass_generic(*K, H, p, blo, bhi);
                double parts = box_mass_generic(*K, H, p, blo, hi1) + box_mass_generic(*K, H, p, lo2, bhi);
                worst_part = std::max(worst_part, std::fabs(whole - parts));
                // a 2x2x2 grid covering the whole support sums to 1
                double mid[3], tot = 0;
                for (int d = 0; d < 3; ++d) { lo[d] = p[d] - 1.3 * H; hi[d] = p[d] + 1.3 * H; mid[d] = lo[d] + (hi[d] - lo[d]) * (0.2 + 0.6 * U(rng)); }
                for (int c = 0; c < 8; ++c) {
                    double cl[3], ch[3];
                    for (int d = 0; d < 3; ++d) { bool up = (c >> d) & 1; cl[d] = up ? mid[d] : lo[d]; ch[d] = up ? hi[d] : mid[d]; }
                    tot += box_mass_generic(*K, H, p, cl, ch);
                }
                worst_grid = std::max(worst_grid, std::fabs(tot - 1.0));
            }
            std::printf("D2. %-4s box mass: |enclosing box - 1| %8.2e, |partition - whole| %8.2e, |2x2x2 grid - 1| %8.2e\n",
                        K->key, worst_one, worst_part, worst_grid);
            if (worst_one > 1e-10 || worst_part > 1e-10 || worst_grid > 1e-10) { std::printf("  FAIL\n"); ++failures; }
        }
    }

    std::printf(failures ? "FAILED (%d)\n" : "all checks passed\n", failures);
    return failures ? 1 : 0;
}
