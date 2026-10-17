// Regression test for the generic-kernel analytic vertex integral (kiwi/kernel_integral_generic.hpp).
//
//  A. cubic spline: generic_full_integral() must reproduce the hand-written full_integral_new()
//     (H = 2h) at random (r0, R_0, phi);
//  B. all six kernels: generic_full_integral() must agree with brute-force Gauss-Legendre quadrature
//     of the kernel taken straight from its definition (the C++ twin of tools/kernel_tables/reference.py);
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

    std::printf(failures ? "FAILED (%d)\n" : "all checks passed\n", failures);
    return failures ? 1 : 0;
}
