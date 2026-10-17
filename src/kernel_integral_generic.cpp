#include "kiwi/kernel_integral_generic.hpp"

#include <cmath>

namespace kiwi {

namespace {

constexpr int NI = KernelSpec::NI;

// I_i(phi) = integral_0^phi mu^i dphi' for i = 1, 0, -1, ..., -(NI-2), stored at out[1 - i], where
// mu = cos(phi) / sqrt(a^2 + cos^2(phi))   (a = R_0 / r0).  Only i >= imin are filled.
//
//  even i = -2k :  mu^-2 = 1 + a^2 sec^2,  so  I_-2k = sum_l C(k,l) a^(2l) J_2l,
//                  J_m = integral_0^phi sec^m  (recurrence J_m = sec^(m-2) tan / (m-1) + (m-2)/(m-1) J_(m-2)).
//  odd  i        :  with u = a sin(phi) / sqrt(a^2 + cos^2), I_1 = atan(u / a) and
//                  I_(-2m-1) = I_(-2m+1) + a (1 + a^2)^m K_(m+1),  K_n = integral du / (1 - u^2)^n
//                  (K_1 = atanh u,  K_(n+1) = u / (2n (1-u^2)^n) + (2n-1)/(2n) K_n);  I_-1 = a K_1 + I_1.
void basis_integrals(double a, double a2, double phi, int imin, double* out) {
    const double sinp = std::sin(phi), cosp = std::cos(phi);
    const double den = a2 + cosp * cosp;
    const double u = a * sinp / std::sqrt(den);
    const double one_minus_u2 = cosp * cosp * (1.0 + a2) / den; // 1 - u^2 without cancellation

    out[0] = std::atan(u / a);  // I_1
    out[1] = phi;               // I_0

    // odd negative powers
    const int mmax = (-imin - 1) / 2;       // I_(-2m-1) for m = 0..mmax
    double Ioddprev = out[0];
    double K = std::atanh(u);
    double pow1a2 = 1.0;
    for (int m = 0; m <= mmax && -2 * m - 1 >= imin; ++m) {
        if (m > 0) {
            const int n = m;                                  // K_(m+1) from K_m
            K = u / (2.0 * n * std::pow(one_minus_u2, n)) + (2.0 * n - 1.0) / (2.0 * n) * K;
            pow1a2 *= (1.0 + a2);
        }
        Ioddprev += a * pow1a2 * K;
        out[1 + (2 * m + 1)] = Ioddprev;                      // index 1 - i with i = -(2m+1)
    }

    // even negative powers
    const int kmax = -imin / 2;
    if (kmax >= 1) {
        double J[NI];                                          // J[2l] for l = 0..kmax
        const double t = sinp / cosp, s = 1.0 / cosp;
        J[0] = phi;
        J[2] = t;
        double spow = 1.0;                                     // sec^(m-2)
        for (int l = 2; l <= kmax; ++l) {
            const int m = 2 * l;
            spow *= s * s;                                     // sec^(m-2) for m = 2l  (starts at sec^2 for l=2)
            J[m] = spow * t / (m - 1.0) + (m - 2.0) / (m - 1.0) * J[m - 2];
        }
        for (int k = 1; k <= kmax; ++k) {
            double sum = J[0], apow = 1.0, b = 1.0;            // l = 0 term: C(k,0) a^0 J_0
            for (int l = 1; l <= k; ++l) {
                b = b * (k - l + 1) / l;                       // C(k,l)
                apow *= a2;
                sum += b * apow * J[2 * l];
            }
            out[1 + 2 * k] = sum;                              // i = -2k
        }
    }
}

int region_of(const KernelSpec& K, double q) {
    for (int k = 0; k < K.nreg - 1; ++k)
        if (q <= K.c[k]) return k + 1;
    return K.nreg;
}

// Most negative power present in region k's coefficients.
int imin_of(const GenericWallConstants& w, int k) {
    for (int idx = NI - 1; idx >= 0; --idx)
        if (w.coef[k - 1][idx] != 0.0) return 1 - idx;
    return 1;
}

// G_k(phi) = sum_i coef_i I_i(phi) + bphi * phi  (the antiderivative of dI_P/dphi in region k, zero at phi = 0).
double G(const GenericWallConstants& w, const GenericEdgeConstants& e, int k, double phi) {
    const int imin = imin_of(w, k);
    double I[NI];
    basis_integrals(e.a, e.a2, phi, imin, I);
    double sum = w.bphi[k - 1] * phi;
    for (int idx = 0; idx < NI; ++idx) {
        const double c = w.coef[k - 1][idx];
        if (c != 0.0) sum += c * I[idx];
    }
    return sum;
}

} // namespace

GenericWallConstants generic_wall_constants(const KernelSpec& K, double H, double r0) {
    GenericWallConstants w;
    w.K = &K;
    w.H = H;
    w.r0 = r0;
    if (r0 == 0.0) return w;
    w.valid = true;
    w.q0 = r0 / H;
    w.j = region_of(K, w.q0);
    const double A3 = K.A3_pi / M_PI;
    for (int k = 1; k <= K.nreg; ++k) {
        for (int idx = 0; idx < NI; ++idx) {
            const double a = K.a[k - 1][idx];
            w.coef[k - 1][idx] = (a == 0.0) ? 0.0 : A3 * a * std::pow(w.q0, idx); // q0^(1-i) with i = 1 - idx
        }
        double b = 0.0;
        if (k >= w.j)
            for (int p = KernelSpec::NB - 1; p >= 0; --p) b = b * w.q0 + K.B[w.j - 1][k - 1][p];
        w.bphi[k - 1] = A3 * b;
    }
    return w;
}

GenericEdgeConstants generic_edge_constants(const GenericWallConstants& w, double R_0) {
    GenericEdgeConstants e;
    e.R_0 = R_0;
    if (!w.valid || R_0 == 0.0) return e;
    e.valid = true;
    e.a = R_0 / w.r0;
    e.a2 = e.a * e.a;
    const KernelSpec& K = *w.K;
    const double qR = R_0 / w.H;
    const double r_start = std::sqrt(w.q0 * w.q0 + qR * qR); // vertex radius at phi = 0, in units of H
    int k = region_of(K, r_start);
    e.nseg = 0;
    for (; k <= K.nreg; ++k) {
        e.region[e.nseg] = k;
        if (k == K.nreg) {
            e.phi_end[e.nseg] = 0.5 * M_PI;
        } else {
            const double c = K.c[k - 1];
            e.phi_end[e.nseg] = std::acos(std::fmin(1.0, qR / std::sqrt(c * c - w.q0 * w.q0)));
        }
        ++e.nseg;
    }
    // Offsets: I_P(phi) = offset[s] + G_k(phi); continuity at each segment boundary.
    e.offset[0] = 0.0;
    for (int s = 0; s + 1 < e.nseg; ++s) {
        const double pb = e.phi_end[s];
        e.offset[s + 1] = e.offset[s] + G(w, e, e.region[s], pb) - G(w, e, e.region[s + 1], pb);
    }
    return e;
}

double generic_vertex_integral(const GenericWallConstants& w, const GenericEdgeConstants& e, double phi) {
    if (!w.valid || !e.valid) return 0.0;
    if (phi == 0.0) return 0.0;
    int s = 0;
    while (s + 1 < e.nseg && phi > e.phi_end[s]) ++s;
    return e.offset[s] + G(w, e, e.region[s], phi);
}

double generic_full_integral(const KernelSpec& K, double H, double phi, double r0, double R_0) {
    if (r0 == 0.0 || R_0 == 0.0 || phi == 0.0) return 0.0;
    const GenericWallConstants w = generic_wall_constants(K, H, r0);
    const GenericEdgeConstants e = generic_edge_constants(w, R_0);
    return generic_vertex_integral(w, e, phi);
}

} // namespace kiwi
