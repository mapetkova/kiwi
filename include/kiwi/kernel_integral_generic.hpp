// Generic-kernel counterpart of kernel_integral.hpp's analytic vertex integral.
//
// Same three-tier structure as WallIntegralConstants / EdgeIntegralConstants / vertex_integral()
// (wall = once per face, edge = once per edge, vertex = twice per edge), but driven by the
// coefficient tables in KernelSpec instead of hand-written cubic-spline constants, so it works for
// every kernel in kernel_spec.hpp.  The cubic-spline fast path in kernel_integral.hpp is a separate
// implementation and serves as the cross-check for this one (tools/test_generic_kernel.cpp).
//
// Maths (thesis ch. 4 and appendix A; derivation in tools/kernel_tables/derive.py): the kernel mass
// of the "vertex pyramid" spanned by the particle, the foot of the perpendicular on the wall plane
// (distance r0), the foot on the wall edge (distance R_0 from it) and the edge vertex at azimuth phi is
//
//   I_P(phi) = integral_0^phi dI_P/dphi' dphi',
//   dI_P/dphi = A_3 * [ sum_{i != 0} a_i^(k) mu^i q0^(1-i) + B^(j,k)(q0) ],   mu = r0 / r,
//
// where k is the region of the vertex radius r = sqrt(r0^2 + R_0^2 / cos^2 phi) and j that of r0.
// Integrating term by term needs I_i(phi) = integral_0^phi mu^i dphi' for i = 1, 0, -1, ..., -13,
// which have closed forms (see generic_basis_integrals() in the .cpp).  Where the vertex radius
// crosses a breakpoint c_k H the region changes; the edge tier records those angles and the offsets
// that keep I_P continuous across them.
#pragma once

#include "kiwi/kernel_spec.hpp"

namespace kiwi {

struct GenericWallConstants {
    bool valid = false; // r0 != 0
    const KernelSpec* K = nullptr;
    double H = 0.0, r0 = 0.0, q0 = 0.0;
    int j = 0;          // region of q0 (1-based)
    // Per vertex region k (index k-1): the table coefficients with A_3 and the powers of q0 folded in.
    double coef[KernelSpec::MAX_REGIONS][KernelSpec::NI]; // coef[k][1 - i] multiplies I_i
    double bphi[KernelSpec::MAX_REGIONS];                  // multiplies I_0 = phi
};

struct GenericEdgeConstants {
    bool valid = false; // wall.valid && R_0 != 0
    double R_0 = 0.0, a = 0.0, a2 = 0.0;
    int nseg = 0;       // number of regions the edge passes through
    int region[KernelSpec::MAX_REGIONS];     // region k of each segment
    double phi_end[KernelSpec::MAX_REGIONS]; // segment s covers (phi_end[s-1], phi_end[s]]; the last is pi/2
    double offset[KernelSpec::MAX_REGIONS];  // I_P(phi) = offset[s] + G_k(phi) inside segment s
};

// Once per wall: r0 is the (unsigned) perpendicular distance, H the kernel support radius.
GenericWallConstants generic_wall_constants(const KernelSpec& K, double H, double r0);

// Once per edge.
GenericEdgeConstants generic_edge_constants(const GenericWallConstants& wall, double R_0);

// Once per vertex; phi in [0, pi/2).
double generic_vertex_integral(const GenericWallConstants& wall, const GenericEdgeConstants& edge, double phi);

// All three tiers in one call, mirroring full_integral_new(h, phi, r0, R_0) (there, H = 2h).
double generic_full_integral(const KernelSpec& K, double H, double phi, double r0, double R_0);

} // namespace kiwi
