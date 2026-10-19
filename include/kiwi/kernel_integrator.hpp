// Picks, once per stage, which face-integral implementation a run uses for its configured kernel:
// the original cubic-spline fast path (analytic or tabulated, kernel_integral.hpp) or the generic
// analytic one (kernel_integral_generic.hpp). Stages hold one KernelIntegrator and call face();
// everything kernel-specific (support radius, which implementation) stays in here.
#pragma once

#include <vector>

#include "kiwi/config.hpp"
#include "kiwi/kernel_integral.hpp"

namespace kiwi {

struct KernelIntegrator {
    const KernelSpec* generic = nullptr; // null: cubic-spline fast path
    double support_in_h = 2.0;           // kernel support radius in units of the particle's smoothing length h
    bool tabulated = false;              // cubic fast path only
    const DensityLookupTable* table = nullptr;

    // Kernel mass of the polygonal face `vertices` for a particle at (x0,y0,z0) with smoothing length h;
    // same winding requirement and out_vert meaning as compute_face_integral().
    double face(double x0, double y0, double z0, double h, const std::vector<double>& vertices,
                int* out_vert = nullptr) const {
        if (generic) return compute_face_integral_generic(x0, y0, z0, *generic, support_in_h * h, vertices, out_vert);
        return compute_face_integral(x0, y0, z0, h, tabulated, table, vertices, out_vert);
    }

    // Support radius of a particle with smoothing length h.
    double support(double h) const { return support_in_h * h; }
};

// Builds the integrator for cfg's kernel settings. `tabulated`/`table` are the stage's own
// use_tabulated_integral setting and lookup table (ignored, and never set, for the generic path --
// Config::load_from_file already rejects that combination).
inline KernelIntegrator make_kernel_integrator(const Config& cfg, bool tabulated, const DensityLookupTable* table) {
    KernelIntegrator ki;
    ki.support_in_h = cfg.kernel_support_in_h();
    if (cfg.use_generic_kernel_path()) {
        ki.generic = &cfg.kernel_spec();
    } else {
        ki.tabulated = tabulated;
        ki.table = table;
    }
    return ki;
}

} // namespace kiwi
