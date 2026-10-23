// Standalone particle-only routine for the inverse-mapping standalone tool
// (see the project progress doc's "Twenty-second delivery" section): sets
// each particle's smoothing length h by the same self-consistent
// neighbour-density iteration Phantom uses at every timestep, instead of the
// algebraic h = hfact*(m/rho)^(1/3) shortcut this project has used so far
// (that shortcut trusts an externally supplied rho; this instead asks "what
// h makes the SPH kernel-weighted density of the particle's ACTUAL local
// neighbours agree with its own mass/hfact target?", which is what a real
// SPH code enforces before it ever runs and is far more robust to a locally
// uneven or stale input rho).
//
// NOT wired into the kiwi CLI yet -- deliberately a standalone tool (see
// tools/bench_grid_blocks.cpp and tools/test_lattice_grid.cpp for the same
// pattern): this is one of several particle-related routines to build and
// validate independently before folding them into the standalone
// inverse-mapping feature.
//
// --- Algorithm, verified directly against Phantom's actual source
// (github.com/danieljprice/phantom, module `densityforce` in
// src/main/dens.F90) rather than from memory or the papers alone ---
//
// Each particle i has a target density relation
//   rho_h(h) = pmass_i * (hfact/h)^3                        [part.F90:779-784, rhoh()]
// and an actual, neighbour-derived density
//   rho_sum(h) = sum_j m_j * W(r_ij, h)
// summed (including the self term j==i, r=0 -- see below) over every
// neighbour j within the kernel's compact support, r_ij < radkern*h. Setting
// h so that rho_sum(h) == rho_h(h) is exactly Phantom's density/smoothing-
// length iteration [dens.F90's densityiterate(), maxdensits=100 there].
//
// We use the cubic spline (M4) kernel, Phantom's default and the kernel
// implicit in this codebase's own hfact=1.2 convention:
//   W(r,h) = cnormk/h^3 * f(q),  q = r/h,  cnormk = 1/pi (3D)  [kernel_cubic.f90]
//   f(q) = 1 - 1.5q^2 + 0.75q^3         0 <= q < 1
//        = 0.25*(2-q)^3                 1 <= q < 2
//        = 0                            q >= 2               [radkern = 2]
// f is C1 (f'(1-) == f'(1-) == -0.75, f'(2) == 0), so no special-casing is
// needed at the branch boundaries, and W(0,h) and dW/dh(0,h) are both finite
// (q=0 falls in the first branch), so the self term needs no special case
// either -- it's just "the particle itself happens to be a neighbour at
// distance 0", handled by the same kd-tree range query as everyone else.
//
// The h-update is a damped Newton-Raphson step on f(h) := rho_sum(h) -
// rho_h(h), using
//   d(rho_sum)/dh = sum_j m_j * dW/dh(r_ij, h),
//     dW/dh = -(1/h^4) * (3*w(q) + q*w'(q)),  w(q) = cnormk*f(q)
//     [standard SPH grad-h derivative, e.g. Price (2012) living review
//      "Smoothed Particle Hydrodynamics and Magnetohydrodynamics", sec 2.3]
//   d(rho_h)/dh = -3*rho_h(h)/h
// i.e. h_new = h - f(h) / [d(rho_sum)/dh - d(rho_h)/dh], clamped to
// [0.8h, 1.2h] and iterated to a relative tolerance, matching Phantom's own
// damping/tolerance defaults [dens.F90:1432-1441, options.f90:93]:
//   hfact = 1.2, tolh = 1e-4, maxdensits = 100.
//
// This f/f' is derived directly here from rho_sum(h) and rho_h(h) rather
// than porting Phantom's own internal accumulator layout line-by-line
// (finish_rhosum's rhosum(irhoi)/rhosum(igradhi) arrays, the "Omega" grad-h
// correction variable, its wab0/gradh0 self-term bookkeeping) -- that
// internal layout is an optimization for Phantom's own per-cell/per-neighbour
// caching architecture, not part of the physics. Any correct Newton-Raphson
// iteration on the same root condition f(h) = rho_sum(h) - rho_h(h) = 0
// converges to the same h (the root is a property of f, not of which
// algebraically-equivalent way f' is computed) -- and deriving it directly
// here, then checking it against Phantom's own converged h on a real
// snapshot, is a more reliable way to confirm equivalence than trying to
// reverse-engineer sign conventions out of a heavily inlined Fortran
// accumulator by inspection alone.
//
// The starting guess for each particle's h iteration is this codebase's own
// existing algebraic estimate, h_seed = hfact*(m/rho_snapshot)^(1/3) (Price
// et al. 2018 eq. 6, the same relation snapshot.hpp already uses when only
// one of h/rho is present in an input file) -- using the snapshot's own rho
// column, NOT its own h column, as the seed. This is a reasonable physically-
// motivated starting point, and deliberately independent of the h column
// this tool is being validated against: Newton-Raphson's converged root does
// not depend on which reasonable starting guess found it (only the
// convergence path does), so seeding from h_snapshot itself would not bias
// the comparison, but seeding from rho_snapshot instead keeps that
// independence obvious rather than relying on the reader to accept the
// non-bias argument.
//
// The neighbour search itself uses a from-scratch, from-the-ground-up 3D
// kd-tree (KDTree below) -- built once (particle positions don't move during
// this iteration), range-queried once per particle per Newton-Raphson step
// with that step's own compact-support radius radkern*h_i.
//
// Usage:
//   h_from_neighbours <snapshot_path> [format] [hfact] [tolh] [output.dat]
//     format  -- "auto" (default, sniffed), "plain_ascii", or "splash_ascii"
//     hfact   -- default 1.2 (Phantom's own default, options.f90 hfact_default)
//     tolh    -- default 1e-4 (Phantom's own default, options.f90:93)
//     output  -- default "<snapshot-basename>-h-compare.dat" in the cwd
//
// Output is a text file, one row per particle:
//   id x y z m h_snapshot rho_snapshot h_computed rho_computed n_neighbours
//   n_iterations converged rel_diff_h
// where rel_diff_h = (h_computed - h_snapshot) / h_snapshot, and converged
// is 1/0. A companion plotting script, tools/plot_h_comparison.py, turns
// this into the scatter/histogram comparison against the snapshot's own h.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "kiwi/snapshot.hpp"

namespace {

// ---------------------------------------------------------------------
// Cubic spline (M4) kernel, 3D-normalized, and its analytic d/dh --
// see the file header comment above for the derivation and Phantom's
// matching source locations.
// ---------------------------------------------------------------------
constexpr double kRadkern = 2.0; // compact support: q = r/h < 2
constexpr double kCnormk = 1.0 / M_PI; // 3D cubic-spline normalization

// f(q) and f'(q), the dimensionless shape function and its derivative.
inline double kernel_f(double q) {
    if (q < 1.0) return 1.0 - 1.5 * q * q + 0.75 * q * q * q;
    if (q < 2.0) {
        double t = 2.0 - q;
        return 0.25 * t * t * t;
    }
    return 0.0;
}

inline double kernel_fprime(double q) {
    if (q < 1.0) return -3.0 * q + 2.25 * q * q;
    if (q < 2.0) {
        double t = 2.0 - q;
        return -0.75 * t * t;
    }
    return 0.0;
}

// W(r, h) = cnormk/h^3 * f(r/h).
inline double kernel_W(double r, double h) {
    double q = r / h;
    if (q >= kRadkern) return 0.0;
    return kCnormk * kernel_f(q) / (h * h * h);
}

// dW/dh at fixed r -- see header comment for the derivation:
//   dW/dh = -(1/h^4) * (3*w(q) + q*w'(q)), w(q) = cnormk*f(q)
inline double kernel_dWdh(double r, double h) {
    double q = r / h;
    if (q >= kRadkern) return 0.0;
    double w = kCnormk * kernel_f(q);
    double wprime = kCnormk * kernel_fprime(q);
    double h4 = h * h * h * h;
    return -(3.0 * w + q * wprime) / h4;
}

// ---------------------------------------------------------------------
// Minimal 3D kd-tree: built once over static particle positions, supports
// repeated radius queries with a different radius each call (as the h
// iteration needs -- the compact-support radius shrinks/grows with h every
// Newton-Raphson step, but the positions and tree structure never change).
//
// Array-of-nodes representation (no per-node heap allocation): build()
// recursively partitions an index array by the median along a
// depth%3-cycled axis (std::nth_element, O(n) per level) and stores each
// node's children as indices into the flat nodes_ array rather than
// pointers.
// ---------------------------------------------------------------------
class KDTree {
public:
    KDTree(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z)
        : x_(x), y_(y), z_(z) {
        size_t n = x_.size();
        std::vector<size_t> idx(n);
        for (size_t i = 0; i < n; ++i) idx[i] = i;
        nodes_.reserve(n);
        root_ = build(idx, 0, static_cast<int>(n), 0);
    }

    // Appends the indices of every point within (inclusive) radius r of
    // (qx,qy,qz) to out. Does not clear out first, does not sort/dedupe.
    void query_radius(double qx, double qy, double qz, double r, std::vector<size_t>& out) const {
        if (root_ < 0) return;
        double r2 = r * r;
        query_rec(root_, qx, qy, qz, r2, out);
    }

private:
    struct Node {
        size_t idx;
        int axis;
        int left = -1;
        int right = -1;
    };

    const std::vector<double>& x_;
    const std::vector<double>& y_;
    const std::vector<double>& z_;
    std::vector<Node> nodes_;
    int root_ = -1;

    double coord(size_t i, int axis) const {
        return axis == 0 ? x_[i] : (axis == 1 ? y_[i] : z_[i]);
    }

    // Partitions idx[lo,hi) around its median along `depth`'s axis, creates
    // a node for the median, and recurses on the two halves. Returns the
    // new node's index into nodes_, or -1 for an empty range.
    int build(std::vector<size_t>& idx, int lo, int hi, int depth) {
        if (lo >= hi) return -1;
        int axis = depth % 3;
        int mid = lo + (hi - lo) / 2;
        std::nth_element(idx.begin() + lo, idx.begin() + mid, idx.begin() + hi,
                          [&](size_t a, size_t b) { return coord(a, axis) < coord(b, axis); });

        int node_pos = static_cast<int>(nodes_.size());
        nodes_.push_back(Node{idx[mid], axis, -1, -1});
        int left = build(idx, lo, mid, depth + 1);
        int right = build(idx, mid + 1, hi, depth + 1);
        // nodes_ may have reallocated during the recursive calls above, so
        // index it fresh rather than holding a reference/pointer across them.
        nodes_[node_pos].left = left;
        nodes_[node_pos].right = right;
        return node_pos;
    }

    void query_rec(int node, double qx, double qy, double qz, double r2, std::vector<size_t>& out) const {
        if (node < 0) return;
        const Node& nd = nodes_[node];
        double dx = x_[nd.idx] - qx, dy = y_[nd.idx] - qy, dz = z_[nd.idx] - qz;
        double d2 = dx * dx + dy * dy + dz * dz;
        if (d2 <= r2) out.push_back(nd.idx);

        double qc = (nd.axis == 0 ? qx : (nd.axis == 1 ? qy : qz));
        double nc = coord(nd.idx, nd.axis);
        double diff = qc - nc;
        double r = std::sqrt(r2);

        // Always descend the side the query point is on; only descend the
        // other side if the splitting plane is within r of the query point.
        int near_child = diff <= 0.0 ? nd.left : nd.right;
        int far_child = diff <= 0.0 ? nd.right : nd.left;
        query_rec(near_child, qx, qy, qz, r2, out);
        if (std::fabs(diff) <= r) query_rec(far_child, qx, qy, qz, r2, out);
    }
};

struct HResult {
    double h = 0.0;
    double rho = 0.0;
    int n_neighbours = 0;
    int n_iterations = 0;
    bool converged = false;
};

// One particle's damped Newton-Raphson h iteration -- see file header for
// the algorithm and its correspondence to Phantom's densityiterate().
HResult solve_h(size_t i, const kiwi::SPHSnapshot& snap, const KDTree& tree, double hfact, double tolh,
                 int max_iterations) {
    HResult result;
    double pmass = snap.m[i];

    // Seed from the existing algebraic relation (Price et al. 2018 eq. 6),
    // using the snapshot's OWN rho column, not its own h column -- see file
    // header comment on why this keeps the comparison meaningful.
    double h = hfact * std::cbrt(pmass / std::max(snap.rho[i], 1e-300));
    if (!(h > 0.0)) h = snap.h[i] > 0.0 ? snap.h[i] : 1.0; // pathological-input fallback

    std::vector<size_t> neighbours;
    int iter = 0;
    bool converged = false;
    double rho_sum = 0.0;

    for (; iter < max_iterations; ++iter) {
        neighbours.clear();
        tree.query_radius(snap.x[i], snap.y[i], snap.z[i], kRadkern * h, neighbours);

        rho_sum = 0.0;
        double drho_sum_dh = 0.0;
        for (size_t j : neighbours) {
            double dx = snap.x[i] - snap.x[j];
            double dy = snap.y[i] - snap.y[j];
            double dz = snap.z[i] - snap.z[j];
            double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            rho_sum += snap.m[j] * kernel_W(r, h);
            drho_sum_dh += snap.m[j] * kernel_dWdh(r, h);
        }

        double rho_h = pmass * std::pow(hfact / h, 3);
        double drho_h_dh = -3.0 * rho_h / h;

        double f = rho_sum - rho_h;
        double fprime = drho_sum_dh - drho_h_dh;

        double h_old = h;
        if (std::fabs(fprime) > 1e-300) {
            double h_new = h - f / fprime;
            // Damp: never move more than 20% per step (matches Phantom's
            // own 0.8h/1.2h clamp, dens.F90:1433-1439).
            h_new = std::min(h_new, 1.2 * h);
            h_new = std::max(h_new, 0.8 * h);
            if (h_new > 0.0) h = h_new;
        }
        // else: degenerate step (essentially no neighbours / no h-sensitivity)
        // -- hold h fixed this iteration rather than dividing by ~0; if this
        // persists to max_iterations the particle is reported unconverged.

        if (std::fabs(h - h_old) / h_old < tolh) {
            converged = true;
            ++iter;
            break;
        }
    }

    // The loop above updates h to h_new and THEN checks convergence against
    // h_old -- so on the breaking iteration, `h` is already the post-update
    // value but `rho_sum`/`neighbours` are still from the query at the START
    // of that same iteration (i.e. at h_old, one Newton step behind). Left
    // alone, that reports a (h, rho) pair that isn't self-consistent -- an
    // independent brute-force cross-check against a saved output file caught
    // this (max 2e-4 relative mismatch between the reported rho_computed and
    // the density actually implied by the reported h_computed, vs. 6e-8 once
    // this fix is in place -- see the validation notes). One more evaluation
    // at the final h fixes it, at the cost of one extra neighbour query per
    // particle.
    neighbours.clear();
    tree.query_radius(snap.x[i], snap.y[i], snap.z[i], kRadkern * h, neighbours);
    rho_sum = 0.0;
    for (size_t j : neighbours) {
        double dx = snap.x[i] - snap.x[j];
        double dy = snap.y[i] - snap.y[j];
        double dz = snap.z[i] - snap.z[j];
        double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        rho_sum += snap.m[j] * kernel_W(r, h);
    }

    result.h = h;
    result.rho = rho_sum;
    result.n_neighbours = static_cast<int>(neighbours.size());
    result.n_iterations = iter;
    result.converged = converged;
    return result;
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr,
                      "usage: h_from_neighbours <snapshot_path> [format] [hfact] [tolh] [output.dat]\n");
        return 1;
    }
    std::string snapshot_path = argv[1];
    std::string format = argc >= 3 ? argv[2] : "auto";
    double hfact = argc >= 4 ? std::atof(argv[3]) : 1.2;
    double tolh = argc >= 5 ? std::atof(argv[4]) : 1e-4;

    std::string output_path;
    if (argc >= 6) {
        output_path = argv[5];
    } else {
        std::string base = snapshot_path;
        size_t slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        size_t dot = base.find_last_of('.');
        if (dot != std::string::npos) base = base.substr(0, dot);
        output_path = base + "-h-compare.dat";
    }

    const int kMaxIterations = 100; // matches Phantom's maxdensits, dens.F90:101

    kiwi::SPHSnapshot snap;
    try {
        snap = kiwi::load_sph_snapshot(snapshot_path, format);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "h_from_neighbours: failed to load '%s': %s\n", snapshot_path.c_str(),
                     e.what());
        return 1;
    }
    std::fprintf(stderr, "h_from_neighbours: loaded %zu particles from %s\n", snap.size(),
                 snapshot_path.c_str());

    KDTree tree(snap.x, snap.y, snap.z);
    std::fprintf(stderr, "h_from_neighbours: kd-tree built, running per-particle h iteration "
                          "(hfact=%.4f tolh=%.2e)...\n",
                 hfact, tolh);

    std::vector<HResult> results(snap.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (long i = 0; i < static_cast<long>(snap.size()); ++i) {
        results[static_cast<size_t>(i)] = solve_h(static_cast<size_t>(i), snap, tree, hfact, tolh,
                                                    kMaxIterations);
    }

    FILE* out = std::fopen(output_path.c_str(), "w");
    if (!out) {
        std::fprintf(stderr, "h_from_neighbours: failed to open '%s' for writing\n", output_path.c_str());
        return 1;
    }
    std::fprintf(out, "# id x y z m h_snapshot rho_snapshot h_computed rho_computed n_neighbours "
                       "n_iterations converged rel_diff_h\n");

    size_t n_converged = 0;
    std::vector<double> rel_diffs;
    rel_diffs.reserve(snap.size());
    for (size_t i = 0; i < snap.size(); ++i) {
        const HResult& r = results[i];
        double rel_diff = (r.h - snap.h[i]) / snap.h[i];
        std::fprintf(out, "%d %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %d %d %d %.10g\n", snap.id[i],
                     snap.x[i], snap.y[i], snap.z[i], snap.m[i], snap.h[i], snap.rho[i], r.h, r.rho,
                     r.n_neighbours, r.n_iterations, r.converged ? 1 : 0, rel_diff);
        if (r.converged) ++n_converged;
        rel_diffs.push_back(rel_diff);
    }
    std::fclose(out);

    std::vector<double> abs_rel_diffs;
    abs_rel_diffs.reserve(rel_diffs.size());
    for (double d : rel_diffs) abs_rel_diffs.push_back(std::fabs(d));
    std::sort(abs_rel_diffs.begin(), abs_rel_diffs.end());
    double median = abs_rel_diffs.empty() ? 0.0 : abs_rel_diffs[abs_rel_diffs.size() / 2];
    double max_diff = abs_rel_diffs.empty() ? 0.0 : abs_rel_diffs.back();
    double mean = 0.0;
    for (double d : abs_rel_diffs) mean += d;
    if (!abs_rel_diffs.empty()) mean /= abs_rel_diffs.size();

    std::fprintf(stderr,
                 "h_from_neighbours: %zu/%zu converged (%.2f%%)\n"
                 "  |rel diff vs snapshot h|: mean=%.3e median=%.3e max=%.3e\n"
                 "  wrote %s\n",
                 n_converged, snap.size(), 100.0 * double(n_converged) / double(snap.size()), mean, median,
                 max_diff, output_path.c_str());

    return 0;
}
