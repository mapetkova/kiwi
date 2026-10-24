#include "kiwi/mapping.hpp"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <sys/stat.h>
#include <unordered_map>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "voro++.hh"
#include "kiwi/kernel_integral.hpp"
#include "kiwi/kernel_integrator.hpp"
#include "kiwi/point_set.hpp"
#include "kiwi/snapshot.hpp"
#include "kiwi/voro_grid.hpp"

namespace {
bool file_exists(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0;
}
} // namespace

namespace kiwi {

namespace {

// The face-vertex geometry a cell needs for the kernel integral, extracted
// once from voro++'s voronoicell_neighbor and kept as plain data afterward.
//
// voro++'s voronoicell_neighbor carries substantial internal bookkeeping
// beyond this (per-vertex-order edge/back-pointer arrays, duplicated
// neighbor-tracking arrays, scratch capacity from its growth strategy) --
// several times the footprint of the geometry alone. Retaining ~450k live
// voronoicell_neighbor objects (one per grid cell, needed simultaneously
// since the BFS from different particles revisits the same cells
// repeatedly) would push memory use well past what fits in a constrained
// environment. Extracting just f_vert/f_coord/n_list/volume once
// per cell and discarding the voro++ object cuts that down to the data
// compute_cell_integral() actually reads.
struct CellFaceGeometry {
    std::vector<int> f_vert;     // voro++ face_vertices() format: per face, [count, idx0, idx1, ...]
    std::vector<double> f_coord; // voro++ vertices() format: vertex coords relative to the site
    std::vector<int> n_list;     // one neighbor id per face, same order as f_vert's faces
    double volume = 0;
};

// Per-particle (thread-local, cleared before each particle's BFS) cache of
// already-evaluated face contributions, keyed by the unordered pair of cell
// indices on either side of the face: (min(a,b) << 32) | max(a,b).
//
// Within one particle's BFS, h and the particle position are fixed, so when
// two BFS-visited cells share a face, the two cells' Msum contributions from
// that shared face are exact negatives of each other -- not merely similar,
// but provably so: the shared face's normal flips sign between the two
// cells' own local vertex orderings (voro++ stores each cell's faces
// outward-facing relative to its own site), which flips r0's sign once;
// meanwhile the per-edge R_0 and the {phi1, phi2} pair are geometric
// (direction-independent) so their magnitudes are unaffected, and s1*s2
// picks up two sign flips (one from the reversed face winding via s2, one
// per edge from the reversed vertex order via s1) that cancel. So the only
// net sign change between the two sides is the one from r0 -- which is
// exactly the M = (s1*s2*r0<=0) flip. Whichever cell reaches the face first
// pays for the per-edge integral evaluation and caches the raw (signed, as
// seen from that cell) contribution; the second cell to reach the same face
// just negates it and skips the edge loop entirely.
//
// This also makes the BFS-discovery step redundant on the second visit: a
// cache hit means the neighbor across this face already computed and cached
// it, which can only happen if that neighbor was already dequeued and
// processed -- i.e. already marked visited -- so there is nothing left for
// the second side to discover through this face.
//
// Not applicable to domain-wall faces (negative neighbor ids, no real cell
// on the other side) or to a neighbor id that doesn't resolve through
// reverse_id_list -- those are computed directly every time, same as before
// this cache existed. Also assumes general position: two Voronoi cells
// share at most one face, true for a genuine (non-degenerate) point set,
// which real SPH particle positions are.
using FaceCache = std::unordered_map<int64_t, std::pair<int, double>>;

int64_t face_cache_key(int a, int b) {
    int lo = a < b ? a : b;
    int hi = a < b ? b : a;
    return (static_cast<int64_t>(lo) << 32) | static_cast<uint32_t>(hi);
}

// The shared evaluation of the SPH kernel line integral for one
// (particle, cell) pair. A slightly-negative Msum (a face-marching numerical
// artefact) is clamped to zero. Since this one value feeds both density and
// velocity, the clamp is applied unconditionally -- an unclamped small
// negative contribution would otherwise corrupt the density too.
//
// While walking this cell's faces, any unvisited neighbor still reachable by
// this particle's kernel (out_vert < number of face vertices) is appended to
// the BFS queue -- this is how the walk spreads from the particle's home
// cell out to every cell its smoothing kernel overlaps.
double compute_cell_integral(
    int ind, int cell_idx, double x0, double y0, double z0, double h,
    const KernelIntegrator& kernel,
    const std::vector<double>& x_list, const std::vector<double>& y_list, const std::vector<double>& z_list,
    const std::vector<CellFaceGeometry>& cell_geom,
    const std::vector<int>& reverse_id_list,
    std::vector<int>& kernel_list,
    std::vector<int>& next_list, int& next_list_top,
    FaceCache& face_cache) {

    const CellFaceGeometry& cell = cell_geom[cell_idx];
    const std::vector<int>& f_vert = cell.f_vert;
    const std::vector<double>& f_coord = cell.f_coord;
    const std::vector<int>& n_list = cell.n_list;
    const double x = x_list[cell_idx], y = y_list[cell_idx], z = z_list[cell_idx];

    double Msum = 0.0;
    int i = 0;

    auto it = f_vert.begin();
    do {
        auto iv = it;
        int face_vertex_count = *iv;

        // This face's neighbor and, if it resolves to a real cell in this
        // tessellation, the key the *other* side of this face would use to
        // look up this same evaluation -- see FaceCache's comment above for
        // why the two sides' contributions are provably exact negatives of
        // each other.
        int cell_ind = *(n_list.begin() + i);
        int r_cell_ind = -1;
        int64_t key = 0;
        bool cacheable = false;
        if (cell_ind >= 0) { // negative ids are domain walls, per voro++ convention
            r_cell_ind = reverse_id_list[cell_ind];
            if (r_cell_ind != -1) {
                cacheable = true;
                key = face_cache_key(cell_idx, r_cell_ind);
            }
        }

        if (cacheable) {
            auto cached = face_cache.find(key);
            if (cached != face_cache.end()) {
                // Second side: the owning cell already paid for this face's
                // edge loop, and (since a cache entry can only exist once
                // the owning cell has been visited) that neighbor is
                // already visited too, so there's nothing left to discover
                // through this face either. Skip straight to the next face.
                Msum += (cached->second.first == cell_idx) ? cached->second.second : -cached->second.second;
                it += face_vertex_count + 1;
                ++i;
                continue;
            }
        }

        // Extract this face's corners (voro++'s vertex-index indirection,
        // resolved to absolute coordinates) into the plain flattened form
        // compute_face_integral() takes -- the same per-face math this loop
        // used to do inline now lives there, shared with cartesian_grid.cpp.
        // Vertex order is preserved exactly (vertex k of the face, for
        // k = 0..face_vertex_count-1, same order voro++'s face_vertices()
        // already produces) -- unchanged from what this loop always fed
        // into the inline math before the refactor, so whatever convention
        // voro++'s actual output follows here is undisturbed either way.
        std::vector<double> face_verts;
        face_verts.reserve(static_cast<size_t>(face_vertex_count) * 3);
        for (int vi = 0; vi < face_vertex_count; ++vi) {
            int coord_idx = *(it + 1 + vi);
            face_verts.push_back(f_coord[3 * coord_idx] + x);
            face_verts.push_back(f_coord[3 * coord_idx + 1] + y);
            face_verts.push_back(f_coord[3 * coord_idx + 2] + z);
        }

        int out_vert = 0;
        double face_Msum = kernel.face(x0, y0, z0, h, face_verts, &out_vert);

        Msum += face_Msum;
        if (cacheable) {
            // Cache miss confirmed above (we only reach here when the
            // lookup found nothing), so this key isn't in the map yet.
            face_cache.emplace(key, std::make_pair(cell_idx, face_Msum));
        }

        if (out_vert != face_vertex_count) {
            if (cell_ind >= 0 && r_cell_ind != -1 && kernel_list[r_cell_ind] != ind) {
                next_list[next_list_top] = r_cell_ind;
                kernel_list[r_cell_ind] = ind;
                ++next_list_top;
            }
        }
        it += face_vertex_count + 1;
        ++i;
    } while (it < f_vert.end());

    if (Msum < 0.0) {
        if (Msum <= -0.01)
            std::cerr << "[kiwi] warning: negative mass fraction " << Msum
                      << " (particle " << ind << ", cell " << cell_idx << ")\n";
        Msum = 0.0;
    }

    return Msum;
}

// Tessellates `points` with voro++ and reduces each cell to the plain-data
// CellFaceGeometry compute_cell_integral() needs -- the setup both
// run_mapping_stage() and run_inverse_mapping_test_stage() share (geometry
// isn't persisted between stages, so each has to rebuild it fresh). Same
// deterministic container/block setup every time given the same `points`
// and cfg.grid_blocks, but callers should still look cells up by their
// voro++ site id via reverse_id_list (not assume index-for-index equality
// with any other run's own build of the same points).
//
// Returns the voro::container itself (heap-allocated -- it's not a movable
// or default-constructible type, so it can't be a plain return value or an
// out-param filled in place) because callers also need it afterward for
// find_voronoi_cell() in their own particle BFS loop; ownership passes to
// the caller.
std::unique_ptr<voro::container> build_cell_geometry(
        const Config& cfg, const PointSet& points,
        std::vector<int>& id_list,
        std::vector<double>& x_list, std::vector<double>& y_list, std::vector<double>& z_list,
        std::vector<CellFaceGeometry>& cell_geom,
        std::vector<int>& reverse_id_list) {
    using namespace voro;

    const Domain& d = cfg.domain;
    pre_container pcon(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max, false, false, false);
    for (size_t i = 0; i < points.size(); ++i) {
        pcon.put(points.id[i], points.x[i], points.y[i], points.z[i]);
    }
    int nx, ny, nz;
    guess_scaled_blocks(pcon, cfg.grid_blocks.block_scale_factor, nx, ny, nz);

    auto con = std::make_unique<container>(d.x_min, d.x_max, d.y_min, d.y_max, d.z_min, d.z_max,
                                            nx, ny, nz, false, false, false, 8);
    pcon.setup(*con);

    int max_id = 0;
    for (int pid : points.id) max_id = std::max(max_id, pid);
    reverse_id_list.assign(max_id + 1, -1);

    voronoicell_neighbor c;
    c_loop_all clo(*con);
    double x0, y0, z0;
    if (clo.start()) do {
        if (con->compute_cell(c, clo)) {
            int idx = static_cast<int>(cell_geom.size());
            id_list.push_back(clo.pid());
            reverse_id_list[clo.pid()] = idx;
            clo.pos(x0, y0, z0);
            x_list.push_back(x0);
            y_list.push_back(y0);
            z_list.push_back(z0);

            CellFaceGeometry geom;
            c.face_vertices(geom.f_vert);
            c.vertices(geom.f_coord);
            c.neighbors(geom.n_list);
            geom.volume = c.volume();
            cell_geom.push_back(std::move(geom));
        }
    } while (clo.inc());

    return con;
}

// Loads mapping.density_table from disk if set AND actually there;
// otherwise auto-generates the same canonical table in memory rather than
// failing the run -- covers both "not set at all" and "set but the file
// isn't found". Shared by run_mapping_stage() and
// run_inverse_mapping_test_stage() so both use the identical kernel
// evaluation method (tabulated vs. analytic) -- that's not just tidiness
// here: run_inverse_mapping_test_stage()'s per-particle coverage sum is only
// meaningful as a check on the SAME K_{a,i} values the forward mapping used
// to build M_i, so the two BFS passes have to evaluate the kernel integral
// the same way.
DensityLookupTable load_or_generate_density_table(const MappingConfig& mcfg) {
    DensityLookupTable table;
    if (!mcfg.use_tabulated_integral) {
        std::cerr << "[kiwi] using full_integral_new (analytic) -- no lookup table needed\n";
        return table;
    }
    if (!mcfg.density_table.empty() && file_exists(mcfg.density_table)) {
        std::cerr << "[kiwi] loading density lookup table from " << mcfg.density_table << "...\n";
        table = load_density_table(mcfg.density_table);
    } else {
        if (mcfg.density_table.empty()) {
            std::cerr << "[kiwi] mapping.density_table not set -- generating the kernel-integral "
                         "lookup table in memory instead of loading one from disk...\n";
        } else {
            std::cerr << "[kiwi] warning: mapping.density_table (" << mcfg.density_table
                       << ") not found on disk -- generating it in memory instead...\n";
        }
        table = generate_density_table();
    }
    std::cerr << "[kiwi] density lookup table ready (nr0=" << table.nr0 << " nR_0=" << table.nR_0
               << " nphi=" << table.nphi << ")\n";
    return table;
}

} // namespace

void save_mapping(const MappingResult& result, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("could not open mapping file for writing: " + path);
    f << std::setprecision(15);
    for (const auto& c : result.cells) {
        f << c.id << " " << c.x << " " << c.y << " " << c.z << " " << c.volume << " "
          << c.density << " " << c.vx << " " << c.vy << " " << c.vz << "\n";
    }
}

MappingResult load_mapping(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open mapping file: " + path);

    MappingResult result;
    CellMapping c;
    while (f >> c.id >> c.x >> c.y >> c.z >> c.volume >> c.density >> c.vx >> c.vy >> c.vz) {
        result.cells.push_back(c);
    }
    if (result.cells.empty())
        throw std::runtime_error("mapping file contained no rows (wrong format?): " + path);
    return result;
}

MappingResult run_mapping_stage(const Config& cfg) {
    using namespace voro;

    std::string final_path = cfg.stage_path("mapping.dat");
    if (file_exists(final_path)) {
        return load_mapping(final_path);
    }

    DensityLookupTable table = load_or_generate_density_table(cfg.mapping);
    const KernelIntegrator kernel = make_kernel_integrator(cfg, cfg.mapping.use_tabulated_integral, &table);

    PointSet points = run_point_set_stage(cfg);
    SPHSnapshot snap = load_sph_snapshot(cfg.input_snapshot, cfg.input_snapshot_format, cfg.kernel_hfact(), cfg.phantom_units,
                                          cfg.allow_missing_velocity);
    std::cerr << "[kiwi] loaded " << points.size() << " grid points, " << snap.size() << " SPH particles\n";

    // One voronoicell_neighbor is computed at a time and immediately reduced
    // to a CellFaceGeometry (plain data) -- see that struct's comment for
    // why we don't keep ~450k live voro++ cell objects in memory at once.
    std::vector<int> id_list;
    std::vector<double> x_list, y_list, z_list;
    std::vector<CellFaceGeometry> cell_geom;
    id_list.reserve(points.size());
    x_list.reserve(points.size());
    y_list.reserve(points.size());
    z_list.reserve(points.size());
    cell_geom.reserve(points.size());
    std::vector<int> reverse_id_list;

    std::unique_ptr<container> con_ptr =
        build_cell_geometry(cfg, points, id_list, x_list, y_list, z_list, cell_geom, reverse_id_list);
    container& con = *con_ptr;
    int max_id = static_cast<int>(reverse_id_list.size()) - 1;

    const size_t n_cells = cell_geom.size();

    // The particle loop below is embarrassingly parallel: each particle's
    // BFS walk only reads shared, read-only state (x_list/y_list/z_list,
    // cell_geom, reverse_id_list, con, table) and writes to per-cell
    // accumulators (dens_cells/v{x,y,z}_cells) that different particles'
    // kernels can legitimately overlap on. Rather than lock those shared
    // accumulators (contention would erase most of the parallel speedup,
    // since nearby particles hit the same handful of cells), each thread
    // gets its own full-size accumulator + its own kernel_list/next_list
    // BFS scratch space; buffers are summed together after the parallel
    // region. Cost: O(n_cells * n_threads) extra memory -- a few hundred MB
    // at most for realistic thread counts, trivial next to the input data.
    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    std::cerr << "[kiwi] mapping: running with " << n_threads << " thread(s)\n";

    std::vector<std::vector<double>> dens_cells_t(n_threads, std::vector<double>(n_cells, 0.0));
    std::vector<std::vector<double>> vx_cells_t(n_threads, std::vector<double>(n_cells, 0.0));
    std::vector<std::vector<double>> vy_cells_t(n_threads, std::vector<double>(n_cells, 0.0));
    std::vector<std::vector<double>> vz_cells_t(n_threads, std::vector<double>(n_cells, 0.0));
    std::vector<std::vector<int>> kernel_list_t(n_threads, std::vector<int>(n_cells, 0));
    std::vector<std::vector<int>> next_list_t(n_threads, std::vector<int>(n_cells, 0));
    // One FaceCache per thread, cleared before each particle's BFS -- see
    // FaceCache's comment. Reused (cleared, not reallocated) across
    // particles within a thread, same spirit as kernel_list/next_list.
    std::vector<FaceCache> face_cache_t(n_threads);

    double vx_mean = 0, vy_mean = 0, vz_mean = 0;
    long unresolved = 0; // particle position didn't resolve to a known cell -- skipped
    long unmapped = 0;   // particle's kernel never reached beyond its own home cell
    const long n_particles = static_cast<long>(snap.size());

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256) reduction(+:vx_mean,vy_mean,vz_mean,unresolved,unmapped)
#endif
    for (long idx = 0; idx < n_particles; ++idx) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<double>& dens_cells = dens_cells_t[tid];
        std::vector<double>& vx_cells = vx_cells_t[tid];
        std::vector<double>& vy_cells = vy_cells_t[tid];
        std::vector<double>& vz_cells = vz_cells_t[tid];
        std::vector<int>& kernel_list = kernel_list_t[tid];
        std::vector<int>& next_list = next_list_t[tid];
        FaceCache& face_cache = face_cache_t[tid];

        int ind = snap.id[idx];
        double x0 = snap.x[idx], y0 = snap.y[idx], z0 = snap.z[idx];
        double vx = snap.vx[idx], vy = snap.vy[idx], vz = snap.vz[idx];
        double m = snap.m[idx], h = snap.h[idx];

        vx_mean += vx;
        vy_mean += vy;
        vz_mean += vz;

        if (h <= 0.0) continue; // degenerate/invalid smoothing length -- no contribution

        double fx, fy, fz;
        int site_id;
        // If a particle lies outside the domain box, find_voronoi_cell()
        // returns false *without* setting site_id, so the return value must
        // be checked before indexing reverse_id_list[] (an uninitialized
        // site_id would be undefined behavior). Particles outside the box
        // genuinely have no cell to contribute to, so they are skipped.
        bool found = con.find_voronoi_cell(x0, y0, z0, fx, fy, fz, site_id);
        int k = (found && site_id >= 0 && site_id <= max_id) ? reverse_id_list[site_id] : -1;
        if (k == -1) {
            ++unresolved;
            continue;
        }

        next_list[0] = k;
        kernel_list[k] = ind;
        int next_list_top = 1;
        int next_list_bottom = 0;
        face_cache.clear(); // fresh per particle -- h and (x0,y0,z0) just changed

        while (next_list_bottom < next_list_top) {
            int cell_idx = next_list[next_list_bottom];
            double integral = compute_cell_integral(ind, cell_idx, x0, y0, z0, h, kernel,
                                                      x_list, y_list, z_list, cell_geom, reverse_id_list,
                                                      kernel_list, next_list, next_list_top, face_cache);
            dens_cells[cell_idx] += integral * m;
            vx_cells[cell_idx] += integral * m * vx;
            vy_cells[cell_idx] += integral * m * vy;
            vz_cells[cell_idx] += integral * m * vz;

            ++next_list_bottom;
        }
        if (next_list_top == 1) ++unmapped;
    }

    std::cerr << "[kiwi] mapping: reducing per-thread accumulators...\n";
    std::vector<double> dens_cells(n_cells, 0.0), vx_cells(n_cells, 0.0), vy_cells(n_cells, 0.0), vz_cells(n_cells, 0.0);
    for (int t = 0; t < n_threads; ++t) {
        for (size_t i = 0; i < n_cells; ++i) {
            dens_cells[i] += dens_cells_t[t][i];
            vx_cells[i] += vx_cells_t[t][i];
            vy_cells[i] += vy_cells_t[t][i];
            vz_cells[i] += vz_cells_t[t][i];
        }
    }

    if (snap.size() > 0) {
        vx_mean /= static_cast<double>(snap.size());
        vy_mean /= static_cast<double>(snap.size());
        vz_mean /= static_cast<double>(snap.size());
    }

    std::cerr << "[kiwi] mapping: " << unresolved << " particle(s) did not resolve to a known cell (skipped), "
              << unmapped << " particle(s)' kernel never reached beyond their home cell\n";

    MappingResult result;
    result.cells.reserve(n_cells);
    for (size_t i = 0; i < n_cells; ++i) {
        CellMapping cm;
        cm.id = id_list[i];
        cm.x = x_list[i];
        cm.y = y_list[i];
        cm.z = z_list[i];
        cm.volume = cell_geom[i].volume;
        cm.density = (cm.volume > 0.0) ? dens_cells[i] / cm.volume : 0.0;
        if (dens_cells[i] > 0.0) {
            cm.vx = vx_cells[i] / dens_cells[i];
            cm.vy = vy_cells[i] / dens_cells[i];
            cm.vz = vz_cells[i] / dens_cells[i];
        } else {
            cm.vx = vx_mean;
            cm.vy = vy_mean;
            cm.vz = vz_mean;
        }
        result.cells.push_back(cm);
    }

    save_mapping(result, final_path);
    return result;
}

void save_inverse_mapping_test(const InverseMappingTestResult& result, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("could not open inverse-mapping-test file for writing: " + path);
    f << std::setprecision(15);
    f << "# total_original_mass total_grid_mass total_recovered_mass\n";
    f << result.total_original_mass << " " << result.total_grid_mass << " " << result.total_recovered_mass << "\n";
    for (const auto& p : result.particles) {
        f << p.id << " " << p.m_original << " " << p.m_recovered << " " << p.coverage_fraction << "\n";
    }
}

InverseMappingTestResult load_inverse_mapping_test(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("could not open inverse-mapping-test file: " + path);

    std::string header_line;
    if (!std::getline(f, header_line) || header_line.empty() || header_line[0] != '#')
        throw std::runtime_error("inverse-mapping-test file missing expected '#' totals header: " + path);

    InverseMappingTestResult result;
    if (!(f >> result.total_original_mass >> result.total_grid_mass >> result.total_recovered_mass))
        throw std::runtime_error("could not read totals line from inverse-mapping-test file: " + path);

    InverseMappingTestRow row;
    while (f >> row.id >> row.m_original >> row.m_recovered >> row.coverage_fraction) {
        result.particles.push_back(row);
    }
    if (result.particles.empty())
        throw std::runtime_error("inverse-mapping-test file contained no particle rows (wrong format?): " + path);
    return result;
}

InverseMappingTestResult run_inverse_mapping_test_stage(const Config& cfg) {
    using namespace voro;

    std::string final_path = cfg.stage_path("inverse-mapping-test.dat");
    if (file_exists(final_path)) {
        return load_inverse_mapping_test(final_path);
    }

    // Reuses the forward mapping's own cached output (or triggers it to be
    // computed, same as any other stage chaining from an earlier one) as the
    // per-cell mass denominator -- this is the "testing step" case, so no
    // second BFS pass is needed to build that denominator from scratch.
    MappingResult forward = run_mapping_stage(cfg);

    std::unordered_map<int, double> cell_mass_by_id;
    cell_mass_by_id.reserve(forward.cells.size() * 2);
    double total_grid_mass = 0.0;
    for (const auto& c : forward.cells) {
        double m_i = c.density * c.volume;
        cell_mass_by_id.emplace(c.id, m_i);
        total_grid_mass += m_i;
    }

    DensityLookupTable table = load_or_generate_density_table(cfg.mapping);
    const KernelIntegrator kernel = make_kernel_integrator(cfg, cfg.mapping.use_tabulated_integral, &table);

    PointSet points = run_point_set_stage(cfg);
    SPHSnapshot snap = load_sph_snapshot(cfg.input_snapshot, cfg.input_snapshot_format, cfg.kernel_hfact(), cfg.phantom_units,
                                          cfg.allow_missing_velocity);
    std::cerr << "[kiwi] inverse-mapping-test: loaded " << points.size() << " grid points, "
              << snap.size() << " SPH particles\n";

    // Fresh tessellation -- see build_cell_geometry()'s comment for why cell
    // identity is looked up by voro++ site id (cell_mass_by_id, keyed by
    // c.id above) rather than assumed to line up index-for-index with
    // `forward`'s own cell ordering.
    std::vector<int> id_list;
    std::vector<double> x_list, y_list, z_list;
    std::vector<CellFaceGeometry> cell_geom;
    id_list.reserve(points.size());
    x_list.reserve(points.size());
    y_list.reserve(points.size());
    z_list.reserve(points.size());
    cell_geom.reserve(points.size());
    std::vector<int> reverse_id_list;

    std::unique_ptr<container> con_ptr =
        build_cell_geometry(cfg, points, id_list, x_list, y_list, z_list, cell_geom, reverse_id_list);
    container& con = *con_ptr;
    int max_id = static_cast<int>(reverse_id_list.size()) - 1;

    int n_threads = 1;
#ifdef _OPENMP
    n_threads = omp_get_max_threads();
#endif
    std::cerr << "[kiwi] inverse-mapping-test: running with " << n_threads << " thread(s)\n";

    std::vector<std::vector<int>> kernel_list_t(n_threads, std::vector<int>(cell_geom.size(), 0));
    std::vector<std::vector<int>> next_list_t(n_threads, std::vector<int>(cell_geom.size(), 0));
    std::vector<FaceCache> face_cache_t(n_threads);

    const long n_particles = static_cast<long>(snap.size());
    std::vector<InverseMappingTestRow> particles(n_particles);
    long unresolved = 0;         // particle position didn't resolve to a known cell -- skipped, same as forward
    long missing_cell_mass = 0;  // BFS-visited cell not found in the forward mapping's own cell set -- staleness/cache-mismatch signal

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 256) reduction(+:unresolved,missing_cell_mass)
#endif
    for (long idx = 0; idx < n_particles; ++idx) {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        std::vector<int>& kernel_list = kernel_list_t[tid];
        std::vector<int>& next_list = next_list_t[tid];
        FaceCache& face_cache = face_cache_t[tid];

        int ind = snap.id[idx];
        double x0 = snap.x[idx], y0 = snap.y[idx], z0 = snap.z[idx];
        double m = snap.m[idx], h = snap.h[idx];

        InverseMappingTestRow row;
        row.id = ind;
        row.m_original = m;
        row.m_recovered = 0.0;
        row.coverage_fraction = 0.0;

        if (h <= 0.0) { // degenerate/invalid smoothing length -- forward mapping skipped this particle too
            particles[idx] = row;
            continue;
        }

        double fx, fy, fz;
        int site_id;
        bool found = con.find_voronoi_cell(x0, y0, z0, fx, fy, fz, site_id);
        int k = (found && site_id >= 0 && site_id <= max_id) ? reverse_id_list[site_id] : -1;
        if (k == -1) {
            ++unresolved;
            particles[idx] = row;
            continue;
        }

        next_list[0] = k;
        kernel_list[k] = ind;
        int next_list_top = 1;
        int next_list_bottom = 0;
        face_cache.clear();

        double recovered = 0.0;
        while (next_list_bottom < next_list_top) {
            int cell_idx = next_list[next_list_bottom];
            double integral = compute_cell_integral(ind, cell_idx, x0, y0, z0, h, kernel,
                                                      x_list, y_list, z_list, cell_geom, reverse_id_list,
                                                      kernel_list, next_list, next_list_top, face_cache);
            // weight_{a,i} * P_i, with P_i = M_i (see this stage's header
            // comment for why that makes M_i cancel algebraically -- kept
            // explicit here anyway, both for fidelity to the general
            // CMacIonize-derived formula and because a future property other
            // than mass itself would need the real division).
            auto it = cell_mass_by_id.find(id_list[cell_idx]);
            if (it != cell_mass_by_id.end() && it->second > 0.0) {
                double m_i = it->second;
                recovered += (integral * m / m_i) * m_i;
            } else if (it == cell_mass_by_id.end()) {
                ++missing_cell_mass;
            }
            ++next_list_bottom;
        }

        row.m_recovered = recovered;
        row.coverage_fraction = (m > 0.0) ? recovered / m : 0.0;
        particles[idx] = row;
    }

    if (missing_cell_mass > 0) {
        std::cerr << "[kiwi] warning: inverse-mapping-test hit " << missing_cell_mass
                  << " BFS-visited cell(s) not present in the forward mapping's own cell set -- "
                     "mapping.dat and points-final.dat under this run directory may be out of sync "
                     "(check_manifest_consistency() should normally catch this; consider deleting the "
                     "run directory and re-running from scratch if this looks wrong)\n";
    }
    std::cerr << "[kiwi] inverse-mapping-test: " << unresolved
              << " particle(s) did not resolve to a known cell (skipped)\n";

    InverseMappingTestResult result;
    result.particles = std::move(particles);
    result.total_grid_mass = total_grid_mass;
    result.total_original_mass = 0.0;
    result.total_recovered_mass = 0.0;
    for (const auto& p : result.particles) {
        result.total_original_mass += p.m_original;
        result.total_recovered_mass += p.m_recovered;
    }

    std::cerr << "[kiwi] inverse-mapping-test: total original mass = " << result.total_original_mass
              << ", total grid mass = " << result.total_grid_mass
              << ", total recovered mass = " << result.total_recovered_mass << "\n";

    save_inverse_mapping_test(result, final_path);
    return result;
}

} // namespace kiwi
