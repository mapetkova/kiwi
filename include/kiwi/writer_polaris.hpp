// Stage 5a: POLARIS binary grid writer. Reproduces the exact binary format
// grid-file.cc wrote (see tests/decode_polaris_grid.py for the byte layout
// this was reverse-engineered from and validated against).
//
// Differs from the legacy writer in one deliberate way: neighbor ids are
// translated through an explicit id -> output-row-index map built from the
// cells actually being written, rather than the legacy code's `id - 1`
// assumption (which silently produces wrong neighbor references whenever
// any site in [1, max_id] failed to produce a tessellated cell -- rare, but
// confirmed to happen). A neighbor reference to a
// missing cell is dropped and counted, never miswired to the wrong row.
#pragma once

#include <string>

#include "kiwi/config.hpp"

namespace kiwi {

struct PolarisWriteStats {
    size_t cells_written = 0;
    size_t missing_cells = 0;       // ids present in the point set but absent from the tessellation/mapping join
    size_t dropped_neighbor_refs = 0; // positive neighbor ids that pointed at a missing cell
};

// Runs (or loads, cached) the point-set, tessellation, and mapping stages,
// joins tessellation geometry with mapping density/velocity by cell id, and
// writes the POLARIS binary grid to cfg.stage_path("grid-polaris.dat").
PolarisWriteStats write_polaris_grid(const Config& cfg);

} // namespace kiwi
