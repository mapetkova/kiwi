// Generic plain-text Voronoi export. Not tied to POLARIS or any
// particular downstream consumer -- one row per cell, human-readable, easy
// to load with numpy/pandas for a quick look or a custom pipeline.
//
// Format: a single header line
//   # id x y z volume density vx vy vz n_neighbors neighbor_ids...
// followed by one line per cell:
//   <id> <x> <y> <z> <volume> <density> <vx> <vy> <vz> <n_neighbors> <id...>
// Positions/volume/velocity are in the run's native (unscaled) units --
// whatever cfg.input_snapshot used -- since this format has no fixed
// downstream unit convention to match, unlike the POLARIS writer.
// Neighbor ids are the *cell ids* (not row indices) of neighboring cells,
// or a negative voro++ domain-wall marker; only ids present in this file's
// own output are meaningful as cross-references. A neighbor id that has no
// corresponding row (missing tessellation or mapping data for that
// particle) is dropped and counted, same as the POLARIS writer.
#pragma once

#include <string>

#include "kiwi/config.hpp"

namespace kiwi {

struct VoronoiWriteStats {
    size_t cells_written = 0;
    size_t missing_cells = 0;
    size_t dropped_neighbor_refs = 0;
};

// Runs (or loads, cached) the point-set, tessellation, and mapping stages,
// joins tessellation geometry with mapping density/velocity by cell id, and
// writes the result to cfg.stage_path("grid-voronoi.txt").
VoronoiWriteStats write_voronoi_text(const Config& cfg);

} // namespace kiwi
