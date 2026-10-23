"""Decode a legacy POLARIS Voronoi grid file (grid-file.cc's binary format)
for validation purposes. Not part of the tool itself -- a throwaway
cross-check script.

Format (see grid-file.cc / GridVoronoi.cpp analysis):
  header: unsigned short sim, unsigned short params, unsigned short plist[params],
          double cells_d, double len
  per cell (rows are implicitly ordered by ascending original particle id,
  skipping any id that failed to produce a cell -- there is no id field):
          float x0, float y0, float z0, double vol,
          float <params values, in plist order>,
          int neigh_count, int neigh_count * neighbor_id
"""
import struct
import sys


def decode(path):
    with open(path, "rb") as f:
        data = f.read()

    off = 0
    sim, params = struct.unpack_from("<HH", data, off); off += 4
    plist = struct.unpack_from(f"<{params}H", data, off); off += 2 * params
    cells_d, length_scale = struct.unpack_from("<dd", data, off); off += 16

    cells = []
    neighbors = []  # per cell: neighbour row indices (>= 0) or domain-wall markers (< 0)
    n = 0
    while off < len(data):
        x0, y0, z0 = struct.unpack_from("<fff", data, off); off += 12
        (vol,) = struct.unpack_from("<d", data, off); off += 8
        vals = struct.unpack_from(f"<{params}f", data, off); off += 4 * params
        (neigh,) = struct.unpack_from("<i", data, off); off += 4
        neighbor_ids = struct.unpack_from(f"<{neigh}i", data, off); off += 4 * neigh
        cells.append((x0, y0, z0, vol, vals, neigh))
        neighbors.append(neighbor_ids)
        n += 1

    return {"sim": sim, "params": params, "plist": plist, "cells_d": cells_d,
            "length_scale": length_scale, "cells": cells, "neighbors": neighbors}


if __name__ == "__main__":
    info = decode(sys.argv[1])
    print("sim:", info["sim"], "params:", info["params"], "plist:", info["plist"])
    print("cells_d:", info["cells_d"], "length_scale (box x-extent, m):", info["length_scale"])
    print("num cells decoded:", len(info["cells"]))
    print("first cell:", info["cells"][0])
