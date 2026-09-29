"""Cross-check a POLARIS grid file written by `kiwi write-polaris` against:
  1. a legacy grid-file.cc binary output, for geometry only (position,
     volume, neighbor count/topology). Density/velocity are NOT compared
     against the legacy file here, because grid-file.cc's toggle-by-mutation
     history means you can't always trust that its embedded density/velocity
     payload was actually produced from the same run as its geometry (this
     is exactly the bug class the whole project started from -- see
     manifest.json/README for the tool's own guard against it). Geometry
     only depends on the coordinate file, which is easier to trust.
  2. this tool's own mapping.dat from the same run directory, for
     density/velocity, scaled by the run's polaris_writer config.

Usage:
  python3 validate_writer_polaris.py <legacy_grid.dat> <new_grid-polaris.dat> <mapping.dat> \\
      [length_scale_m] [velocity_scale] [density_scale]

The three scale factors default to this tool's polaris_writer defaults (see
config.hpp) -- pass your own if you ran write-polaris with a non-default
polaris_writer config.

Not part of the tool itself -- a standalone validation script.
"""
import sys

sys.path.insert(0, "tests")
from decode_polaris_grid import decode

DEFAULT_LENGTH_SCALE = 3.08567758e16
DEFAULT_VEL_SCALE = 1000.0


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        sys.exit(1)

    legacy_path, new_path, mapping_path = sys.argv[1], sys.argv[2], sys.argv[3]
    length_scale = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_LENGTH_SCALE
    vel_scale = float(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_VEL_SCALE
    dens_scale = float(sys.argv[6]) if len(sys.argv) > 6 else 1.9891e30 / length_scale**3

    legacy = decode(legacy_path)
    new = decode(new_path)

    print(f"legacy cells: {len(legacy['cells'])}, new cells: {len(new['cells'])}")
    assert legacy["cells_d"] == new["cells_d"] == len(legacy["cells"]) == len(new["cells"]), \
        "cell counts don't match -- these two files are probably from different runs"

    len_relerr = abs(legacy["length_scale"] - new["length_scale"]) / legacy["length_scale"]
    print(f"header length_scale relative difference: {len_relerr:.6e}")
    assert legacy["plist"] == new["plist"]
    assert legacy["sim"] == new["sim"]

    # --- geometry: position + volume, row for row ---
    max_pos_err = 0.0
    max_vol_relerr = 0.0
    neigh_count_mismatches = 0
    for lc, nc in zip(legacy["cells"], new["cells"]):
        lx0, ly0, lz0, lvol, _lvals, lneigh = lc
        nx0, ny0, nz0, nvol, _nvals, nneigh = nc
        max_pos_err = max(max_pos_err, abs(lx0 - nx0), abs(ly0 - ny0), abs(lz0 - nz0))
        if lvol != 0:
            max_vol_relerr = max(max_vol_relerr, abs(lvol - nvol) / abs(lvol))
        if lneigh != nneigh:
            neigh_count_mismatches += 1

    print(f"max abs position error (m): {max_pos_err:.6e}")
    print(f"max relative volume error: {max_vol_relerr:.6e}")
    print(f"neighbor-count mismatches: {neigh_count_mismatches} / {len(legacy['cells'])}")

    # --- density/velocity: writer output vs. mapping.dat, scaled ---
    mapping_by_id = {}
    with open(mapping_path) as f:
        for line in f:
            parts = line.split()
            pid = int(parts[0])
            _vol, dens, vx, vy, vz = map(float, parts[4:9])
            mapping_by_id[pid] = (dens, vx, vy, vz)

    # The writer's row order is strictly ascending cell id, and mapping.dat's
    # ids are exactly the tessellated cell ids it was joined against, so
    # sorting mapping.dat's ids reconstructs the same sequence -- valid as
    # long as the writer's own "missing cells" count was 0 for this run
    # (check the write-polaris stderr log if this assertion trips).
    ordered_ids = sorted(mapping_by_id.keys())
    assert len(ordered_ids) == len(new["cells"]), (
        f"{len(ordered_ids)} mapping.dat ids vs {len(new['cells'])} written cells -- "
        "the writer likely dropped some cells (see its stderr log); can't align by "
        "position without a matching count"
    )

    max_dens_relerr = 0.0
    max_vel_abserr = 0.0
    for pid, (_x0, _y0, _z0, _vol, vals, _neigh) in zip(ordered_ids, new["cells"]):
        dens, vx, vy, vz = mapping_by_id[pid]
        _Tg, wvx, wvy, wvz, wdens = vals
        expected_dens = dens * dens_scale
        if expected_dens != 0:
            max_dens_relerr = max(max_dens_relerr, abs(wdens - expected_dens) / abs(expected_dens))
        max_vel_abserr = max(
            max_vel_abserr,
            abs(wvx - vx * vel_scale), abs(wvy - vy * vel_scale), abs(wvz - vz * vel_scale),
        )

    print(f"max relative density error (writer vs mapping.dat, scaled): {max_dens_relerr:.6e}")
    print(f"max abs velocity error (m/s, writer vs mapping.dat, scaled): {max_vel_abserr:.6e}")
    print("OK")


if __name__ == "__main__":
    main()
