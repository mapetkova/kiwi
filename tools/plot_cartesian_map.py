#!/usr/bin/env python3
"""Load and plot a kiwi Cartesian grid: cartesian_map.dat (from `kiwi cartesian-map`) or
cartesian_grid.dat (from `kiwi sph-to-grid`); the format is recognised from the file.

Usage:
    plot_cartesian_map.py <cartesian_map.dat | cartesian_grid.dat> [--field density|vx|vy|vz] [--axis x|y|z]
                           [--mode slice|projection] [--index N] [--out out.png] [--log]
                           [--cmap NAME] [--vmin V] [--vmax V]

Two display modes (--mode, default "slice"):

  slice (default): a single-pixel-thick slice through the resampled
    Cartesian grid at --index along --axis (mid-plane by default). --index
    is ignored (with a note) in projection mode.

  projection: collapses the whole grid along --axis instead of picking one
    slice. For --field density this is the column density -- density summed
    along --axis and multiplied by that axis's cell width, i.e. sum(rho)*dz
    (the standard line-of-sight integral, not a bare unweighted sum). For
    --field {vx,vy,vz} it's the density-weighted average of that velocity
    component along --axis, i.e. sum(rho*v)/sum(rho) per pixel (nan where a
    column has zero density everywhere).

Colormap: defaults to 'coolwarm' (diverging, reads naturally around zero)
for --field {vx,vy,vz}, and matplotlib's own default ('viridis') for
density. Override either with --cmap.

Colorbar range: --vmin/--vmax fix the colorbar's lower/upper bound explicitly
(either or both -- the other end is left to matplotlib's own auto-scaling).
Applied to the plotted data as-is, i.e. AFTER --log if both are given -- pass
log10-scaled bounds (e.g. --vmin -2 for a floor of 1e-2) when combining the
two.

This is meant as a quick look, not a finished figure -- read the format
below and pull the arrays into your own plotting code (or a Jupyter
notebook) for anything more involved.

Two file formats (binary, little-endian on any machine this is likely to run on), both with the
density and velocity arrays flattened x-major, (ix*ny+iy)*nz+iz:

cartesian_map.dat, magic "VGTCMAP1" (cubic grid, one resolution for all three axes):
    8 bytes   magic "VGTCMAP1"
    int32     resolution (cells per axis, same for x/y/z)
    6 x f64   domain bounds: x_min, x_max, y_min, y_max, z_min, z_max
    N x f64   density  (N = resolution**3, flattened (ix*res+iy)*res+iz)
    N x f64   vx
    N x f64   vy
    N x f64   vz
    N x i32   voronoi_cell_id (the Voronoi cell each Cartesian cell resolved
              to; 0 means the search found no enclosing cell there)

cartesian_grid.dat, magic "VGTCGRD2" (per-axis resolution, no cell ids):
    8 bytes   magic "VGTCGRD2"
    3 x int32 nx, ny, nz
    6 x f64   domain bounds, as above
    N x f64   density, vx, vy, vz  (N = nx*ny*nz)

A minimal loader for cartesian_map.dat, with no dependency beyond numpy, looks like this:

    import numpy as np

    def load_cartesian_map(path):
        with open(path, "rb") as f:
            magic = f.read(8)
            assert magic == b"VGTCMAP1", "not a kiwi cartesian_map.dat file"
            res = int(np.fromfile(f, dtype="<i4", count=1)[0])
            bounds = np.fromfile(f, dtype="<f8", count=6)
            n = res ** 3
            density = np.fromfile(f, dtype="<f8", count=n).reshape(res, res, res)
            vx = np.fromfile(f, dtype="<f8", count=n).reshape(res, res, res)
            vy = np.fromfile(f, dtype="<f8", count=n).reshape(res, res, res)
            vz = np.fromfile(f, dtype="<f8", count=n).reshape(res, res, res)
            cell_id = np.fromfile(f, dtype="<i4", count=n).reshape(res, res, res)
        return res, bounds, density, vx, vy, vz, cell_id
"""
import argparse
import sys


def load_cartesian_map(path):
    """Read either format. Returns (dims, bounds, density, vx, vy, vz, cell_id) where dims is
    (nx, ny, nz), the arrays have shape dims, and cell_id is None for a cartesian_grid.dat file."""
    import numpy as np

    with open(path, "rb") as f:
        magic = f.read(8)
        if magic == b"VGTCMAP1":
            res = int(np.fromfile(f, dtype="<i4", count=1)[0])
            dims = (res, res, res)
        elif magic == b"VGTCGRD2":
            dims = tuple(int(v) for v in np.fromfile(f, dtype="<i4", count=3))
        else:
            raise ValueError(f"{path}: not a kiwi cartesian_map.dat / cartesian_grid.dat file (bad magic)")
        bounds = np.fromfile(f, dtype="<f8", count=6)
        n = dims[0] * dims[1] * dims[2]
        density, vx, vy, vz = (np.fromfile(f, dtype="<f8", count=n).reshape(dims) for _ in range(4))
        cell_id = np.fromfile(f, dtype="<i4", count=n).reshape(dims) if magic == b"VGTCMAP1" else None
    return dims, bounds, density, vx, vy, vz, cell_id


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cartesian_map", help="path to cartesian_map.dat or cartesian_grid.dat")
    ap.add_argument("--field", choices=["density", "vx", "vy", "vz"], default="density")
    ap.add_argument("--axis", choices=["x", "y", "z"], default="z",
                     help="axis to slice/project along (default: z, i.e. an x-y view)")
    ap.add_argument("--mode", choices=["slice", "projection"], default="slice",
                     help="'slice' (default): a single-pixel-thick slice at --index. "
                          "'projection': collapse along --axis -- column density for "
                          "--field density, density-weighted average velocity for vx/vy/vz "
                          "(see module docstring). --index is ignored in this mode.")
    ap.add_argument("--index", type=int, default=None,
                     help="[--mode slice] grid index to slice at along --axis (default: the middle)")
    ap.add_argument("--out", default=None, help="output PNG (default: shows interactively)")
    ap.add_argument("--log", action="store_true", help="log10-scale the field (density only, "
                     "clips non-positive values)")
    ap.add_argument("--cmap", default=None,
                     help="matplotlib colormap name (default: 'coolwarm' for --field "
                          "vx/vy/vz, matplotlib's own default for density)")
    ap.add_argument("--vmin", type=float, default=None,
                     help="colorbar lower bound (default: auto-scaled from the plotted data). "
                          "Applied after --log if both are given.")
    ap.add_argument("--vmax", type=float, default=None,
                     help="colorbar upper bound (default: auto-scaled from the plotted data). "
                          "Applied after --log if both are given.")
    args = ap.parse_args()

    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        print(f"plot_cartesian_map.py: --vmin ({args.vmin}) must be less than --vmax "
              f"({args.vmax})", file=sys.stderr)
        return 2

    try:
        import numpy as np
        import matplotlib
        if args.out:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"plot_cartesian_map.py: numpy and matplotlib are required: {e}", file=sys.stderr)
        return 3

    dims, bounds, density, vx, vy, vz, cell_id = load_cartesian_map(args.cartesian_map)
    x_min, x_max, y_min, y_max, z_min, z_max = bounds
    axis_bounds = {"x": (x_min, x_max), "y": (y_min, y_max), "z": (z_min, z_max)}

    field_by_name = {"density": density, "vx": vx, "vy": vy, "vz": vz}
    field = field_by_name[args.field]

    axis_idx = {"x": 0, "y": 1, "z": 2}[args.axis]
    n_axis = dims[axis_idx]  # cells along the axis being sliced/projected

    if axis_idx == 0:
        extent = (y_min, y_max, z_min, z_max)
        xlabel, ylabel = "y", "z"
    elif axis_idx == 1:
        extent = (x_min, x_max, z_min, z_max)
        xlabel, ylabel = "x", "z"
    else:
        extent = (x_min, x_max, y_min, y_max)
        xlabel, ylabel = "x", "y"

    if args.mode == "slice":
        index = args.index if args.index is not None else n_axis // 2
        if not (0 <= index < n_axis):
            print(f"plot_cartesian_map.py: --index {index} out of range [0, {n_axis})", file=sys.stderr)
            return 2

        if axis_idx == 0:
            slice_2d = field[index, :, :]
        elif axis_idx == 1:
            slice_2d = field[:, index, :]
        else:
            slice_2d = field[:, :, index]

        plot_data = slice_2d.T  # imshow expects (row, col) = (y-like, x-like)
        title = f"{{label}}, {args.axis}-slice at index {index}/{n_axis}"
    else:
        if args.index is not None:
            print("plot_cartesian_map.py: --index is ignored in --mode projection",
                  file=sys.stderr)

        axis_lo, axis_hi = axis_bounds[args.axis]
        dz = (axis_hi - axis_lo) / n_axis

        if args.field == "density":
            proj_2d = density.sum(axis=axis_idx) * dz
        else:
            weighted_sum = (density * field).sum(axis=axis_idx)
            density_sum = density.sum(axis=axis_idx)
            with np.errstate(divide="ignore", invalid="ignore"):
                proj_2d = np.where(density_sum > 0, weighted_sum / density_sum, np.nan)

        plot_data = proj_2d.T
        if args.field == "density":
            title = f"{{label}}, {args.axis}-projection"
        else:
            title = f"{{label}}, {args.axis}-projection (density-weighted avg)"

    label = args.field
    if args.log and args.field == "density":
        with np.errstate(divide="ignore"):
            plot_data = np.where(plot_data > 0, np.log10(plot_data), np.nan)
        label = "log10(density)" if args.mode == "slice" else "log10(column density)"
    elif args.mode == "projection" and args.field == "density":
        label = "column density"

    cmap = args.cmap
    if cmap is None and args.field in ("vx", "vy", "vz"):
        cmap = "coolwarm"

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(plot_data, origin="lower", extent=extent, aspect="equal", cmap=cmap,
                    vmin=args.vmin, vmax=args.vmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title.format(label=label))
    fig.colorbar(im, ax=ax, label=label)
    fig.tight_layout()

    if args.out:
        fig.savefig(args.out, dpi=150)
        print(f"plot_cartesian_map.py: wrote {args.out}")
    else:
        plt.show()

    return 0


if __name__ == "__main__":
    sys.exit(main())
