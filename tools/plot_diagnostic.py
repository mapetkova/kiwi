#!/usr/bin/env python3
"""Diagnostic scatter plot for a kiwi point-set stage checkpoint.

Usage:
    plot_diagnostic.py <points_file> <n_real> <out_png> <title>

<points_file> is any file written by kiwi's save_point_set() -- one
"id x y z" row per line (points-buff.dat or points-final.dat). Points with
id <= n_real are drawn as "real" (SPH-particle-derived) sites; points with
id > n_real are drawn as "buffer" sites, matching the id convention every
buffer strategy in point_set.cpp uses (buffer ids start at n_real + 1, see
add_graded_shells_buffer / add_coverage_grid_buffer).

This is invoked automatically by kiwi (via std::system, from diagnostics.cpp)
whenever "diagnostics": {"enabled": true} is set in a run's config -- it is
not part of kiwi's data contract, just a best-effort visual sanity check. Safe
to run by hand too, e.g. to re-plot an existing points-final.dat with a
different n_real.

Requires matplotlib (pip3 install matplotlib). Exits 3 (distinct from the
usage-error exit 2) if matplotlib isn't importable, so the caller can print a
specific "install matplotlib" hint.
"""
import sys


def main():
    if len(sys.argv) != 5:
        print("usage: plot_diagnostic.py <points_file> <n_real> <out_png> <title>",
              file=sys.stderr)
        return 2

    points_file, n_real_s, out_png, title = sys.argv[1:5]
    n_real = int(n_real_s)

    try:
        import matplotlib
        matplotlib.use("Agg")  # headless -- never try to open a GUI window
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"plot_diagnostic.py: matplotlib is required (pip3 install matplotlib): {e}",
              file=sys.stderr)
        return 3

    real_x, real_y = [], []
    buffer_x, buffer_y = [], []
    with open(points_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            pid = int(float(parts[0]))
            x, y = float(parts[1]), float(parts[2])
            if pid <= n_real:
                real_x.append(x)
                real_y.append(y)
            else:
                buffer_x.append(x)
                buffer_y.append(y)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(real_x, real_y, s=1, c="0.35", label=f"real ({len(real_x):,})")
    if buffer_x:
        ax.scatter(buffer_x, buffer_y, s=12, c="crimson",
                   label=f"buffer ({len(buffer_x):,})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="datalim")
    ax.legend(markerscale=4, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
