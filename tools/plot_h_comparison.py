#!/usr/bin/env python3
"""Comparison plot for tools/h_from_neighbours.cpp's output: the smoothing
length recovered by the Phantom-style neighbour-density iteration vs. the
smoothing length already stored in the input snapshot.

Usage:
    plot_h_comparison.py <h_compare.dat> <out_png> [title]

<h_compare.dat> is a file written by h_from_neighbours -- one
"id x y z m h_snapshot rho_snapshot h_computed rho_computed n_neighbours
n_iterations converged rel_diff_h" row per line (# header comment ignored).

Produces a two-panel figure: left is h_computed vs h_snapshot on a log-log
scatter with the y=x line for reference; right is a histogram of the signed
relative difference (h_computed - h_snapshot) / h_snapshot. Unconverged
particles (converged == 0), if any, are drawn in a distinct color in the
scatter panel rather than silently mixed in.

Also prints the same summary statistics h_from_neighbours itself prints to
stderr (mean/median/max |relative diff|, convergence fraction), so this can
be re-run standalone against a saved .dat file without re-running the C++
tool.

Requires matplotlib (pip3 install matplotlib). Exits 3 (distinct from the
usage-error exit 2) if matplotlib isn't importable, matching
plot_diagnostic.py's convention.
"""
import sys


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: plot_h_comparison.py <h_compare.dat> <out_png> [title]", file=sys.stderr)
        return 2

    data_file, out_png = sys.argv[1], sys.argv[2]
    title = sys.argv[3] if len(sys.argv) == 4 else "h_from_neighbours vs snapshot h"

    try:
        import matplotlib
        matplotlib.use("Agg")  # headless -- never try to open a GUI window
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as e:
        print(f"plot_h_comparison.py: matplotlib/numpy are required (pip3 install matplotlib numpy): {e}",
              file=sys.stderr)
        return 3

    h_snap, h_comp, converged, rel_diff = [], [], [], []
    with open(data_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 13:
                continue
            h_snap.append(float(parts[5]))
            h_comp.append(float(parts[7]))
            converged.append(int(parts[11]))
            rel_diff.append(float(parts[12]))

    h_snap = np.array(h_snap)
    h_comp = np.array(h_comp)
    converged = np.array(converged)
    rel_diff = np.array(rel_diff)
    n = len(h_snap)
    if n == 0:
        print(f"plot_h_comparison.py: no data rows found in {data_file}", file=sys.stderr)
        return 2

    abs_rel = np.abs(rel_diff)
    n_converged = int(converged.sum())
    print(f"n particles: {n}")
    print(f"converged: {n_converged}/{n} ({100.0 * n_converged / n:.2f}%)")
    print(f"|rel diff| vs snapshot h: mean={abs_rel.mean():.3e} "
          f"median={np.median(abs_rel):.3e} max={abs_rel.max():.3e}")
    for thresh in (0.001, 0.005, 0.01, 0.05, 0.1):
        frac = np.mean(abs_rel > thresh) * 100.0
        print(f"  frac |rel diff| > {thresh}: {frac:.3f}%")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.5))

    ok = converged == 1
    bad = ~ok
    ax1.scatter(h_snap[ok], h_comp[ok], s=3, c="0.35", alpha=0.5, label=f"converged ({ok.sum():,})")
    if bad.any():
        ax1.scatter(h_snap[bad], h_comp[bad], s=10, c="crimson", label=f"not converged ({bad.sum():,})")
    lo = min(h_snap.min(), h_comp.min())
    hi = max(h_snap.max(), h_comp.max())
    ax1.plot([lo, hi], [lo, hi], c="steelblue", lw=1, ls="--", label="y = x")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("h (snapshot)")
    ax1.set_ylabel("h (neighbour-density iteration)")
    ax1.set_title("recovered vs. snapshot h")
    ax1.legend(markerscale=3, fontsize=8, loc="upper left")

    ax2.hist(rel_diff, bins=100, color="0.35")
    ax2.set_xlabel("(h_computed - h_snapshot) / h_snapshot")
    ax2.set_ylabel("count")
    ax2.set_title("signed relative difference")
    ax2.axvline(0.0, c="steelblue", lw=1, ls="--")

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main())
