"""Cumulative diversity curve and 90% saturation N for each cut (tracker task P02).

Reads the cut labels written by P01 (results/provisional/clusters.csv). For each cut, 100 random
orderings of the 351 entries (one seeded generator, seed 42) give the number of distinct clusters seen
after N entries. Writes mean and SD for N = 1..351 (results/provisional/diversity_curve.csv), and per cut
the total cluster count and n_90: the smallest N whose mean reaches 90% of the total
(results/provisional/saturation.csv). Writes figures/provisional_diversity_curve.png with the 90% line.

Run from anywhere: python scripts/provisional_saturation.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provisional_common import CUTS, N_ORDERINGS, PROVISIONAL, ROOT, SEED, saturation  # noqa: E402

FIGURE = ROOT / "figures/provisional_diversity_curve.png"
INK, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # categorical slots 1-4, light mode, cut order


def main():
    clusters = pd.read_csv(PROVISIONAL / "clusters.csv", dtype=str, keep_default_na=False)
    curve_rows, summary_rows = [], []
    for cut in CUTS:
        label = f"{cut:.1f}"
        labels = clusters[f"cut_{label}"].astype(int).to_numpy()
        mean, sd, total, n_90 = saturation(labels)
        for position in range(len(labels)):
            curve_rows.append({"cut": label, "n": position + 1, "mean_clusters": mean[position],
                               "sd_clusters": sd[position], "n_clusters_total": total,
                               "seed": SEED, "n_orderings": N_ORDERINGS})
        summary_rows.append({"cut": label, "n_90": n_90, "n_clusters_total": total,
                             "threshold_90pct": 0.9 * total, "mean_at_n90": mean[n_90 - 1],
                             "seed": SEED, "n_orderings": N_ORDERINGS, "n_entries": len(labels)})
    curves, summary = pd.DataFrame(curve_rows), pd.DataFrame(summary_rows)
    curves.to_csv(PROVISIONAL / "diversity_curve.csv", index=False)
    summary.to_csv(PROVISIONAL / "saturation.csv", index=False)
    figure(curves, summary)
    print(summary.to_string(index=False))
    print(f"wrote {PROVISIONAL / 'diversity_curve.csv'}, {PROVISIONAL / 'saturation.csv'} and {FIGURE}")


def figure(curves, summary):
    fig, ax = plt.subplots(figsize=(7.0, 4.3), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for colour, row in zip(SLOTS, summary.itertuples(index=False)):
        part = curves[curves.cut == row.cut]
        share = part.mean_clusters.astype(float).to_numpy() / row.n_clusters_total
        ax.plot(part.n.astype(int), share, color=colour, lw=1.2,
                label=f"{row.cut} Å: 90% at N={row.n_90}")
        ax.plot([row.n_90], [0.9], marker="o", ms=6, color=colour, mec=SURFACE, mew=1, ls="none")
        ax.annotate(f"{row.cut} Å", (100, share[99]), textcoords="offset points", xytext=(0, 5),
                    ha="center", color=INK, fontsize=8)
    ax.axhline(0.9, color=MUTED, lw=0.8)
    ax.text(351, 0.9, " 90%", color=MUTED, fontsize=8, va="center")
    ax.set_xlim(0, 351)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("Entries added (N)", color=MUTED, fontsize=9)
    ax.set_ylabel("Mean share of total clusters seen", color=MUTED, fontsize=9)
    ax.set_title(f"Cluster discovery over {N_ORDERINGS} random orderings (seed {SEED})", loc="left",
                 fontsize=10, color=INK)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.legend(frameon=False, loc="lower right", fontsize=8, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(FIGURE, facecolor=SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    main()
