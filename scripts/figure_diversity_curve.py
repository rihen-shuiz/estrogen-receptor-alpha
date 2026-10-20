"""Tracker R03: cumulative diversity figure, mean +/- 1 SD band per cut and a dashed line at 90% of total diversity.

Input: results/diversity_curve.csv (R02) and results/saturation.csv (R04, the n_90 markers).
The y axis is the share of the cut's cluster total seen, so all four cuts share one 90% line.
Output: figures/diversity_curve.png. Light surface; categorical slots 1-4 in cut order.
Table twins: results/diversity_curve.csv, results/saturation.csv.
Run from the repo root: python scripts/figure_diversity_curve.py
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/diversity_curve.png"
INK, INK_SECONDARY, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]


def main():
    curve = pd.read_csv(ROOT / "results/diversity_curve.csv", dtype={"cut": str})
    summary = pd.read_csv(ROOT / "results/saturation.csv", dtype={"cut": str})
    n_chains = int(curve.n.max())
    fig, ax = plt.subplots(figsize=(8, 4.8), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for colour, row in zip(SLOTS, summary.itertuples(index=False)):
        part = curve[curve.cut == row.cut]
        mean = part.mean_clusters / row.n_clusters_total
        sd = part.sd_clusters / row.n_clusters_total
        ax.fill_between(part.n, (mean - sd).clip(lower=0), (mean + sd).clip(upper=1), color=colour, alpha=0.15, lw=0)
        ax.plot(part.n, mean, color=colour, lw=1.4,
                label=f"{row.cut} Å ({row.n_clusters_total} clusters): 90% at N={row.n_90}")
        ax.plot([row.n_90], [0.9], marker="o", ms=6, color=colour, mec=SURFACE, mew=1, ls="none")
    ax.axhline(0.9, color=INK_SECONDARY, lw=1.0, ls="--")
    ax.text(n_chains, 0.9, " 90%", color=INK_SECONDARY, fontsize=8, va="center")
    ax.set_xlim(0, n_chains)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel(f"Chains added (N, of {n_chains})", color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel("Share of total clusters seen (mean ± 1 SD)", color=INK_SECONDARY, fontsize=9)
    seed, orderings = int(summary.seed.iloc[0]), int(summary.n_orderings.iloc[0])
    ax.set_title(f"Cumulative structural diversity, {orderings} random orderings (seed {seed})",
                 loc="left", fontsize=11, color=INK, fontweight="bold")
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.legend(frameon=False, loc="lower right", fontsize=8, labelcolor=INK)
    fig.tight_layout()
    fig.savefig(OUT, facecolor=SURFACE)
    print(f"R03: {OUT.relative_to(ROOT)} written")


if __name__ == "__main__":
    main()
