"""Tracker S14: threshold sensitivity of the selected-set clustering, cut from 1.0 to 2.5 A.

Input: results/linkage_entry.npy (S08 linkage), results/cluster_stats.csv (S09, the four reported cuts),
primary_structural_threshold_A from analysis_settings.json.
Method: flat clusters from the saved linkage every 0.05 A (scipy fcluster, distance criterion, same as S08).
Steps show the value at each sampled cut. The values at the four reported cuts are checked against
results/cluster_stats.csv before drawing.
Output: figures/structural_sensitivity.png. Two panels, one series each (no dual axis): cluster count and
share of entries in the three largest clusters. Reported cuts are marked and labeled. Light surface only.
Run from the repo root: python scripts/figure_sensitivity.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.cluster.hierarchy import fcluster  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
OUT = ROOT / "figures/structural_sensitivity.png"

# Reference palette (docs palette.md), light mode.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SERIES = "#2a78d6"  # categorical slot 1, one series per panel
REPORTED = [1.0, 1.5, 2.0, 2.5]


def stats_at(linkage_matrix, cut, n_entries):
    labels = fcluster(linkage_matrix, t=cut, criterion="distance")
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    return int(sizes.size), float(sizes.head(3).sum() / n_entries)


def main():
    primary = json.loads(SETTINGS.read_text())["primary_structural_threshold_A"]
    linkage_matrix = np.load(ROOT / "results/linkage_entry.npy")
    n_entries = linkage_matrix.shape[0] + 1
    reported = pd.read_csv(ROOT / "results/cluster_stats.csv")

    grid = np.round(np.arange(1.0, 2.5 + 1e-9, 0.05), 2)
    counts, shares = np.array([stats_at(linkage_matrix, cut, n_entries) for cut in grid]).T
    for cut in REPORTED:
        row = reported[np.isclose(reported.cut, cut)].iloc[0]
        at = int(np.flatnonzero(np.isclose(grid, cut))[0])
        if counts[at] != row.n_clusters or not np.isclose(shares[at], row.top3_fraction):
            raise ValueError(f"Sweep at {cut} A disagrees with results/cluster_stats.csv")

    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.4), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    panels = [
        (left, counts, "Clusters", (0, 20), range(0, 21, 4), "{:d}", True),
        (right, shares * 100, "Entries in the 3 largest clusters (%)", (80, 100), range(80, 101, 5), "{:.0f}%", False),
    ]
    for ax, values, title, ylim, yticks, fmt, label_above in panels:
        ax.set_facecolor(SURFACE)
        ax.plot(grid, values, color=SERIES, linewidth=2, drawstyle="steps-post", solid_joinstyle="round",
                solid_capstyle="round")
        for cut in REPORTED:
            at = int(np.flatnonzero(np.isclose(grid, cut))[0])
            ax.plot(cut, values[at], "o", color=SERIES, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.5)
            text = fmt.format(int(values[at])) if label_above else f"{values[at]:.1f}%"
            ax.annotate(text, (cut, values[at]), textcoords="offset points", xytext=(0, 9 if label_above else -14),
                        ha="center", color=INK_SECONDARY, fontsize=9,
                        bbox=dict(boxstyle="square,pad=0.15", facecolor=SURFACE, edgecolor="none"))
        ax.axvline(primary, color=INK_SECONDARY, linestyle="--", linewidth=1)
        ax.set_xlim(0.95, 2.55)
        ax.set_ylim(*ylim)
        ax.set_yticks(list(yticks))
        ax.set_xticks([1.0, 1.5, 2.0, 2.5])
        ax.set_title(title, loc="left", color=INK, fontsize=11, fontweight="bold")
        ax.set_xlabel("Distance cut (Å)", color=INK_SECONDARY)
        ax.tick_params(colors=MUTED)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(BASELINE)

    left.annotate(f"primary {primary:.1f} Å", (primary, 20), textcoords="offset points", xytext=(4, -2),
                  ha="left", va="top", color=INK_SECONDARY, fontsize=9)
    fig.suptitle(f"Threshold sensitivity, selected ESR1 set ({n_entries} entries, one chain each)",
                 x=0.02, ha="left", color=INK, fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT, facecolor=SURFACE, bbox_inches="tight")
    print(f"S14: {OUT.relative_to(ROOT)} written; counts {int(counts[0])}->{int(counts[-1])} over 1.0-2.5 A")


if __name__ == "__main__":
    main()
