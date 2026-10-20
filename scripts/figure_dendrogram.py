"""Tracker S10: dendrogram of the all-eligible-chain average linkage (keep-all run), with the primary cut drawn.

Input: results/linkage_entry.npy (S08 linkage), results/entry_clusters.csv (S08, cluster count at the cut),
primary_structural_threshold_A from analysis_settings.json.
Output: figures/dendrogram.png. One series in categorical slot 1 (blue); the cut is a dashed reference line.
Leaves (chains) are not labeled; the chain labels are in results/entry_clusters.csv. Light surface only.
Run from the repo root: python scripts/figure_dendrogram.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.cluster.hierarchy import dendrogram  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
OUT = ROOT / "figures/dendrogram.png"

# Reference palette (docs palette.md), light mode.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
SERIES = "#2a78d6"  # categorical slot 1


def main():
    cut = json.loads(SETTINGS.read_text())["primary_structural_threshold_A"]
    linkage_matrix = np.load(ROOT / "results/linkage_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    n_chains = len(clusters)
    n_entries = clusters.entry_id.nunique()
    n_clusters = int(clusters[f"cut_{cut:.1f}"].nunique())

    fig, ax = plt.subplots(figsize=(10, 5.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    dendrogram(linkage_matrix, no_labels=True, link_color_func=lambda _: SERIES, ax=ax)
    ax.axhline(cut, color=INK_SECONDARY, linestyle="--", linewidth=1.2)
    ax.text(0.99, cut, f"cut {cut:.1f} Å: {n_clusters} clusters", transform=ax.get_yaxis_transform(),
            ha="right", va="bottom", color=INK_SECONDARY, fontsize=9)

    ax.set_title("Average-linkage dendrogram, all eligible ESR1 chains", loc="left", color=INK, fontsize=12, fontweight="bold")
    ax.set_ylabel("C-alpha RMSD linkage distance (Å)", color=INK_SECONDARY)
    ax.set_xlabel(f"{n_chains} chains from {n_entries} entries (leaves unlabeled)", color=INK_SECONDARY)
    ax.set_xticks([])
    ax.tick_params(axis="y", colors=MUTED)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(BASELINE)

    fig.savefig(OUT, facecolor=SURFACE, bbox_inches="tight")
    print(f"S10: {OUT.relative_to(ROOT)} written; {n_chains} leaves ({n_entries} entries); cut {cut} A gives {n_clusters} clusters")


if __name__ == "__main__":
    main()
