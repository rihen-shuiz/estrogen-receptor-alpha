"""Tracker X04: ESR2 PCA landscape, same method as S11 / block 12, ESR2 only.

Input: results/rmsd_esr2.npy and results/rmsd_esr2_index.txt (X02), results/esr2_clusters.csv (X03).
Method: esr_analysis.structural_pca.structural_pca (sklearn PCA, 2 components, unscaled RMSD rows). Signs are fixed so
that the largest X03 cluster has a negative mean PC1 and PC2 (display only).
Colour: X03 cluster at the primary cut, categorical slots in fixed order by cluster size. The ESR1 S11 figure colours
by ligand class; ESR2 ligand classes are not labelled in this repository, so ESR2 uses structural clusters.
Outputs: figures/esr2_pca_landscape.png and results/esr2_pca_coords.csv (chain_key, entry_id, pc1, pc2, cluster) with
the explained variance in results/esr2_pca_variance.json.
Run from the repo root: PYTHONPATH=. python scripts/esr2_pca.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from esr_analysis.structural_pca import structural_pca  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = f'cut_{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE, INK, INK_SECONDARY, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def main():
    keys = (ROOT / "results/rmsd_esr2_index.txt").read_text().split()
    clusters = pd.read_csv(ROOT / "results/esr2_clusters.csv").set_index("chain_key").loc[keys]
    coords, variance = structural_pca(np.load(ROOT / "results/rmsd_esr2.npy"))
    order = clusters[CUT].value_counts().index.tolist()
    largest = (clusters[CUT] == order[0]).to_numpy()
    for k in range(2):
        if coords[largest, k].mean() > 0:
            coords[:, k] *= -1
    table = pd.DataFrame({"chain_key": keys, "entry_id": clusters.entry_id.to_numpy(), "pc1": coords[:, 0],
                          "pc2": coords[:, 1], "cluster": clusters[CUT].to_numpy()})
    table.to_csv(ROOT / "results/esr2_pca_coords.csv", index=False)
    (ROOT / "results/esr2_pca_variance.json").write_text(json.dumps(
        {"explained_variance_ratio": variance.tolist(), "features": "unscaled RMSD rows", "n_chains": len(keys)}, indent=2) + "\n")

    fig, ax = plt.subplots(figsize=(7, 5.2), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for slot, cluster in enumerate(order):
        part = table[table.cluster == cluster]
        ax.scatter(part.pc1, part.pc2, s=40, color=SLOTS[slot], edgecolor=SURFACE, linewidth=1.5,
                   label=f"cluster {cluster} ({len(part)} chains, {part.entry_id.nunique()} entries)")
    ax.set_xlabel(f"PC1 ({variance[0]:.1%} of variance)", color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel(f"PC2 ({variance[1]:.1%} of variance)", color=INK_SECONDARY, fontsize=9)
    ax.set_title(f"ESR2 conformational landscape, {len(keys)} chains (PCA of RMSD rows)", loc="left",
                 fontsize=11, color=INK, fontweight="bold")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, title=f"X03 cluster at {CUT[4:]} Å", title_fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/esr2_pca_landscape.png", facecolor=SURFACE)
    print(f"X04: {len(keys)} chains; variance {variance.round(4).tolist()}; clusters {order}")


if __name__ == "__main__":
    main()
