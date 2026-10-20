"""Tracker X12: PCA of P2Rank top-pocket descriptors across all eligible chains.

Input: results/pockets.csv (X09, one row per chain), results/entry_clusters.csv (S08 structural cluster at the
primary cut). Descriptors: pocket_score, pocket_probability, sas_points, surf_atoms, n_residues, q18_residue_jaccard,
standardised (z-scores) before sklearn PCA (2 components).
Colour: S08 structural cluster, categorical slots in fixed order by cluster size. The task asks whether pockets
cluster by ligand class; ligand-class labels (Q10) do not exist yet, so ligand class cannot be shown.
Outputs: figures/pocket_pca.png, results/pocket_pca_coords.csv (chain_key, pc1, pc2, cluster) and
results/pocket_pca_loadings.csv (variance ratio and loadings per descriptor).
Run from the repo root: python scripts/pocket_pca.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.decomposition import PCA  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = f'cut_{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
DESCRIPTORS = ["pocket_score", "pocket_probability", "sas_points", "surf_atoms", "n_residues", "q18_residue_jaccard"]
SLOTS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SURFACE, INK, INK_SECONDARY, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"


def main():
    pockets = pd.read_csv(ROOT / "results/pockets.csv").merge(
        pd.read_csv(ROOT / "results/entry_clusters.csv")[["chain_key", CUT]], on="chain_key", validate="one_to_one")
    data = pockets.dropna(subset=DESCRIPTORS)
    pca = PCA(n_components=2)
    coords = pca.fit_transform(StandardScaler().fit_transform(data[DESCRIPTORS]))
    out = data[["entry_id", "chain_key", CUT]].rename(columns={CUT: "cluster"}).assign(pc1=coords[:, 0], pc2=coords[:, 1])
    out.to_csv(ROOT / "results/pocket_pca_coords.csv", index=False)
    loadings = pd.DataFrame(pca.components_.T, index=DESCRIPTORS, columns=["pc1", "pc2"])
    loadings.loc["explained_variance_ratio"] = pca.explained_variance_ratio_
    loadings.to_csv(ROOT / "results/pocket_pca_loadings.csv")

    order = out.cluster.value_counts().index.tolist()
    fig, ax = plt.subplots(figsize=(7.5, 5.4), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for slot, cluster in enumerate(order):
        part = out[out.cluster == cluster]
        ax.scatter(part.pc1, part.pc2, s=14 if len(part) > 50 else 34, color=SLOTS[slot % len(SLOTS)],
                   alpha=0.55 if len(part) > 50 else 0.95, edgecolor=SURFACE, linewidth=0.6,
                   label=f"cluster {cluster} ({len(part)} chains)")
    var = pca.explained_variance_ratio_
    ax.set_xlabel(f"PC1 ({var[0]:.1%} of variance)", color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel(f"PC2 ({var[1]:.1%} of variance)", color=INK_SECONDARY, fontsize=9)
    ax.set_title(f"P2Rank top-pocket descriptors, {len(out)} ESR1 chains (PCA, z-scored)", loc="left",
                 fontsize=11, color=INK, fontweight="bold")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, title=f"S08 structural cluster at {CUT[4:]} Å", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/pocket_pca.png", facecolor=SURFACE)
    print(f"X12: {len(out)} chains; variance {var.round(3).tolist()}\n{loadings.round(3).to_string()}")


if __name__ == "__main__":
    main()
