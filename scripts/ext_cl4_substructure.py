"""Tracker Y45: inside cl4 - what splits the antagonist-like cluster?

Input: the cl4 rows/columns of results/rmsd_entry.npy (329 chains), chain_table(), Y02 chain metadata, Y08 chemistry,
Y09 descriptors, Y11 structural features, Y44 H12 contacts, S03 coordinates.
Sub-clustering: the frozen S08 labels at the 1.0 A cut restricted to cl4 (sub-clusters of >= 10 chains kept, smaller
ones pooled as 'small'). A first attempt (average linkage on the cl4 sub-matrix, K = 2..5 by silhouette) picked a
degenerate 327 + 2 split and is only reported (silhouettes in the numbers file). Sub-clusters are described by H12 geometry, H12 completeness,
lattice contacts at H12, peptide, construct mutation, ligand basic amine / size / chemistry cluster and space group, and
by a per-residue distance profile between the two largest sub-clusters (mean C-alpha positions).
Outputs: results/ext/cl4_subclusters.csv, results/ext/cl4_subcluster_summary.csv, figures/ext/Y45_cl4_substructure.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_cl4_substructure.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.manifold import MDS
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (CATEGORICAL, GREY, INK, INK2, RES, SEED, VIOLET, YELLOW, chain_table, coords, helix_bands,
                        mpl_setup, rmsd, save_fig, write_numbers)

SUBCOL = [VIOLET, YELLOW, CATEGORICAL[5], CATEGORICAL[4], CATEGORICAL[7]]


def main():
    d, idx = rmsd()
    t = chain_table()
    m = (t.cluster == 4).values
    D = d[np.ix_(m, m)]
    Z = linkage(squareform(D, checks=False), "average")
    sil = {k: silhouette_score(D, fcluster(Z, k, "maxclust"), metric="precomputed") for k in range(2, 6)}
    c10 = t.loc[m, "cut_1.0"].values
    vc = pd.Series(c10).value_counts()
    big = [c for c in vc.index if vc[c] >= 10]
    remap = {c: i + 1 for i, c in enumerate(big)}
    lab = np.array([remap.get(c, len(big) + 1) for c in c10])
    K = int(lab.max())
    sil_used = silhouette_score(D, lab, metric="precomputed")
    sub = t[m][["chain_key", "entry_id", "ligand_id", "h12_label", "space_group", "cut_1.0"]].copy()
    sub["sub_cluster"] = lab
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status", "peptide_present"]]
    sf = pd.read_csv(RES / "structural_features.csv")[["chain_key", "h12_shift_A", "h12_h3_angle_deg", "k362_e542_A",
                                                        "h12_n_ca", "ligand_heavy_atoms", "ligand_basic_amine"]]
    hc = pd.read_csv(RES / "h12_crystal_contacts.csv")[["chain_key", "h12_contact_atoms"]]
    ch = pd.read_csv(RES / "chem_labels_best.csv")[["ligand_id", "chem_cluster_auto"]]
    sub = sub.merge(cm, on="chain_key").merge(sf, on="chain_key").merge(hc, on="chain_key").merge(ch, on="ligand_id", how="left")
    sub.to_csv(RES / "cl4_subclusters.csv", index=False)
    summ = sub.groupby("sub_cluster").agg(
        n_chains=("chain_key", "size"), n_entries=("entry_id", "nunique"), h12_shift_A=("h12_shift_A", "median"),
        h12_h3_angle_deg=("h12_h3_angle_deg", "median"), k362_e542_A=("k362_e542_A", "median"),
        h12_complete=("h12_label", lambda s: (s == "complete").mean()), h12_contacts=("h12_contact_atoms", "median"),
        L536S=("mutation_status", lambda s: (s == "L536S").mean()), Y537S=("mutation_status", lambda s: (s == "Y537S").mean()),
        basic_amine=("ligand_basic_amine", lambda s: s.fillna(False).astype(bool).mean()),
        ligand_heavy_atoms=("ligand_heavy_atoms", "median"),
        top_space_group=("space_group", lambda s: s.value_counts().index[0]),
        top_space_group_share=("space_group", lambda s: s.value_counts(normalize=True).iloc[0]),
        top_chem=("chem_cluster_auto", lambda s: s.value_counts().index[0] if s.notna().any() else np.nan))
    summ.to_csv(RES / "cl4_subcluster_summary.csv")
    x, pos = coords()
    xc = x[m]
    prof = np.linalg.norm(np.nanmean(xc[lab == 1], 0) - np.nanmean(xc[lab == 2], 0), axis=1)

    plt = mpl_setup()
    fig = plt.figure(figsize=(18, 9.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 0.8])
    ax = fig.add_subplot(gs[0, 0])
    z = MDS(2, dissimilarity="precomputed", random_state=SEED, n_init=4, normalized_stress="auto").fit_transform(D)
    for c in sorted(set(lab)):
        mm = lab == c
        nm = f"cl4-{c}" if c < K else "cl4-small"
        ax.scatter(z[mm, 0], z[mm, 1], s=12, c=SUBCOL[c - 1] if c < K else GREY, lw=0, label=f"{nm} (n={mm.sum()})")
    ax.legend(fontsize=8)
    ax.set(xlabel="MDS 1 (A)", ylabel="MDS 2 (A)", title=f"cl4 sub-clusters = S08 1.0 A cut\n({K - 1} groups of >= 10 chains + small; silhouette {sil_used:.2f})")
    ax = fig.add_subplot(gs[0, 1:])
    feats = ["h12_shift_A", "h12_h3_angle_deg", "k362_e542_A", "h12_complete", "h12_contacts", "L536S", "Y537S",
             "basic_amine", "ligand_heavy_atoms", "top_space_group_share"]
    S = summ[feats].astype(float)
    Zs = (S - S.mean()) / S.std().replace(0, 1)
    im = ax.imshow(Zs.values, cmap="esr_div", vmin=-1.5, vmax=1.5, aspect="auto")
    for i in range(S.shape[0]):
        for j in range(S.shape[1]):
            v = S.values[i, j]
            txt = f"{v:.0%}" if feats[j] in ("h12_complete", "L536S", "Y537S", "basic_amine", "top_space_group_share") else f"{v:.1f}"
            if feats[j] == "top_space_group_share":
                txt = f"{summ.top_space_group.iloc[i]} {v:.0%}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=8, color=INK)
    ax.set_xticks(range(len(feats)), feats, rotation=25, ha="right")
    ax.set_yticks(range(S.shape[0]), [f"cl4-{c} (n={n})" for c, n in zip(summ.index, summ.n_chains)])
    fig.colorbar(im, ax=ax, fraction=0.025, label="z-score across sub-clusters (cells show values)")
    ax.set_title("What differs between the cl4 sub-clusters?")
    ax = fig.add_subplot(gs[1, :])
    helix_bands(ax)
    ax.plot(pos, prof, color=VIOLET, lw=1.5)
    ax.set_yscale("symlog", linthresh=1)
    top = pd.Series(prof, index=pos).sort_values(ascending=False).head(3)
    ax.set(xlabel="ESR1 residue", ylabel="cl4-1 vs cl4-2 mean C-alpha distance (A)",
           title=f"Where cl4-1 and cl4-2 differ: largest at {', '.join(map(str, top.index))}")
    a1, a2 = summ.loc[1], summ.loc[2]
    fig.suptitle(f"Inside the antagonist-like cluster: cl4-1 (n={int(a1.n_chains)}, L536S {a1.L536S:.0%}, {a1.top_space_group}) vs "
                 f"cl4-2 (n={int(a2.n_chains)}, L536S {a2.L536S:.0%}, {a2.top_space_group}) differ in the H11-H12 loop "
                 f"({', '.join(map(str, top.index[:2]))}) and the H2-H3 loop, not in H12", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y45_cl4_substructure", section="C. Why two clusters", task="Y45",
             shows="MDS of the cl4 sub-matrix coloured by sub-cluster; feature heatmap per sub-cluster; per-residue profile between the two largest sub-clusters",
             takeaway=(f"At 1.0 A cl4 splits into {a1.n_chains:.0f} + {a2.n_chains:.0f} chains (+ small groups); see the heatmap for "
                       f"which variables separate them."),
             source="results/ext/cl4_subclusters.csv, results/ext/cl4_subcluster_summary.csv", script="scripts/ext_cl4_substructure.py")
    write_numbers("Y45", dict(K=int(K), silhouette_used=float(sil_used), silhouettes_avg_maxclust={int(k): float(v) for k, v in sil.items()},
                              summary=summ.reset_index().round(3).astype(str).to_dict("records"),
                              top_positions=top.index.tolist()))
    print(sil, summ.T)


if __name__ == "__main__":
    main()
