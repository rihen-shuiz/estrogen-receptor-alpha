"""Tracker Y15: multi-method outlier flags on the 788 x 788 RMSD matrix. Every chain is kept and labelled.

Inputs: results/rmsd_entry.npy, results/ext/zoo_mds10.npy (10-D classical MDS of the same matrix, Y04), S08 clusters,
results/ext/cluster_displacement_profile.csv (complete-H12 medoids, Y10).
Methods and flag rules (fixed before looking at which chains get flagged):
  knn      mean RMSD to the 5 nearest chains; flag robust z > 3.5 (z = (x - median) / (1.4826 MAD))
  medoid   RMSD to the nearer of the two complete-H12 main-cluster medoids (cl7 5TM1|A|1, cl4 7UJ7|C|1); robust z > 3.5
  lof      local outlier factor (sklearn, precomputed RMSD, 20 neighbours); flag LOF > 1.5
  iforest  isolation forest on the 10-D MDS vectors (300 trees, contamination 'auto', seed 42); flag = predicted -1
  mahal    robust Mahalanobis distance (MinCovDet) on the 5-D MDS vectors within the chain's main cluster (minor
           clusters: the nearer main cluster); flag chi-square(5) p < 0.001
n_methods = number of flags (0-5); consensus outlier = n_methods >= 3.
Outputs: results/ext/outlier_flags.csv, figures/ext/Y15_outliers.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_outliers.py
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2
from sklearn.covariance import MinCovDet
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (GREY, INK, INK2, RED, RES, SEED, VIOLET, YELLOW, chain_table, mpl_setup, rmsd, save_fig,
                        scatter_clusters, write_numbers)


def robust_z(x):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    return (x - med) / mad


def main():
    d, idx = rmsd()
    t = chain_table()
    x10 = np.load(RES / "zoo_mds10.npy")
    out = pd.DataFrame({"chain_key": t.chain_key, "entry_id": t.entry_id, "cluster": t.cluster})
    dd = d.copy()
    np.fill_diagonal(dd, np.inf)
    knn = np.sort(dd, 1)[:, :5].mean(1)
    out["knn5_mean_A"] = knn
    out["knn_z"] = robust_z(knn)
    out["knn_flag"] = out.knn_z > 3.5
    prof = pd.read_csv(RES / "cluster_displacement_profile.csv")
    m7 = int(np.where(idx.chain_key == prof.medoid_cl7.iloc[0])[0][0])
    m4 = int(np.where(idx.chain_key == prof.medoid_cl4.iloc[0])[0][0])
    md = np.minimum(d[:, m7], d[:, m4])
    out["medoid_dist_A"] = md
    out["medoid_z"] = robust_z(md)
    out["medoid_flag"] = out.medoid_z > 3.5
    lof = LocalOutlierFactor(n_neighbors=20, metric="precomputed").fit(d)
    out["lof"] = -lof.negative_outlier_factor_
    out["lof_flag"] = out.lof > 1.5
    iso = IsolationForest(n_estimators=300, random_state=SEED).fit(x10)
    out["iforest_score"] = -iso.score_samples(x10)
    out["iforest_flag"] = iso.predict(x10) == -1
    x5 = x10[:, :5]
    home = np.where(np.isin(t.cluster, (7, 4)), t.cluster, np.where(d[:, m7] < d[:, m4], 7, 4))
    out["home_cluster"] = home
    mah = np.full(len(t), np.nan)
    for c in (7, 4):
        m = (t.cluster.values == c)
        mcd = MinCovDet(random_state=SEED).fit(x5[m])
        sel = home == c
        mah[sel] = mcd.mahalanobis(x5[sel])
    out["mahal_d2"] = mah
    out["mahal_p"] = chi2.sf(mah, 5)
    out["mahal_flag"] = out.mahal_p < 0.001
    flags = ["knn_flag", "medoid_flag", "lof_flag", "iforest_flag", "mahal_flag"]
    out["n_methods"] = out[flags].sum(1)
    out["consensus_outlier"] = out.n_methods >= 3
    out.to_csv(RES / "outlier_flags.csv", index=False)

    emb = pd.read_csv(RES / "embedding_rmsd.csv")
    plt = mpl_setup()
    fig, axs = plt.subplots(2, 4, figsize=(19, 9))
    axs = axs.ravel()
    names = {"knn_flag": "kNN distance (k=5), robust z > 3.5", "medoid_flag": "distance to main medoid, robust z > 3.5",
             "lof_flag": "LOF (20 nb) > 1.5", "iforest_flag": "isolation forest (MDS-10)",
             "mahal_flag": "robust Mahalanobis p < 0.001"}
    for ax, f in zip(axs, flags):
        m = out[f].values
        ax.scatter(emb.mds1[~m], emb.mds2[~m], s=4, c=GREY, lw=0)
        ax.scatter(emb.mds1[m], emb.mds2[m], s=14, c=RED, lw=0.3, edgecolor="white")
        ax.set_title(f"{names[f]}: {m.sum()} flagged", fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    ax = axs[5]
    scatter_clusters(ax, emb.mds1, emb.mds2, t.cluster, s=5, legend=False)
    m = out.consensus_outlier.values
    ax.scatter(emb.mds1[m], emb.mds2[m], s=60, facecolor="none", edgecolor=INK, lw=1.2,
               label=f"consensus (>= 3 methods): {m.sum()}")
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title("Consensus outliers on the shared MDS map", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
    ax = axs[6]
    vc = out.n_methods.value_counts().sort_index()
    ax.bar(vc.index, vc.values, color=[GREY if i < 3 else RED for i in vc.index])
    for i, v in vc.items():
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set(xlabel="number of methods flagging the chain", ylabel="chains (log)", title="How many methods agree?")
    ax = axs[7]
    jac = np.zeros((5, 5))
    for i, a in enumerate(flags):
        for j, b in enumerate(flags):
            A, B = out[a].values, out[b].values
            jac[i, j] = (A & B).sum() / max((A | B).sum(), 1)
    im = ax.imshow(jac, cmap="esr_seq", vmin=0, vmax=1)
    lab = [f.replace("_flag", "") for f in flags]
    ax.set_xticks(range(5), lab, rotation=40); ax.set_yticks(range(5), lab)
    for i in range(5):
        for j in range(5):
            ax.text(j, i, f"{jac[i, j]:.2f}", ha="center", va="center", fontsize=8, color="white" if jac[i, j] > .6 else INK)
    ax.set_title("Jaccard overlap of flagged sets", fontsize=9)
    cons = out[out.consensus_outlier]
    fig.suptitle(f"Five outlier detectors flag {int((out.n_methods > 0).sum())} chains at least once; "
                 f"{len(cons)} chains ({cons.entry_id.nunique()} entries) are flagged by >= 3 methods",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y15_outliers", section="D. Outliers and quality", task="Y15",
             shows="flags of five outlier methods on the shared MDS map, consensus, agreement counts and Jaccard overlap",
             takeaway=(f"{len(cons)} consensus outliers in {cons.entry_id.nunique()} entries; "
                       f"{int(cons.cluster.isin((7, 4)).sum())} of them sit inside a main cluster, the rest are the minor clusters."),
             source="results/ext/outlier_flags.csv", script="scripts/ext_outliers.py")
    write_numbers("Y15", dict(n_any=int((out.n_methods > 0).sum()), n_consensus=len(cons),
                              n_consensus_entries=int(cons.entry_id.nunique()),
                              per_method={f: int(out[f].sum()) for f in flags},
                              consensus_chains=cons.chain_key.tolist()))
    print({f: int(out[f].sum()) for f in flags})
    print(cons[["chain_key", "cluster", "n_methods", "knn5_mean_A", "lof"]].to_string())


if __name__ == "__main__":
    main()
