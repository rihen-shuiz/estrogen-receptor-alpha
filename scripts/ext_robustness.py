"""Tracker Y28: robustness grid - clustering method x alignment core / position set x quality subset.

Matrices ("core" axis): core v3 (S04 results/rmsd_entry.npy), old core (K09 results/rmsd_old_core.npy), and the Y42
position-restricted matrices on the core-v3 superposition: fixed positions observed in >= 99% of chains
(results/ext/rmsd_fixedpos.npy) and no H11-H12 loop / H12 / tail (results/ext/rmsd_noh12.npy) - the last two are
negative controls (they remove the region that defines the split).
Methods: average linkage cut at 1.5 A (the S08 method), average linkage K=3, spectral K=2 (sigma 2 A), k-means K=2
and Ward K=2 on the 10-D classical MDS of the matrix, GMM K=2 on the same MDS.
Subsets (every subset is a slice of the full matrix, never a replacement): all 788 chains, tier A, R-free <= 0.25,
one chain per entry (first matrix chain), complete H12 only.
Metrics on the chains of the subset: ARI against the S08 1.5 A labels; main-split agreement = after mapping each
method cluster to the S08 main cluster (cl7 / cl4) holding most of its main-cluster chains, the fraction of the subset's
cl7/cl4 chains placed on the right side.
Outputs: results/ext/robustness_grid.csv, figures/ext/Y28_robustness_grid.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_robustness.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import KMeans, SpectralClustering
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import INK, RES, ROOT, SEED, chain_table, mpl_setup, rmsd, save_fig, write_numbers
from ext_zoo import classical_mds

warnings.filterwarnings("ignore")


def run(method, d):
    if method == "average cut 1.5 A (S08)":
        return fcluster(linkage(squareform(d, checks=False), "average"), 1.5, "distance")
    if method == "average K=3":
        return fcluster(linkage(squareform(d, checks=False), "average"), 3, "maxclust")
    if method == "spectral K=2":
        return SpectralClustering(2, affinity="precomputed", random_state=SEED,
                                  assign_labels="cluster_qr").fit_predict(np.exp(-d ** 2 / 8.0))
    x = classical_mds(d, 10)
    if method == "k-means K=2":
        return KMeans(2, n_init=20, random_state=SEED).fit_predict(x)
    if method == "Ward K=2":
        return fcluster(linkage(x, "ward"), 2, "maxclust")
    if method == "GMM K=2":
        return GaussianMixture(2, random_state=SEED, n_init=3, reg_covar=1e-4).fit_predict(x)
    raise ValueError(method)


def split_agreement(lab, ref):
    main = np.isin(ref, (7, 4))
    ok = 0
    for c in np.unique(lab):
        m = (lab == c) & main
        if m.any():
            maj = 7 if (ref[m] == 7).sum() >= (ref[m] == 4).sum() else 4
            ok += (ref[m] == maj).sum()
    return ok / main.sum()


def main():
    d, idx = rmsd()
    t = chain_table()
    mats = {"core v3": d, "old core": np.load(ROOT / "results/rmsd_old_core.npy"),
            "fixed positions (control)": np.load(RES / "rmsd_fixedpos.npy"),
            "no H12 528-550 (control)": np.load(RES / "rmsd_noh12.npy")}
    subsets = {"all chains": np.ones(len(t), bool), "tier A": (t.tier == "A").values,
               "R-free <= 0.25": (t.r_free <= 0.25).values, "one chain per entry": ~t.entry_id.duplicated().values,
               "complete H12": (t.h12_label == "complete").values}
    methods = ["average cut 1.5 A (S08)", "average K=3", "spectral K=2", "k-means K=2", "Ward K=2", "GMM K=2"]
    ref = t.cluster.values
    rows = []
    for cname, D in mats.items():
        for sname, m in subsets.items():
            Dm = D[np.ix_(m, m)]
            for meth in methods:
                lab = run(meth, Dm)
                rows.append(dict(method=meth, core=cname, subset=sname, n_chains=int(m.sum()),
                                 n_clusters=len(set(lab)), ari_vs_ref=adjusted_rand_score(ref[m], lab),
                                 main_split_agreement=split_agreement(lab, ref[m])))
    g = pd.DataFrame(rows)
    g.to_csv(RES / "robustness_grid.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 2, figsize=(19, 9))
    for ax, val, ttl in ((axs[0], "main_split_agreement", "cl7/cl4 split agreement"), (axs[1], "ari_vs_ref", "ARI vs S08 1.5 A")):
        piv = g.pivot_table(index=["subset", "method"], columns="core", values=val, sort=False)[list(mats)]
        im = ax.imshow(piv.values, cmap="esr_seq", vmin=0 if val == "ari_vs_ref" else 0.5, vmax=1, aspect="auto")
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                v = piv.values[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7, color="white" if v > 0.8 else INK)
        ax.set_xticks(range(piv.shape[1]), piv.columns, rotation=15)
        ax.set_yticks(range(piv.shape[0]), [f"{s} | {m}" for s, m in piv.index], fontsize=7)
        for k in range(1, len(subsets)):
            ax.axhline(k * len(methods) - 0.5, color="white", lw=2)
        ax.axvline(1.5, color="white", lw=3)
        fig.colorbar(im, ax=ax, fraction=0.03, label=ttl)
        ax.set_title(ttl)
    real = g[g.core.isin(["core v3", "old core"])]
    ctrl = g[~g.core.isin(["core v3", "old core"])]
    fig.suptitle(f"Robustness: with either real core the split is recovered by every method in every subset (agreement "
                 f"{real.main_split_agreement.min():.2f}-{real.main_split_agreement.max():.2f}); remove H12 and it falls to "
                 f"{ctrl.main_split_agreement.min():.2f}-{ctrl.main_split_agreement.max():.2f}", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y28_robustness_grid", section="M. Robustness", task="Y28",
             shows="heatmaps of main-split agreement and ARI vs S08 for 6 methods x 4 matrices (2 cores + 2 negative controls) x 5 subsets",
             takeaway=(f"Main-split agreement {real.main_split_agreement.min():.2f}-{real.main_split_agreement.max():.2f} for every "
                       f"method/subset on core v3 and the old core; the H12-free controls drop to "
                       f"{ctrl.main_split_agreement.min():.2f}-{ctrl.main_split_agreement.max():.2f}."),
             source="results/ext/robustness_grid.csv", script="scripts/ext_robustness.py")
    write_numbers("Y28", dict(n_cells=len(g), real_min=float(real.main_split_agreement.min()),
                              real_max=float(real.main_split_agreement.max()), ctrl_min=float(ctrl.main_split_agreement.min()),
                              ctrl_max=float(ctrl.main_split_agreement.max()), real_ari_min=float(real.ari_vs_ref.min())))
    print(g.groupby("core")[["main_split_agreement", "ari_vs_ref"]].agg(["min", "max"]))


if __name__ == "__main__":
    main()
