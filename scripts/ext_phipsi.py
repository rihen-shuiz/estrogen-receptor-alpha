"""Tracker Y12: backbone dihedral (phi/psi) PCA of the 788 matrix chains.

Input: data/clean/<entry>_<chain>.pdb (S17 cleaned chains, author numbering), results/lbd_extract/residue_mapping.csv
(author residue -> ESR1 position), chain_table().
Dihedrals: Bio.PDB PPBuilder per chain; phi and psi per residue, mapped to ESR1 positions 305-550.
Features: (sin phi, cos phi, sin psi, cos psi) per position, for positions with both angles in >= 90% of chains
(circular encoding, so -179 and 179 deg are close). Missing values are imputed with the column mean (and counted).
PCA (sklearn, 10 components). Per-position loading = sum of squared loadings of its 4 features on dPC1.
Rigid-body superposition plays no role here (dihedrals are internal), so this is an alignment-free check of the split.
Outputs: results/ext/phipsi_pca.csv (chain_key, dpc1..dpc3, n_imputed), results/ext/phipsi_loadings.csv,
figures/ext/Y12_phipsi_pca.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_phipsi.py
"""
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (INK2, RES, ROOT, SEED, VIOLET, chain_table, helix_bands, mpl_setup, save_fig, scatter_clusters,
                        write_numbers)

warnings.filterwarnings("ignore")
POS = np.arange(305, 551)


def dihedrals(path, amap):
    s = PDBParser(QUIET=True).get_structure("x", path)
    out = {}
    for pp in PPBuilder().build_peptides(s[0]):
        for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            key = (res.id[1], res.id[2].strip())
            p = amap.get(key)
            if p is not None and phi is not None and psi is not None:
                out[p] = (phi, psi)
    return out


def main():
    t = chain_table()
    lc = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False)[["chain_key", "chain_id"]]
    t = t.merge(lc, on="chain_key", how="left")
    rm = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "auth_seq_id", "insertion_code", "esr1_position", "has_residue"])
    rm = rm[rm.chain_key.isin(t.chain_key) & rm.has_residue & rm.esr1_position.notna() & rm.auth_seq_id.notna()]
    rm["icode"] = rm.insertion_code.fillna("").astype(str).replace({"?": "", ".": "", "nan": ""})
    maps = {k: {(int(a), i): int(p) for a, i, p in zip(g.auth_seq_id, g.icode, g.esr1_position)}
            for k, g in rm.groupby("chain_key")}
    ang = np.full((len(t), len(POS), 2), np.nan)
    for i, r in t.iterrows():
        f = ROOT / f"data/clean/{r.entry_id}_{r.chain_id}.pdb"
        if not f.exists():
            continue
        for p, (phi, psi) in dihedrals(f, maps.get(r.chain_key, {})).items():
            if 305 <= p <= 550:
                ang[i, p - 305] = (phi, psi)
    obs = np.isfinite(ang[..., 0]).mean(0)
    keep = obs >= 0.90
    a = ang[:, keep]
    feats = np.concatenate([np.sin(a[..., :1]), np.cos(a[..., :1]), np.sin(a[..., 1:]), np.cos(a[..., 1:])], axis=2)
    xf = feats.reshape(len(t), -1)
    n_imp = np.isnan(xf).sum(1) // 4
    col_mean = np.nanmean(xf, 0)
    xf = np.where(np.isnan(xf), col_mean, xf)
    pca = PCA(10, random_state=SEED).fit(xf)
    z = pca.transform(xf)
    out = pd.DataFrame({"chain_key": t.chain_key, "dpc1": z[:, 0], "dpc2": z[:, 1], "dpc3": z[:, 2],
                        "n_imputed_positions": n_imp, "cluster": t.cluster})
    out.to_csv(RES / "phipsi_pca.csv", index=False)
    load = (pca.components_[0].reshape(-1, 4) ** 2).sum(1)
    load2 = (pca.components_[1].reshape(-1, 4) ** 2).sum(1)
    ld = pd.DataFrame({"position": POS[keep], "load1": load, "load2": load2, "obs_frac": obs[keep]})
    # per-position circular mean difference of psi between clusters
    for c in (7, 4):
        m = t.cluster.values == c
        for k, nm in ((0, "phi"), (1, "psi")):
            v = ang[m][:, keep, k]
            ld[f"{nm}_mean_cl{c}"] = np.degrees(np.arctan2(np.nanmean(np.sin(v), 0), np.nanmean(np.cos(v), 0)))
    for nm in ("phi", "psi"):
        dlt = ld[f"{nm}_mean_cl7"] - ld[f"{nm}_mean_cl4"]
        ld[f"d{nm}_deg"] = (dlt + 180) % 360 - 180
    ld.to_csv(RES / "phipsi_loadings.csv", index=False)
    main_ = t.cluster.isin((7, 4)).values
    sil = silhouette_score(z[main_, :2], t.cluster.values[main_])

    plt = mpl_setup()
    fig = plt.figure(figsize=(17, 6.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 2.1], height_ratios=[1, 1])
    ax = fig.add_subplot(gs[:, 0])
    scatter_clusters(ax, z[:, 0], z[:, 1], t.cluster, s=9)
    ev = pca.explained_variance_ratio_
    ax.set(xlabel=f"dihedral PC1 ({ev[0] * 100:.1f}%)", ylabel=f"dihedral PC2 ({ev[1] * 100:.1f}%)",
           title=f"phi/psi PCA ({keep.sum()} positions): silhouette of cl7/cl4 = {sil:.2f}")
    ax = fig.add_subplot(gs[0, 1])
    helix_bands(ax)
    ax.bar(ld.position, ld.load1, color=VIOLET, width=0.9)
    ax.set(ylabel="dPC1 loading (sum sq.)", xlim=(304, 551))
    top = ld.sort_values("load1", ascending=False).head(4)
    for r in top.itertuples():
        ax.annotate(str(r.position), (r.position, r.load1), xytext=(3, 2), textcoords="offset points", fontsize=7.5,
                    color=INK2)
    ax = fig.add_subplot(gs[1, 1])
    helix_bands(ax, label=False)
    ax.bar(ld.position, np.abs(ld.dpsi_deg), color=INK2, width=0.9, label="|psi cl7 - psi cl4|")
    ax.set(ylabel="|mean psi difference| (deg)", xlabel="ESR1 residue", xlim=(304, 551))
    ax.legend(fontsize=8, loc="upper left")
    fig.suptitle(f"Alignment-free check: backbone dihedrals alone separate cl7 from cl4 (silhouette {sil:.2f} in 2-D); "
                 f"dPC1 is almost one residue: a psi flip of {abs(ld.set_index('position').loc[537, 'dpsi_deg']):.0f} deg at the "
                 f"H11-H12 hinge residue 537", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y12_phipsi_pca", section="C. Why two clusters", task="Y12",
             shows="PCA of circular-encoded phi/psi over positions observed in >= 90% of chains; dPC1 loadings and mean psi differences along the sequence",
             takeaway=(f"dPC1 explains {ev[0] * 100:.0f}% of dihedral variance; largest loadings at "
                       f"{', '.join(map(str, top.position))}; cl7/cl4 silhouette {sil:.2f} without any superposition."),
             source="results/ext/phipsi_pca.csv, results/ext/phipsi_loadings.csv", script="scripts/ext_phipsi.py")
    write_numbers("Y12", dict(n_positions=int(keep.sum()), ev1=float(ev[0]), ev2=float(ev[1]), silhouette_2d=float(sil),
                              top_positions=top.position.tolist(), median_imputed=float(np.median(n_imp)),
                              n_chains_with_angles=int((np.isfinite(ang[..., 0]).sum((1)) > 0).sum())))
    print(keep.sum(), ev[:3], sil, top.position.tolist(), np.median(n_imp), (np.isfinite(ang[..., 0]).sum(1) == 0).sum())


if __name__ == "__main__":
    main()
