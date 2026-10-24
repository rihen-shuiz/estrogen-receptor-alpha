"""Tracker Y20: PCA of the superposed C-alpha coordinates (principal modes of the ESR1 LBD ensemble).

Input: S03 superposed C-alpha coordinates (results/alignment/aligned_coordinates.npz, core-v3 fit; 788 chains x 246
positions x 3, NaN where unobserved). Positions observed in >= 95% of chains are used (this keeps H12 538-544);
missing coordinates are imputed with the position's mean over all chains (count per chain reported).
PCA on the centred 3N-vector (sklearn, 10 components). Per-residue loading of a mode = length of its 3-D loading
vector at that position (sum over x, y, z squared, square-rooted). Each mode is converted to an amplitude in A:
the RMS displacement of the C-alpha set along the mode for one score SD.
Comparison with the RMSD-based maps (Y03): Spearman between cPC1 and PCA-of-RMSD PC1 / metric MDS 1, and between
pairwise distances in cPC1-3 space and the S04 RMSD.
Outputs: results/ext/coord_pca_scores.csv, results/ext/coord_pca_loadings.csv, results/ext/coord_pca_variance.csv,
figures/ext/Y20_coord_pca.png, figures/ext/Y20_coord_pca_loadings.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_coord_pca.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK2, ORANGE, RES, SEED, VIOLET, YELLOW, chain_table, coords, helix_bands,
                        mpl_setup, rmsd, save_fig, scatter_clusters, write_numbers)


def main():
    x, pos = coords()
    d, idx = rmsd()
    t = chain_table()
    obs = np.isfinite(x[..., 0])
    keep = obs.mean(0) >= 0.95
    xs = x[:, keep]
    n_imp = (~np.isfinite(xs[..., 0])).sum(1)
    mu = np.nanmean(xs, 0)
    xs = np.where(np.isfinite(xs), xs, mu)
    flat = xs.reshape(len(xs), -1)
    pca = PCA(10, random_state=SEED).fit(flat)
    z = pca.transform(flat)
    comp = pca.components_.reshape(10, keep.sum(), 3)
    amp = np.sqrt(pca.explained_variance_) / np.sqrt(keep.sum())  # RMS displacement per SD of score
    sc = pd.DataFrame({"chain_key": t.chain_key, "cpc1": z[:, 0], "cpc2": z[:, 1], "cpc3": z[:, 2], "cpc4": z[:, 3],
                       "n_imputed_positions": n_imp, "cluster": t.cluster})
    sc.to_csv(RES / "coord_pca_scores.csv", index=False)
    ld = pd.DataFrame({"position": pos[keep]})
    for k in range(4):
        ld[f"load{k + 1}"] = np.linalg.norm(comp[k], axis=1)
    ld.to_csv(RES / "coord_pca_loadings.csv", index=False)
    var = pd.DataFrame({"mode": np.arange(1, 11), "explained_variance_ratio": pca.explained_variance_ratio_,
                        "rms_amplitude_per_sd_A": amp})
    var.to_csv(RES / "coord_pca_variance.csv", index=False)
    emb = pd.read_csv(RES / "embedding_rmsd.csv")
    r_pca = spearmanr(z[:, 0], emb.pca1)[0]
    r_mds = spearmanr(z[:, 0], emb.mds1)[0]
    iu = np.triu_indices(len(z), 1)
    e3 = np.sqrt(((z[:, None, :3] - z[None, :, :3]) ** 2).sum(-1))
    r_dist = spearmanr(e3[iu], d[iu])[0]

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.8), gridspec_kw={"width_ratios": [0.8, 1, 1, 1]})
    ax = axs[0]
    ax.bar(var["mode"], var.explained_variance_ratio * 100, color=VIOLET)
    ax.plot(var["mode"], var.explained_variance_ratio.cumsum() * 100, "-o", color=INK2, ms=3, lw=1)
    for m, v in zip(var["mode"][:3], var.explained_variance_ratio[:3]):
        ax.text(m, v * 100 + 1, f"{v * 100:.0f}%", ha="center", fontsize=8)
    ax.set(xlabel="principal mode", ylabel="explained variance (%) / cumulative", title="Scree")
    ev = pca.explained_variance_ratio_
    for ax, (a, b) in zip(axs[1:3], ((0, 1), (1, 2))):
        scatter_clusters(ax, z[:, a], z[:, b], t.cluster, s=8, legend=(a == 0))
        ax.set(xlabel=f"cPC{a + 1} ({ev[a] * 100:.1f}%, A)", ylabel=f"cPC{b + 1} ({ev[b] * 100:.1f}%, A)")
    axs[1].set_title("Coordinate PCA: cPC1 vs cPC2")
    axs[2].set_title("cPC2 vs cPC3 (variation inside the clusters)")
    ax = axs[3]
    scatter_clusters(ax, emb.pca1, z[:, 0], t.cluster, s=8, legend=False)
    ax.set(xlabel="PC1 of RMSD rows (Y03)", ylabel="cPC1 of coordinates",
           title=f"Same axis: Spearman {r_pca:.2f}")
    fig.suptitle(f"PCA of superposed C-alpha coordinates ({keep.sum()} positions): mode 1 ({ev[0] * 100:.0f}%) is the "
                 f"cl7/cl4 switch and tracks the RMSD-based PC1; modes 2-3 describe variation inside the clusters",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y20_coord_pca", section="G. Coordinate PCA", task="Y20",
             shows="scree plot, projections of all 788 chains on coordinate modes 1-3, cPC1 against PC1 of the RMSD rows",
             takeaway=(f"cPC1 explains {ev[0] * 100:.0f}% of coordinate variance and correlates with the RMSD-based PC1 "
                       f"(Spearman {r_pca:.2f}); distances in cPC1-3 space track S04 RMSD (Spearman {r_dist:.2f})."),
             source="results/ext/coord_pca_scores.csv, results/ext/coord_pca_variance.csv", script="scripts/ext_coord_pca.py")
    fig, axs = plt.subplots(3, 1, figsize=(16, 8), sharex=True)
    for k, (ax, col) in enumerate(zip(axs, (VIOLET, YELLOW, GREY))):
        helix_bands(ax, label=(k == 0))
        ax.bar(ld.position, ld[f"load{k + 1}"], color=col, width=0.9)
        top = ld.sort_values(f"load{k + 1}", ascending=False).head(3)
        for r in top.itertuples():
            ax.annotate(str(r.position), (r.position, getattr(r, f"load{k + 1}")), xytext=(2, 2),
                        textcoords="offset points", fontsize=7.5, color=INK2)
        ax.set(ylabel=f"mode {k + 1} loading", xlim=(304, 551))
        ax.text(0.995, 0.9, f"mode {k + 1}: {ev[k] * 100:.1f}% variance, {amp[k]:.2f} A RMS per SD", transform=ax.transAxes,
                ha="right", fontsize=8.5, color=INK2)
    axs[-1].set_xlabel("ESR1 residue")
    fig.suptitle("Per-residue loadings of the first three coordinate modes: all three are H12 motions (mode 1 the "
                 "cl7/cl4 switch, modes 2-3 H12 variation within the states)", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y20_coord_pca_loadings", section="G. Coordinate PCA", task="Y20",
             shows="per-residue loading magnitude of coordinate modes 1-3 along the ESR1 sequence with helix bands",
             takeaway="Modes 1-3 load almost only on H12 (538-544) and the 528-537 loop; the helical core has near-zero loading in every mode.",
             source="results/ext/coord_pca_loadings.csv", script="scripts/ext_coord_pca.py")
    write_numbers("Y20", dict(n_positions=int(keep.sum()), ev=ev[:5].round(4).tolist(), amp=amp[:5].round(3).tolist(),
                              spearman_cpc1_rmsdpc1=float(r_pca), spearman_cpc1_mds1=float(r_mds),
                              spearman_dist3_rmsd=float(r_dist), median_imputed=float(np.median(n_imp)),
                              top_load1=ld.sort_values("load1", ascending=False).head(5).position.tolist(),
                              top_load2=ld.sort_values("load2", ascending=False).head(5).position.tolist(),
                              top_load3=ld.sort_values("load3", ascending=False).head(5).position.tolist()))
    print(ev[:5], r_pca, r_mds, r_dist, ld.sort_values("load2", ascending=False).head(5).position.tolist(),
          ld.sort_values("load3", ascending=False).head(5).position.tolist())


if __name__ == "__main__":
    main()
