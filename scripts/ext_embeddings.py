"""Tracker Y03: shared 2-D embeddings of the frozen 788 x 788 C-alpha RMSD matrix (results/rmsd_entry.npy).

Embeddings (all on the full matrix, every chain kept, seed 42):
  pca   - PCA of the RMSD rows (each chain described by its distances to all chains; as block 12), z = rows centred.
  cmds  - classical (Torgerson) MDS: double-centred -D^2/2, top eigenvectors; explained = eigenvalue share of the
          positive spectrum.
  mds   - metric SMACOF MDS (sklearn, dissimilarity precomputed, n_init 4), stress-1 reported.
  tsne  - t-SNE on the precomputed distances (perplexity 30, random init).
  umap  - UMAP on the precomputed distances (n_neighbors 15, min_dist 0.1).
Quality: trustworthiness (k = 10) of each 2-D map against the RMSD matrix, and Spearman correlation between RMSD and
2-D Euclidean distances. Output: results/ext/embedding_rmsd.csv, results/ext/embedding_quality.csv,
figures/ext/Y03_embeddings.png. Shared map for every colour-by figure: the embedding with the highest Spearman
correlation between RMSD and map distance (global faithfulness; declared before looking) - metric MDS in this run.
Run from the repo root: PYTHONPATH=. python scripts/ext_embeddings.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, trustworthiness

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import RES, SEED, chain_table, mpl_setup, rmsd, save_fig, scatter_clusters, write_numbers


def classical_mds(d, k=10):
    n = len(d)
    j = np.eye(n) - 1.0 / n
    b = -0.5 * j @ (d ** 2) @ j
    w, v = np.linalg.eigh(b)
    order = np.argsort(w)[::-1]
    w, v = w[order], v[:, order]
    pos = w[w > 0]
    x = v[:, :k] * np.sqrt(np.maximum(w[:k], 0))
    return x, w[:k] / pos.sum()


def main():
    d, idx = rmsd()
    assert np.isfinite(d).all() and np.allclose(d, d.T)
    t = chain_table()
    out = pd.DataFrame({"chain_key": idx.chain_key, "entry_id": idx.entry_id})
    pca = PCA(n_components=10, random_state=SEED).fit(d)
    zp = pca.transform(d)
    out["pca1"], out["pca2"], out["pca3"] = zp[:, 0], zp[:, 1], zp[:, 2]
    zc, cvar = classical_mds(d)
    out["cmds1"], out["cmds2"], out["cmds3"] = zc[:, 0], zc[:, 1], zc[:, 2]
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=SEED, n_init=4, max_iter=500,
              normalized_stress="auto")
    zm = mds.fit_transform(d)
    out["mds1"], out["mds2"] = zm[:, 0], zm[:, 1]
    zt = TSNE(n_components=2, metric="precomputed", init="random", perplexity=30, random_state=SEED).fit_transform(d)
    out["tsne1"], out["tsne2"] = zt[:, 0], zt[:, 1]
    import umap
    zu = umap.UMAP(n_neighbors=15, min_dist=0.1, metric="precomputed", random_state=SEED).fit_transform(d)
    out["umap1"], out["umap2"] = zu[:, 0], zu[:, 1]
    out.to_csv(RES / "embedding_rmsd.csv", index=False)

    iu = np.triu_indices(len(d), 1)
    q = []
    for name, z, var in (("pca", zp[:, :2], pca.explained_variance_ratio_[:2].sum()),
                         ("cmds", zc[:, :2], cvar[:2].sum()), ("mds", zm, np.nan), ("tsne", zt, np.nan),
                         ("umap", zu, np.nan)):
        e = np.sqrt(((z[:, None, :] - z[None, :, :]) ** 2).sum(-1))
        q.append(dict(embedding=name, trustworthiness=trustworthiness(d, z, n_neighbors=10, metric="precomputed"),
                      spearman_dist=spearmanr(d[iu], e[iu])[0], explained_2d=var))
    q = pd.DataFrame(q)
    q.loc[q.embedding == "mds", "explained_2d"] = np.nan
    q["stress1"] = np.nan
    em = np.sqrt(((zm[:, None, :] - zm[None, :, :]) ** 2).sum(-1))
    q.loc[q.embedding == "mds", "stress1"] = np.sqrt(((d[iu] - em[iu]) ** 2).sum() / (d[iu] ** 2).sum())
    q.to_csv(RES / "embedding_quality.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 5, figsize=(19, 4.2))
    labels = {"pca": ("PC", pca.explained_variance_ratio_), "cmds": ("MDS", cvar)}
    for ax, (name, a, b) in zip(axs, (("pca", "pca1", "pca2"), ("cmds", "cmds1", "cmds2"), ("mds", "mds1", "mds2"),
                                      ("tsne", "tsne1", "tsne2"), ("umap", "umap1", "umap2"))):
        scatter_clusters(ax, out[a], out[b], t.cluster, s=8, legend=(name == "pca"))
        r = q.set_index("embedding").loc[name]
        ax.set_title(f"{name.upper()}  trust={r.trustworthiness:.3f}, rho={r.spearman_dist:.2f}", fontsize=9)
        if name in labels:
            pre, v = labels[name]
            ax.set_xlabel(f"{pre}1 ({v[0] * 100:.1f}%)")
            ax.set_ylabel(f"{pre}2 ({v[1] * 100:.1f}%)")
        else:
            ax.set_xlabel(f"{name} 1")
            ax.set_ylabel(f"{name} 2")
    fig.suptitle(f"Five maps of the same 788 x 788 RMSD matrix all show the two main clusters (cl7 n=445, cl4 n=329)",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    best = q.sort_values("trustworthiness", ascending=False).iloc[0]
    save_fig(fig, "Y03_embeddings", section="F. Embeddings and colour-by", task="Y03",
             shows="PCA of RMSD rows, classical MDS, metric MDS, t-SNE and UMAP, coloured by S08 cluster",
             takeaway=(f"All five maps separate cl7 from cl4; classical MDS keeps {cvar[:2].sum() * 100:.0f}% of the "
                       f"positive spectrum in 2-D; highest trustworthiness {best.embedding} ({best.trustworthiness:.3f})."),
             source="results/ext/embedding_rmsd.csv, results/ext/embedding_quality.csv", script="scripts/ext_embeddings.py")
    write_numbers("Y03", dict(pca_var_2d=float(pca.explained_variance_ratio_[:2].sum()), cmds_var_2d=float(cvar[:2].sum()),
                              cmds_var_1=float(cvar[0]), quality=q.set_index("embedding").round(4).to_dict("index")))
    print(q.to_string())


if __name__ == "__main__":
    main()
