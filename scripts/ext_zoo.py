"""Tracker Y04: structural clustering zoo on the frozen 788 x 788 C-alpha RMSD matrix, automatically evaluated.

Inputs: results/rmsd_entry.npy (+ results/rmsd_entry_index.csv). Methods that need vectors use the 10-D classical
MDS embedding of the same matrix (computed here, deterministic; 10 dims hold most of the positive spectrum).
Families and settings (every setting is one run):
  hier_average / hier_complete / hier_single : SciPy linkage on the RMSD matrix, K = 2..8 (maxclust)
  hier_ward    : Ward linkage on the 10-D MDS vectors, K = 2..8
  kmeans       : k-means (n_init 20) on MDS vectors, K = 2..8
  kmedoids_pam : PAM (BUILD + SWAP, implemented below) on the RMSD matrix, K = 2..8
  spectral     : affinity exp(-d^2 / 2 s^2), s in {0.5, 1, 2} A, K in {2, 3, 4, 6}
  dbscan       : precomputed, eps in {0.75, 1.0, 1.5, 2.0} A, min_samples 5
  hdbscan      : precomputed, min_cluster_size in {5, 10, 20}
  optics       : precomputed, min_samples in {5, 10}, xi 0.05
  affinity_prop: similarity -d^2, preference = median / minimum similarity
  meanshift    : MDS vectors, bandwidth = quantile {0.1, 0.2, 0.3} of distances
  birch        : MDS vectors, n_clusters in {2, 3, 4, 6}, threshold 0.5
  gmm          : full-covariance Gaussian mixture on MDS vectors, K = 2..8 (BIC recorded)
Metrics per run: silhouette on the precomputed RMSD (noise points, label -1, excluded), Davies-Bouldin and
Calinski-Harabasz on the MDS vectors (non-noise), noise fraction, and bootstrap stability: 50 subsamples of 80% of
the chains (seed 42), the same method re-run on the sub-matrix, ARI against the full-run labels on the subsample
(noise kept as its own label); mean and SD reported.
Selection rule (declared in docs/PLAN_EXTENDED.md before running): eligible = n_clusters >= 2 and noise <= 0.10;
score = mean(percentile rank of silhouette, percentile rank of stability ARI); ties -> fewer clusters.
Outputs: results/ext/zoo_runs.csv, results/ext/zoo_labels.csv, results/ext/zoo_ranking.csv,
results/ext/zoo_mds10.npy, figures/ext/Y04_zoo_metrics.png, figures/ext/Y04_zoo_grid.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_zoo.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import (DBSCAN, HDBSCAN, OPTICS, AffinityPropagation, Birch, KMeans, MeanShift,
                             SpectralClustering)
from sklearn.metrics import adjusted_rand_score, calinski_harabasz_score, davies_bouldin_score, silhouette_score
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK2, ORANGE, RES, ROOT, SEED, chain_table, mpl_setup, rmsd,
                        save_fig, write_numbers)

warnings.filterwarnings("ignore")
N_BOOT, FRAC = 50, 0.8


# ---------------------------------------------------------------- PAM
def pam(d: np.ndarray, k: int, max_iter=100) -> np.ndarray:
    """Partitioning Around Medoids (Kaufman & Rousseeuw): greedy BUILD, then best-improvement SWAP until no gain."""
    n = len(d)
    med = [int(np.argmin(d.sum(1)))]
    for _ in range(1, k):
        cur = d[:, med].min(1)
        gain = np.maximum(cur[:, None] - d, 0).sum(0)
        gain[med] = -1
        med.append(int(np.argmax(gain)))
    med = np.array(med)
    for _ in range(max_iter):
        dm = d[:, med]
        order = np.argsort(dm, 1)
        nearest = order[:, 0]
        d1 = dm[np.arange(n), nearest]
        d2 = dm[np.arange(n), order[:, 1]] if k > 1 else np.full(n, np.inf)
        cost = d1.sum()
        best = (0.0, None, None)
        for j in range(k):
            base = np.where(nearest == j, d2, d1)
            new = np.minimum(d, base[:, None]).sum(0) - cost
            new[med] = np.inf
            h = int(np.argmin(new))
            if new[h] < best[0] - 1e-9:
                best = (new[h], j, h)
        if best[1] is None:
            break
        med[best[1]] = best[2]
    return np.argmin(d[:, med], 1)


def classical_mds(d, k=10):
    n = len(d)
    j = np.eye(n) - 1.0 / n
    b = -0.5 * j @ (d ** 2) @ j
    w, v = np.linalg.eigh(b)
    o = np.argsort(w)[::-1][:k]
    return v[:, o] * np.sqrt(np.maximum(w[o], 0))


# ---------------------------------------------------------------- run definitions
def run_list():
    runs = []
    for k in range(2, 9):
        for fam in ("hier_average", "hier_complete", "hier_single"):
            runs.append((fam, {"k": k}))
        runs.append(("hier_ward", {"k": k}))
        runs.append(("kmeans", {"k": k}))
        runs.append(("kmedoids_pam", {"k": k}))
        runs.append(("gmm", {"k": k}))
    for s in (0.5, 1.0, 2.0):
        for k in (2, 3, 4, 6):
            runs.append(("spectral", {"sigma": s, "k": k}))
    for e in (0.3, 0.4, 0.5, 0.75, 1.0, 1.5):
        runs.append(("dbscan", {"eps": e, "min_samples": 5}))
    for m in (5, 10, 20):
        runs.append(("hdbscan", {"min_cluster_size": m}))
    for m in (5, 10):
        runs.append(("optics", {"min_samples": m}))
    for p in ("median", "min"):
        runs.append(("affinity_prop", {"preference": p}))
    for q in (0.1, 0.2, 0.3):
        runs.append(("meanshift", {"quantile": q}))
    for k in (2, 3, 4, 6):
        runs.append(("birch", {"k": k}))
    return runs


def fit(fam: str, p: dict, d: np.ndarray, x: np.ndarray, seed=SEED) -> np.ndarray:
    if fam.startswith("hier_") and fam != "hier_ward":
        z = linkage(squareform(d, checks=False), fam.split("_")[1])
        return fcluster(z, p["k"], "maxclust") - 1
    if fam == "hier_ward":
        return fcluster(linkage(x, "ward"), p["k"], "maxclust") - 1
    if fam == "kmeans":
        return KMeans(p["k"], n_init=20, random_state=seed).fit_predict(x)
    if fam == "kmedoids_pam":
        return pam(d, p["k"])
    if fam == "gmm":
        return GaussianMixture(p["k"], covariance_type="full", random_state=seed, n_init=3, reg_covar=1e-4).fit_predict(x)
    if fam == "spectral":
        a = np.exp(-d ** 2 / (2 * p["sigma"] ** 2))
        return SpectralClustering(p["k"], affinity="precomputed", random_state=seed, assign_labels="cluster_qr").fit_predict(a)
    if fam == "dbscan":
        return DBSCAN(eps=p["eps"], min_samples=p["min_samples"], metric="precomputed").fit_predict(d)
    if fam == "hdbscan":
        return HDBSCAN(min_cluster_size=p["min_cluster_size"], metric="precomputed").fit_predict(d.copy())  # HDBSCAN edits its input
    if fam == "optics":
        return OPTICS(min_samples=p["min_samples"], metric="precomputed", xi=0.05).fit_predict(d)
    if fam == "affinity_prop":
        s = -d ** 2
        pref = np.median(s) if p["preference"] == "median" else s.min()
        return AffinityPropagation(affinity="precomputed", preference=pref, random_state=seed, max_iter=1000,
                                   damping=0.9).fit_predict(s)
    if fam == "meanshift":
        iu = np.triu_indices(len(x), 1)
        e = np.sqrt(((x[:, None] - x[None]) ** 2).sum(-1))[iu]
        return MeanShift(bandwidth=float(np.quantile(e, p["quantile"])), bin_seeding=True).fit_predict(x)
    if fam == "birch":
        return Birch(n_clusters=p["k"], threshold=0.5).fit_predict(x)
    raise ValueError(fam)


def pstr(p):
    return ",".join(f"{k}={v}" for k, v in p.items())


def evaluate(i, fam, p, d, x):
    lab = fit(fam, p, d, x)
    ok = lab >= 0
    n_cl = len(set(lab[ok]))
    row = dict(run_id=f"{fam}|{pstr(p)}", method=fam, params=pstr(p), n_clusters=n_cl, noise_fraction=1 - ok.mean(),
               largest_fraction=np.bincount(lab[ok]).max() / len(lab) if ok.any() else np.nan)
    if n_cl >= 2 and ok.sum() > n_cl:
        row["silhouette"] = silhouette_score(d[np.ix_(ok, ok)], lab[ok], metric="precomputed")
        row["davies_bouldin"] = davies_bouldin_score(x[ok], lab[ok])
        row["calinski_harabasz"] = calinski_harabasz_score(x[ok], lab[ok])
    else:
        row["silhouette"] = row["davies_bouldin"] = row["calinski_harabasz"] = np.nan
    rng = np.random.default_rng(SEED + i)
    aris = []
    for _ in range(N_BOOT):
        sub = np.sort(rng.choice(len(d), int(FRAC * len(d)), replace=False))
        try:
            lb = fit(fam, p, d[np.ix_(sub, sub)], x[sub])
            aris.append(adjusted_rand_score(lab[sub], lb))
        except Exception:
            aris.append(np.nan)
    row["stability_ari"] = np.nanmean(aris)
    row["stability_sd"] = np.nanstd(aris)
    if fam == "gmm":
        row["bic"] = GaussianMixture(p["k"], covariance_type="full", random_state=SEED, n_init=3,
                                     reg_covar=1e-4).fit(x).bic(x)
    return row, lab


def align_to_reference(lab, ref):
    """Map each cluster of `lab` to the S08 cluster it overlaps most (for consistent colours)."""
    out = {}
    for c in sorted(set(lab)):
        if c < 0:
            out[c] = -1
            continue
        vals, cnt = np.unique(ref[lab == c], return_counts=True)
        out[c] = int(vals[np.argmax(cnt)])
    return out


def main():
    d, idx = rmsd()
    x = classical_mds(d, 10)
    np.save(RES / "zoo_mds10.npy", x)
    runs = run_list()
    res = Parallel(n_jobs=-1, verbose=0)(delayed(evaluate)(i, f, p, d, x) for i, (f, p) in enumerate(runs))
    rows = [r for r, _ in res]
    df = pd.DataFrame(rows)
    labels = pd.DataFrame({"chain_key": idx.chain_key})
    for (r, lab) in res:
        labels[r["run_id"]] = lab
    labels.to_csv(RES / "zoo_labels.csv", index=False)
    # selection
    elig = (df.n_clusters >= 2) & (df.noise_fraction <= 0.10) & df.silhouette.notna()
    df["eligible"] = elig
    e = df[elig].copy()
    e["rank_sil"] = e.silhouette.rank(pct=True)
    e["rank_stab"] = e.stability_ari.rank(pct=True)
    e["score"] = (e.rank_sil + e.rank_stab) / 2
    e = e.sort_values(["score", "n_clusters"], ascending=[False, True]).reset_index(drop=True)
    e.insert(0, "rank", np.arange(1, len(e) + 1))
    e[["rank", "run_id", "method", "params", "n_clusters", "score", "silhouette", "stability_ari", "davies_bouldin",
       "calinski_harabasz", "noise_fraction"]].to_csv(RES / "zoo_ranking.csv", index=False)
    df = df.merge(e[["run_id", "score", "rank"]], on="run_id", how="left")
    t = chain_table()
    ref = t.cluster.values
    ref_ari = {r: adjusted_rand_score(ref, labels[r]) for r in df.run_id}
    df["ari_vs_s08"] = df.run_id.map(ref_ari)
    main2 = np.isin(ref, (7, 4))
    df["main_split_ari"] = [adjusted_rand_score(ref[main2], labels[r].values[main2]) for r in df.run_id]
    df.to_csv(RES / "zoo_runs.csv", index=False)
    best = e.iloc[0]
    fam_best = e.sort_values("score", ascending=False).groupby("method").head(1).set_index("method")
    # families with no eligible run (too much noise) are still drawn: their best-silhouette run, marked not eligible
    rest = df[~df.method.isin(fam_best.index)].sort_values(["n_clusters", "silhouette"], ascending=False)
    fam_all = pd.concat([fam_best, rest.groupby("method").head(1).set_index("method")])

    # ---------------- figure 1: metrics
    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(17, 6.2), gridspec_kw={"width_ratios": [1.15, 1, 1]})
    ax = axs[0]
    ax.scatter(df.silhouette, df.stability_ari, s=22, c=GREY, edgecolor="white", lw=0.4, label="all runs (n=%d)" % len(df))
    top = e.head(8)
    ax.scatter(top.silhouette, top.stability_ari, s=46, c=BLUE, edgecolor="white", lw=0.6, label="top 8 by score")
    ax.scatter([best.silhouette], [best.stability_ari], s=120, facecolor="none", edgecolor=ORANGE, lw=2,
               label="selected")
    for _, r in top.iterrows():
        ax.annotate(f"{r['rank']}", (r.silhouette, r.stability_ari), xytext=(4, 2), textcoords="offset points",
                    fontsize=8, color=INK2)
    ax.set(xlabel="silhouette (precomputed RMSD)", ylabel="bootstrap stability (mean ARI, 50 x 80%)",
           title="Silhouette vs stability for every run")
    ax.legend(loc="lower left", fontsize=8)
    ax = axs[1]
    fb = fam_best.sort_values("score")
    yy = np.arange(len(fb))
    ax.barh(yy, fb.score, color=[ORANGE if m == best.method else BLUE for m in fb.index], height=0.6)
    ax.set_yticks(yy, [f"{m}  ({fb.loc[m, 'params']})" for m in fb.index], fontsize=7.5)
    for y, (m, r) in zip(yy, fb.iterrows()):
        ax.text(r.score + 0.01, y, f"K={r.n_clusters}  sil {r.silhouette:.2f}  stab {r.stability_ari:.2f}",
                va="center", fontsize=7, color=INK2)
    ax.set_xlim(0, 1.45)
    ax.set(xlabel="selection score (mean percentile rank)", title="Best setting of each family")
    ax = axs[2]
    kf = ["hier_average", "hier_complete", "hier_single", "hier_ward", "kmeans", "kmedoids_pam", "gmm"]
    for c, f in zip(CATEGORICAL, kf):
        s = df[df.method == f].copy()
        s["k"] = s.params.str.extract(r"k=(\d+)").astype(int)
        s = s.sort_values("k")
        ax.plot(s.k, s.silhouette, "-o", color=c, ms=4, lw=1.5, label=f)
    ax.set(xlabel="number of clusters K", ylabel="silhouette (precomputed RMSD)",
           title="Silhouette peaks at K = 2 for every K-based family")
    ax.legend(fontsize=7.5, ncol=2)
    fig.suptitle(f"Clustering zoo: {len(df)} runs, {df.method.nunique()} families; the declared rule selects {best.run_id} "
                 f"(silhouette {best.silhouette:.2f}, stability {best.stability_ari:.2f})", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    save_fig(fig, "Y04_zoo_metrics", section="A. Structural clustering zoo", task="Y04",
             shows="silhouette vs bootstrap stability for all runs, best setting per family, silhouette vs K",
             takeaway=(f"Selected by the declared rule: {best.run_id} (K={best.n_clusters}, silhouette "
                       f"{best.silhouette:.2f}, stability ARI {best.stability_ari:.2f}); "
                       f"{(e.head(10).n_clusters == 2).sum()} of the top 10 runs have K=2."),
             source="results/ext/zoo_runs.csv, results/ext/zoo_ranking.csv", script="scripts/ext_zoo.py")

    # ---------------- figure 2: grid of best setting per family on the shared MDS map
    emb = pd.read_csv(RES / "embedding_rmsd.csv")
    fams = list(fam_best.sort_values("score", ascending=False).index) + [m for m in fam_all.index if m not in fam_best.index]
    ncol = 5
    nrow = int(np.ceil((len(fams) + 1) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.3 * nrow))
    axs = axs.ravel()
    extra = [AQUA, CATEGORICAL[3], CATEGORICAL[4], CATEGORICAL[5], CATEGORICAL[6], CATEGORICAL[7]]

    def draw(ax, lab, title):
        m = align_to_reference(lab, ref)
        used = {}
        cols = []
        for l in lab:
            tgt = m[l]
            if l < 0:
                cols.append(GREY)
            elif tgt == 7 and 7 not in used.values():
                used[l] = 7; cols.append(BLUE)
            elif tgt == 4 and 4 not in used.values():
                used[l] = 4; cols.append(ORANGE)
            elif l in used:
                cols.append(BLUE if used[l] == 7 else ORANGE)
            else:
                cols.append(None)
        others = sorted({l for l, c in zip(lab, cols) if c is None})
        omap = {l: extra[i % len(extra)] for i, l in enumerate(others)}
        cols = [c if c is not None else omap[l] for l, c in zip(lab, cols)]
        ax.scatter(emb.mds1, emb.mds2, c=cols, s=5, lw=0)
        ax.set_title(title, fontsize=8)
        ax.set_xticks([]); ax.set_yticks([])

    draw(axs[0], ref, f"reference S08 (avg, 1.5 A cut), K={len(set(ref))}")
    for ax, f in zip(axs[1:], fams):
        r = fam_all.loc[f]
        tag = "" if f in fam_best.index else f" [noise {r.noise_fraction:.0%}; not eligible]"
        draw(ax, labels[r.run_id].values, f"{f} {r.params}{tag}\nK={r.n_clusters} sil={r.silhouette:.2f} "
                                         f"stab={r.stability_ari:.2f} ARI(S08)={ref_ari[r.run_id]:.2f}")
    for ax in axs[len(fams) + 1:]:
        ax.axis("off")
    fig.suptitle("Best setting of every family on the shared MDS map: the cl7 (blue) / cl4 (orange) split is found by "
                 "almost all methods (grey = noise)", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    agree = (df[df.n_clusters >= 2].main_split_ari > 0.9).sum()
    save_fig(fig, "Y04_zoo_grid", section="A. Structural clustering zoo", task="Y04",
             shows="best setting of each method family drawn on the shared metric-MDS map, colours matched to S08",
             takeaway=(f"{agree} of {(df.n_clusters >= 2).sum()} runs with >=2 clusters reproduce the cl7/cl4 split "
                       f"with ARI > 0.9 on the 774 main-cluster chains."),
             source="results/ext/zoo_labels.csv, results/ext/embedding_rmsd.csv", script="scripts/ext_zoo.py")
    write_numbers("Y04", dict(n_runs=len(df), n_families=int(df.method.nunique()), n_eligible=int(elig.sum()),
                              best_run=best.run_id, best_k=int(best.n_clusters), best_sil=float(best.silhouette),
                              best_stab=float(best.stability_ari), best_ari_s08=float(ref_ari[best.run_id]),
                              n_main_split_agree=int(agree), n_runs_k2plus=int((df.n_clusters >= 2).sum()),
                              top10=e.head(10).run_id.tolist()))
    print(e.head(15)[["rank", "run_id", "n_clusters", "score", "silhouette", "stability_ari"]].to_string())
    print(fam_best[["params", "n_clusters", "silhouette", "stability_ari", "score"]].to_string())


if __name__ == "__main__":
    main()
