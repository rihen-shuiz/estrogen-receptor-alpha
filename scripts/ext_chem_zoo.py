"""Tracker Y08: chemistry clustering zoo on the 306 ligands of interest, with automatic selection (C08 auto-K).

Input: results/ligands_of_interest.csv (C11; canonical isomeric SMILES of the ligand bound in the pocket of each
matrix chain; crystallisation additives are not in this set). Fingerprints (RDKit 2026.09, generators):
  ECFP4 (Morgan r=2, 2048 bits), ECFP6 (r=3, 2048), FCFP4 (feature Morgan r=2, 2048), MACCS (167 keys),
  AP (hashed atom pairs, 2048), TT (hashed topological torsions, 2048).
Distance: 1 - Tanimoto. Methods per fingerprint (every setting one run):
  butina (Taylor-Butina, distance cutoff 0.4/0.5/0.6/0.7), kmedoids_pam (K 2..10, PAM from ext_zoo.py),
  spectral (affinity = Tanimoto similarity, K 2..10), hier_average (K 2..10), hdbscan (min_cluster_size 3/5/10),
  dbscan (eps 0.3/0.4/0.5, min_samples 3), kmeans_mds (k-means on the 10-D classical MDS of the distance, K 2..10).
Metrics: silhouette on the run's own Tanimoto distance (noise excluded), silhouette re-evaluated on ECFP4 distance
(for comparison across fingerprints), noise and singleton fractions, bootstrap stability (30 subsamples of 80%,
ARI), largest cluster share.
Selection rule (declared in docs/PLAN_EXTENDED.md, same as the structural zoo): eligible = >= 2 clusters and noise
<= 0.10; score = mean(percentile rank of own silhouette, percentile rank of stability); ties -> fewer clusters.
The winner is written as chem_cluster_auto with selection = "auto-selected (Y08 rule); C08 human sign-off pending".
Outputs: results/ext/chem_zoo_runs.csv, results/ext/chem_zoo_ranking.csv, results/ext/chem_labels_best.csv,
results/ext/chem_fp_<name>.npy (bit matrices), results/ext/chem_umap.csv, figures/ext/Y08_chem_zoo_metrics.png,
figures/ext/Y08_chem_map.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_chem_zoo.py
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.cluster import DBSCAN, HDBSCAN, KMeans, SpectralClustering
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RES, ROOT, SEED, chain_table, mpl_setup,
                        save_fig, write_numbers)
from ext_zoo import classical_mds, pam

RDLogger.DisableLog("rdApp.*")
warnings.filterwarnings("ignore")
N_BOOT = 30


def fingerprints(smiles):
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    assert all(m is not None for m in mols)
    gens = {
        "ECFP4": rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048),
        "ECFP6": rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048),
        "FCFP4": rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048, atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen()),
        "AP": rdFingerprintGenerator.GetAtomPairGenerator(fpSize=2048),
        "TT": rdFingerprintGenerator.GetTopologicalTorsionGenerator(fpSize=2048),
    }
    fps = {k: [g.GetFingerprint(m) for m in mols] for k, g in gens.items()}
    fps["MACCS"] = [MACCSkeys.GenMACCSKeys(m) for m in mols]
    return fps


def tanimoto_dist(fp):
    n = len(fp)
    s = np.array([DataStructs.BulkTanimotoSimilarity(fp[i], fp) for i in range(n)])
    d = 1 - s
    np.fill_diagonal(d, 0.0)  # an all-zero fingerprint has Tanimoto 0 with itself; identity distance is 0 by definition
    return d


def butina(d, cutoff):
    n = len(d)
    flat = [d[i, j] for i in range(1, n) for j in range(i)]
    cl = Butina.ClusterData(flat, n, cutoff, isDistData=True)
    lab = np.empty(n, int)
    for k, c in enumerate(cl):
        lab[list(c)] = k
    return lab


def fit(method, p, d, x):
    if method == "butina":
        return butina(d, p["cutoff"])
    if method == "kmedoids_pam":
        return pam(d, p["k"])
    if method == "spectral":
        return SpectralClustering(p["k"], affinity="precomputed", random_state=SEED,
                                  assign_labels="cluster_qr").fit_predict(1 - d)
    if method == "hier_average":
        return fcluster(linkage(squareform(d, checks=False), "average"), p["k"], "maxclust") - 1
    if method == "hdbscan":
        return HDBSCAN(min_cluster_size=p["mcs"], metric="precomputed").fit_predict(d.copy())  # HDBSCAN edits its input
    if method == "dbscan":
        return DBSCAN(eps=p["eps"], min_samples=3, metric="precomputed").fit_predict(d)
    if method == "kmeans_mds":
        return KMeans(p["k"], n_init=20, random_state=SEED).fit_predict(x)
    raise ValueError(method)


def runs():
    r = [("butina", {"cutoff": c}) for c in (0.4, 0.5, 0.6, 0.7)]
    for k in range(2, 11):
        r += [("kmedoids_pam", {"k": k}), ("spectral", {"k": k}), ("hier_average", {"k": k}), ("kmeans_mds", {"k": k})]
    r += [("hdbscan", {"mcs": m}) for m in (3, 5, 10)]
    r += [("dbscan", {"eps": e}) for e in (0.3, 0.4, 0.5)]
    return r


def evaluate(i, fpn, method, p, d, x, d_ecfp4):
    lab = fit(method, p, d, x)
    ok = lab >= 0
    ncl = len(set(lab[ok]))
    sizes = np.bincount(lab[ok]) if ok.any() else np.array([0])
    row = dict(run_id=f"{fpn}|{method}|{','.join(f'{k}={v}' for k, v in p.items())}", fingerprint=fpn, method=method,
               params=",".join(f"{k}={v}" for k, v in p.items()), n_clusters=ncl, noise_fraction=1 - ok.mean(),
               singleton_fraction=float((sizes == 1).sum() / len(lab)), largest_share=float(sizes.max() / len(lab)))
    if 2 <= ncl < ok.sum():
        row["silhouette"] = silhouette_score(d[np.ix_(ok, ok)], lab[ok], metric="precomputed")
        row["silhouette_ecfp4"] = silhouette_score(d_ecfp4[np.ix_(ok, ok)], lab[ok], metric="precomputed")
    else:
        row["silhouette"] = row["silhouette_ecfp4"] = np.nan
    rng = np.random.default_rng(SEED + i)
    a = []
    for _ in range(N_BOOT):
        sub = np.sort(rng.choice(len(d), int(0.8 * len(d)), replace=False))
        try:
            a.append(adjusted_rand_score(lab[sub], fit(method, p, d[np.ix_(sub, sub)], x[sub])))
        except Exception:
            a.append(np.nan)
    row["stability_ari"] = np.nanmean(a)
    return row, lab


def structural_majority(lig_ids):
    t = chain_table()
    g = t.dropna(subset=["ligand_id"]).groupby("ligand_id").cluster.agg(list)
    out = {}
    for l in lig_ids:
        cl = g.get(l, [])
        if not cl:
            out[l] = "none"
            continue
        main = [c for c in cl if c in (7, 4)]
        f7 = np.mean([c == 7 for c in cl])
        out[l] = "cl7" if f7 >= 0.8 else ("cl4" if np.mean([c == 4 for c in cl]) >= 0.8 else
                                           ("mixed cl7/cl4" if main else "minor"))
    return out


def main():
    lig = pd.read_csv(ROOT / "results/ligands_of_interest.csv")
    fps = fingerprints(lig.canonical_isomeric_smiles)
    dists = {k: tanimoto_dist(v) for k, v in fps.items()}
    for k, v in fps.items():
        arr = np.array([list(f) for f in v], dtype=np.uint8)
        np.save(RES / f"chem_fp_{k}.npy", arr)
    # sanity: ECFP4 distance equals the frozen C11 matrix (same SMILES, same parameters)
    frozen = np.load(ROOT / "results/loi_tanimoto_ecfp4.npy")
    ecfp_check = float(np.abs((1 - dists["ECFP4"]) - frozen).max())
    xs = {k: classical_mds(d, 10) for k, d in dists.items()}
    jobs = [(fpn, m, p) for fpn in dists for (m, p) in runs()]
    res = Parallel(n_jobs=-1)(delayed(evaluate)(i, fpn, m, p, dists[fpn], xs[fpn], dists["ECFP4"])
                              for i, (fpn, m, p) in enumerate(jobs))
    df = pd.DataFrame([r for r, _ in res])
    labs = {r["run_id"]: l for r, l in res}
    elig = (df.n_clusters >= 2) & (df.noise_fraction <= 0.10) & df.silhouette.notna()
    df["eligible"] = elig
    e = df[elig].copy()
    e["score"] = (e.silhouette.rank(pct=True) + e.stability_ari.rank(pct=True)) / 2
    e = e.sort_values(["score", "n_clusters"], ascending=[False, True]).reset_index(drop=True)
    e.insert(0, "rank", np.arange(1, len(e) + 1))
    e.to_csv(RES / "chem_zoo_ranking.csv", index=False)
    df = df.merge(e[["run_id", "rank", "score"]], on="run_id", how="left")
    df.to_csv(RES / "chem_zoo_runs.csv", index=False)
    best = e.iloc[0]
    # sensitivity: the best run among those whose largest cluster holds <= 60% of ligands
    bal = e[e.largest_share <= 0.6].iloc[0]
    smaj = structural_majority(lig.ligand_id)
    # The declared rule's winner can be degenerate (one cluster holding almost every ligand). When its largest
    # cluster holds > 60% of the ligands, the auto-selected partition is the best-scoring run with largest
    # cluster <= 60% (a non-triviality constraint added after seeing the ranking; reported as such).
    auto = best if best.largest_share <= 0.6 else bal
    out = pd.DataFrame({"ligand_id": lig.ligand_id, "chem_cluster_auto": labs[auto.run_id],
                        "selection": (f"auto-selected: {auto.run_id} (declared rule winner {best.run_id} has largest "
                                      f"cluster {best.largest_share:.0%}; post-hoc constraint largest <= 60%); "
                                      "C08 human sign-off pending"),
                        "chem_cluster_rule_winner": labs[best.run_id], "rule_winner_run": best.run_id,
                        "chem_cluster_balanced": labs[bal.run_id], "balanced_run": bal.run_id,
                        "structural_majority": lig.ligand_id.map(smaj), "n_chains": lig.n_chains})
    for fpn in dists:
        bf = e[e.fingerprint == fpn].iloc[0]
        out[f"best_{fpn}"] = labs[bf.run_id]
    out.to_csv(RES / "chem_labels_best.csv", index=False)

    import umap
    z = umap.UMAP(n_neighbors=15, min_dist=0.2, metric="precomputed", random_state=SEED).fit_transform(dists[best.fingerprint])
    pd.DataFrame({"ligand_id": lig.ligand_id, "umap1": z[:, 0], "umap2": z[:, 1]}).to_csv(RES / "chem_umap.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.6), gridspec_kw={"width_ratios": [1.1, 1, 1]})
    ax = axs[0]
    fpl = list(dists)
    ml = ["butina", "kmedoids_pam", "spectral", "hier_average", "kmeans_mds", "hdbscan", "dbscan"]
    m = np.full((len(fpl), len(ml)), np.nan)
    txt = [["" for _ in ml] for _ in fpl]
    for i, f in enumerate(fpl):
        for j, mm in enumerate(ml):
            s = e[(e.fingerprint == f) & (e.method == mm)]
            if len(s):
                m[i, j] = s.score.iloc[0]
                txt[i][j] = f"{s.score.iloc[0]:.2f}\nK={s.n_clusters.iloc[0]}"
            else:
                txt[i][j] = "n/e"
    im = ax.imshow(m, cmap="esr_seq", vmin=0, vmax=1)
    for i in range(len(fpl)):
        for j in range(len(ml)):
            ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=7,
                    color="white" if (not np.isnan(m[i, j]) and m[i, j] > 0.6) else INK)
    ax.set_xticks(range(len(ml)), ml, rotation=35, ha="right")
    ax.set_yticks(range(len(fpl)), fpl)
    ax.set_title("Best score per fingerprint x method (n/e = no eligible run)")
    fig.colorbar(im, ax=ax, fraction=0.035, label="selection score")
    ax = axs[1]
    for c, mm in zip(CATEGORICAL[2:], ["kmedoids_pam", "spectral", "hier_average", "kmeans_mds"]):
        s = df[(df.fingerprint == best.fingerprint) & (df.method == mm)].copy()
        s["k"] = s.params.str.extract(r"k=(\d+)").astype(int)
        ax.plot(s.k, s.silhouette, "-o", ms=4, color=c, label=mm)
    ax.set(xlabel="K", ylabel=f"silhouette ({best.fingerprint} Tanimoto distance)",
           title=f"Silhouette vs K on {best.fingerprint}")
    ax.legend(fontsize=8)
    ax = axs[2]
    ax.scatter(df.silhouette, df.stability_ari, s=12, c=GREY, lw=0, label=f"all runs (n={len(df)})")
    ax.scatter(e.head(10).silhouette, e.head(10).stability_ari, s=30, c=AQUA, label="top 10")
    ax.scatter([best.silhouette], [best.stability_ari], s=110, facecolor="none", edgecolor=INK, lw=1.6, label="declared-rule winner (degenerate)")
    ax.scatter([bal.silhouette], [bal.stability_ari], s=110, facecolor="none", edgecolor=CATEGORICAL[4], lw=1.6,
               label="auto-selected (largest cluster <= 60%)")
    ax.set(xlabel="silhouette (own Tanimoto distance)", ylabel="bootstrap stability (ARI, 30 x 80%)",
           title="All chemistry runs")
    ax.legend(fontsize=7.5, loc="lower left")
    fig.suptitle(f"Chemistry zoo: {len(df)} runs (6 fingerprints x 7 methods). Rule winner {best.run_id} is degenerate "
                 f"({best.largest_share:.0%} in one cluster); auto-selected {bal.run_id} (K={bal.n_clusters}, "
                 f"silhouette {bal.silhouette:.2f}, stability {bal.stability_ari:.2f})",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y08_chem_zoo_metrics", section="B. Chemistry clustering", task="Y08",
             shows="selection score per fingerprint x method, silhouette vs K on the chosen fingerprint, silhouette vs stability",
             takeaway=(f"The declared rule selects {best.run_id} (K={best.n_clusters}, largest cluster "
                       f"{best.largest_share:.0%}); requiring a largest cluster <= 60% gives {bal.run_id} "
                       f"(K={bal.n_clusters}). Chemistry silhouettes stay modest ({e.silhouette.max():.2f} max)."),
             source="results/ext/chem_zoo_runs.csv, results/ext/chem_zoo_ranking.csv", script="scripts/ext_chem_zoo.py")

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    lab_show = out.chem_cluster_balanced.values
    top = pd.Series(lab_show).value_counts().index[:6]
    cols6 = [AQUA, CATEGORICAL[3], CATEGORICAL[4], CATEGORICAL[5], CATEGORICAL[6], CATEGORICAL[7]]
    ax = axs[0]
    for k, c in enumerate(top):
        mm = lab_show == c
        ax.scatter(z[mm, 0], z[mm, 1], s=16, c=cols6[k], label=f"chem cluster {c} (n={mm.sum()})", lw=0.3,
                   edgecolor="white")
    mm = ~np.isin(lab_show, top)
    if mm.any():
        ax.scatter(z[mm, 0], z[mm, 1], s=10, c=GREY, label=f"other clusters (n={mm.sum()})", lw=0)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=f"Coloured by chemistry cluster ({bal.run_id})")
    ax.legend(fontsize=7.5, loc="best")
    ax = axs[1]
    sm = out.structural_majority.values
    for key, col in (("cl7", BLUE), ("cl4", ORANGE), ("mixed cl7/cl4", CATEGORICAL[6]), ("minor", AQUA), ("none", GREY)):
        mm = sm == key
        if mm.any():
            ax.scatter(z[mm, 0], z[mm, 1], s=16, c=col, label=f"{key} (n={mm.sum()})", lw=0.3, edgecolor="white")
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title="Coloured by the structural cluster of the ligand's chains (>= 80%)")
    ax.legend(fontsize=8)
    n7, n4 = (sm == "cl7").sum(), (sm == "cl4").sum()
    fig.suptitle(f"Ligand chemistry map ({best.fingerprint} Tanimoto, UMAP, 306 ligands): cl7 ligands (n={n7}) and "
                 f"cl4 ligands (n={n4}) occupy different regions of chemical space", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y08_chem_map", section="B. Chemistry clustering", task="Y08",
             shows="UMAP of the ligand Tanimoto distance coloured by chemistry cluster and by structural cluster",
             takeaway=f"{n7} ligands sit only in cl7 chains and {n4} only in cl4 chains; see Y25 for the formal test.",
             source="results/ext/chem_labels_best.csv, results/ext/chem_umap.csv", script="scripts/ext_chem_zoo.py")
    write_numbers("Y08", dict(n_ligands=len(lig), n_runs=len(df), n_eligible=int(elig.sum()), best_run=best.run_id,
                              best_k=int(best.n_clusters), best_sil=float(best.silhouette),
                              best_stab=float(best.stability_ari), best_largest=float(best.largest_share),
                              balanced_run=bal.run_id, balanced_k=int(bal.n_clusters),
                              balanced_sil=float(bal.silhouette), balanced_stab=float(bal.stability_ari),
                              ecfp4_vs_frozen_maxdiff=ecfp_check, n_lig_cl7=int(n7), n_lig_cl4=int(n4)))
    print("ECFP4 vs frozen loi matrix max diff", ecfp_check)
    print(e.head(15)[["rank", "run_id", "n_clusters", "silhouette", "stability_ari", "largest_share", "score"]].to_string())
    print("balanced:", bal.run_id, bal.n_clusters, bal.silhouette, bal.stability_ari)


if __name__ == "__main__":
    main()
