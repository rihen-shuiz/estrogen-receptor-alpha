"""Rebuild tracker C07 (K-means K=3..7 on all ligands) and C10 (same on non-additive ligands).

Inputs: results/rebuild/fp_ecfp4.npy + fp_index.csv (C02 rebuild); for C10 the additive flag
is_additive from results/ligands_ranked.csv (tracker Q16 output, the input the C10 card names; Q16 itself
is not rebuilt here). A second C10 variant uses only analysis_settings.json solvent_additive_codes.
Method: sklearn KMeans(n_clusters=K, random_state=42, n_init=20) on the 2048-bit ECFP4 matrix (float),
silhouette_score (Euclidean) and inertia per K; labels written 1-based as cluster_k<K>;
chem_cluster = labels of the K with the highest silhouette.
Outputs: results/rebuild/kmeans_scores.csv, ligand_clusters.csv (C07);
results/rebuild/kmeans_scores_noadd.csv, ligand_clusters_noadd.csv (C10, Q16 is_additive);
results/rebuild/kmeans_scores_noadd_settings.csv, ligand_clusters_noadd_settings.csv (settings list only).
Run from the repo root: PYTHONPATH=. python scripts/rebuild_kmeans.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from rebuild_common import OUT, ROOT, SEED, SETTINGS, read_csv

KS = range(3, 8)


def cluster(ids, X, suffix):
    scores, labels = [], {}
    for k in KS:
        model = KMeans(n_clusters=k, random_state=SEED, n_init=20).fit(X)
        labels[k] = model.labels_ + 1
        scores.append({"k": k, "silhouette": silhouette_score(X, model.labels_), "inertia": model.inertia_,
                       "status": "tested", "reason": ""})
    scores = pd.DataFrame(scores)
    best = int(scores.loc[scores.silhouette.idxmax(), "k"])
    table = pd.DataFrame({"ligand_id": ids, "chem_cluster": labels[best],
                          **{f"cluster_k{k}": labels[k] for k in KS}})
    scores.to_csv(OUT / f"kmeans_scores{suffix}.csv", index=False)
    table.to_csv(OUT / f"ligand_clusters{suffix}.csv", index=False)
    print(f"K-means{suffix or ' (all)'}: n={len(ids)}, best K={best}")


def main():
    ecfp = np.load(OUT / "fp_ecfp4.npy").astype(float)
    index = read_csv(OUT / "fp_index.csv")
    ids = index.ligand_id.astype(str)
    cluster(ids.values, ecfp, "")
    ranked = read_csv(ROOT / "results/ligands_ranked.csv")
    additive = set(ranked.loc[ranked.is_additive.astype(str).str.lower() == "true", "ligand_id"].astype(str))
    keep = ~ids.isin(additive)
    cluster(ids[keep].values, ecfp[keep.values], "_noadd")
    keep = ~ids.isin(set(SETTINGS["solvent_additive_codes"]))
    cluster(ids[keep].values, ecfp[keep.values], "_noadd_settings")


if __name__ == "__main__":
    main()
