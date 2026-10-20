"""Tracker S13 (keep-all): QC-tier sensitivity of the clustering, all chains vs tier A only vs tiers A+B.

Replaces the old "drop R-free > 0.30" comparison. The main result stays on ALL chains; this table only shows
whether the conclusion depends on the poorer chains. Nothing is removed from any file.
Inputs: results/rmsd_entry.npy (S04, sliced, never recomputed), results/entry_clusters.csv (S08 all-chain labels and
the K05 qc_label copied per chain).
Method: same as S08, esr_analysis.cluster_structures.structural_groups (average linkage, distance cuts) on the
submatrix of each subset.
Per cut and subset: chains, entries, clusters, top-3 share, all-chain clusters with no member in the subset, and the
adjusted Rand index between the subset labelling and the all-chain labelling of the same chains.
Output: results/qc_clustering_delta.csv.
Run from the repo root: PYTHONPATH=. python scripts/qc_clustering_delta.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
CUTS = [1.0, 1.5, 2.0, 2.5]
SUBSETS = {"all": None, "tier_A": {"A"}, "tier_A+B": {"A", "B"}}


def top3_fraction(labels):
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    return float(sizes.head(3).sum() / len(labels))


def main():
    primary = json.loads(SETTINGS.read_text())["primary_structural_threshold_A"]
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    if clusters.chain_key.tolist() != (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().splitlines():
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")
    rows = []
    for subset, labels_kept in SUBSETS.items():
        mask = np.ones(len(clusters), bool) if labels_kept is None else clusters.qc_label.isin(labels_kept).to_numpy()
        keys = clusters.chain_key[mask].tolist()
        assignments, _ = structural_groups(keys, matrix[np.ix_(mask, mask)], CUTS, primary)
        for cut in CUTS:
            sub = assignments[f"cluster_{cut:.1f}A"].to_numpy()
            full = clusters[f"cut_{cut:.1f}"].to_numpy()
            rows.append({
                "cut": f"{cut:.1f}", "subset": subset, "n_chains": int(mask.sum()),
                "n_entries": int(clusters.entry_id[mask].nunique()),
                "n_clusters": int(len(np.unique(sub))), "top3_fraction": top3_fraction(sub),
                "clusters_all": int(len(np.unique(full))),
                "clusters_without_members": int(len(set(full) - set(full[mask]))),
                "ari_vs_all": float(adjusted_rand_score(full[mask], sub)),
            })
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/qc_clustering_delta.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
