"""Tracker S09: per-cut cluster statistics for the selected-chain clustering.

Input: results/rmsd_entry.npy (S04), results/entry_clusters.csv (S08, same row order as the matrix).
Rows are chains (keep-all run: every eligible chain, an entry with several chains appears several times).
Per cut: number of clusters, share of chains in the three largest clusters, and the median RMSD of
chain pairs inside one cluster and of pairs in different clusters (angstrom, over the S04 matrix).
Output: results/cluster_stats.csv.
Run from the repo root: python scripts/cluster_stats.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CUTS = ["1.0", "1.5", "2.0", "2.5"]


def cut_stats(labels, matrix, upper):
    same = labels[upper[0]] == labels[upper[1]]
    distances = matrix[upper]
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    return {
        "n_clusters": int(sizes.size),
        "top3_fraction": float(sizes.head(3).sum() / labels.size),
        "median_within_A": float(np.median(distances[same])) if same.any() else np.nan,
        "median_between_A": float(np.median(distances[~same])) if (~same).any() else np.nan,
    }


def main():
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    matrix_keys = (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().splitlines()
    if clusters.chain_key.tolist() != matrix_keys:
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")
    upper = np.triu_indices(matrix.shape[0], k=1)
    rows = []
    for cut in CUTS:
        labels = clusters[f"cut_{cut}"].to_numpy()
        rows.append({"cut": cut, "n_chains": labels.size, "n_entries": clusters.entry_id.nunique(),
                     **cut_stats(labels, matrix, upper)})
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/cluster_stats.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
