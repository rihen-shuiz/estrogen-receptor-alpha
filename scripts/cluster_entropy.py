"""Tracker R01: Shannon entropy (bits) of the cluster-size distribution at each cut.

Input: results/entry_clusters.csv (S08; keep-all run, one row per eligible chain).
entropy_bits counts chains. entropy_bits_entry_weighted gives every entry a total weight of 1, split evenly over
its chains, so entries with many copies in the asymmetric unit do not dominate. normalised = entropy / log2(n_clusters)
(1 = clusters of equal size; blank when there is a single cluster). Low entropy = a redundant dataset.
Output: results/entropy.csv.
Run from the repo root: python scripts/cluster_entropy.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CUTS = ["1.0", "1.5", "2.0", "2.5"]


def entropy_bits(weights):
    p = np.asarray(weights, dtype=float)
    p = p[p > 0] / p.sum()
    return float(-(p * np.log2(p)).sum())


def main():
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    entry_weight = 1.0 / clusters.groupby("entry_id").chain_key.transform("size")
    rows = []
    for cut in CUTS:
        labels = clusters[f"cut_{cut}"]
        sizes = labels.value_counts()
        h = entropy_bits(sizes)
        rows.append({
            "cut": cut, "n_chains": len(clusters), "n_entries": clusters.entry_id.nunique(),
            "n_clusters": int(sizes.size), "largest_cluster_chains": int(sizes.max()),
            "entropy_bits": h,
            "entropy_bits_entry_weighted": entropy_bits(entry_weight.groupby(labels).sum()),
            "max_entropy_bits": float(np.log2(sizes.size)),
            "normalised": h / np.log2(sizes.size) if sizes.size > 1 else np.nan,
        })
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/entropy.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
