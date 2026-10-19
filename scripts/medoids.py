"""Tracker S12: medoid entry per cluster at the primary distance cut.

Medoid = the member with the smallest mean RMSD to its cluster mates (ties: first in matrix order).
A singleton cluster is its own medoid, with no mean (blank mean_rmsd_to_mates_A).
Input: results/rmsd_entry.npy (S04), results/rmsd_entry_index.csv (S05), results/entry_clusters.csv (S08),
primary_structural_threshold_A from analysis_settings.json.
Output: results/medoids.csv (cluster, entry_id, chain_key, cluster_size, mean_rmsd_to_mates_A).
Run from the repo root: python scripts/medoids.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"


def main():
    cut = json.loads(SETTINGS.read_text())["primary_structural_threshold_A"]
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    index = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    labels = pd.read_csv(ROOT / "results/entry_clusters.csv")[f"cut_{cut:.1f}"].to_numpy()
    rows = []
    for cluster in sorted(set(labels)):
        members = np.flatnonzero(labels == cluster)
        if members.size == 1:
            best, mean_to_mates = members[0], np.nan
        else:
            block = matrix[np.ix_(members, members)]
            means = (block.sum(axis=1)) / (members.size - 1)  # diagonal is zero, so this excludes self
            pick = int(np.argmin(means))
            best, mean_to_mates = members[pick], float(means[pick])
        rows.append({
            "cluster": int(cluster),
            "entry_id": index.entry_id[best],
            "chain_key": index.chain_key[best],
            "cluster_size": int(members.size),
            "mean_rmsd_to_mates_A": mean_to_mates,
        })
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/medoids.csv", index=False)
    print(f"S12: {len(table)} medoids at cut {cut} A; sizes {table.cluster_size.tolist()}")


if __name__ == "__main__":
    main()
