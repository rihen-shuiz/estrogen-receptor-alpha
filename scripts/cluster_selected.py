"""Tracker S08: average-linkage clusters of the selected-chain RMSD matrix at four distance cuts.

Input: results/rmsd_entry.npy (S04) and results/rmsd_entry_index.csv (S05, same row order).
Method: esr_analysis/cluster_structures.structural_groups (SciPy average linkage, distance criterion).
Outputs: results/entry_clusters.csv (entry_id, chain_key, cut_1.0, cut_1.5, cut_2.0, cut_2.5)
and results/linkage_entry.npy (linkage matrix for the dendrogram).
Run from the repo root: PYTHONPATH=. python scripts/cluster_selected.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
CUTS = [1.0, 1.5, 2.0, 2.5]


def main():
    primary = json.loads(SETTINGS.read_text())["primary_structural_threshold_A"]
    index = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    assignments, linkage_matrix = structural_groups(index.chain_key.tolist(), matrix, CUTS, primary)
    table = pd.DataFrame({"entry_id": index.entry_id, "chain_key": index.chain_key})
    for cut in CUTS:
        table[f"cut_{cut:.1f}"] = assignments[f"cluster_{cut:.1f}A"].to_numpy()
    table.to_csv(ROOT / "results/entry_clusters.csv", index=False)
    np.save(ROOT / "results/linkage_entry.npy", linkage_matrix)
    counts = {f"{cut:.1f}": int(table[f"cut_{cut:.1f}"].nunique()) for cut in CUTS}
    print(f"S08: {len(table)} entries; clusters per cut {counts}; primary cut {primary} A")


if __name__ == "__main__":
    main()
