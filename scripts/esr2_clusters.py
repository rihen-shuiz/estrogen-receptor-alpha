"""Tracker X03: ESR2 clustering, same method as S08 on the X02 ESR2 matrix.

Input: results/rmsd_esr2.npy and results/rmsd_esr2_index.txt (X02, 86 chains, read only), results/esr2_entries.csv
(qc_pass copied for reading only, never used to filter).
Method: esr_analysis.cluster_structures.structural_groups, average linkage, cuts 1.0 / 1.5 / 2.0 / 2.5 A
(primary 1.5 A). Every chain in the matrix gets a label.
Outputs: results/esr2_clusters.csv (entry_id, chain_key, cut_1.0 .. cut_2.5, qc_pass) and results/esr2_linkage.npy.
Run from the repo root: PYTHONPATH=. python scripts/esr2_clusters.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
CUTS = [1.0, 1.5, 2.0, 2.5]


def main():
    primary = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]
    keys = (ROOT / "results/rmsd_esr2_index.txt").read_text().split()
    matrix = np.load(ROOT / "results/rmsd_esr2.npy")
    assignments, linkage_matrix = structural_groups(keys, matrix, CUTS, primary)
    table = pd.DataFrame({"entry_id": [k.split("|")[0] for k in keys], "chain_key": keys})
    for cut in CUTS:
        table[f"cut_{cut:.1f}"] = assignments[f"cluster_{cut:.1f}A"].to_numpy()
    qc = pd.read_csv(ROOT / "results/esr2_entries.csv", usecols=["entry_id", "qc_pass"])
    table = table.merge(qc, on="entry_id", how="left", validate="many_to_one")
    table.to_csv(ROOT / "results/esr2_clusters.csv", index=False)
    np.save(ROOT / "results/esr2_linkage.npy", linkage_matrix)
    counts = {f"{c:.1f}": int(table[f"cut_{c:.1f}"].nunique()) for c in CUTS}
    sizes = table[f"cut_{primary:.1f}"].value_counts().tolist()
    print(f"X03: {len(table)} ESR2 chains, {table.entry_id.nunique()} entries; clusters per cut {counts}; sizes at {primary} A {sizes}")


if __name__ == "__main__":
    main()
