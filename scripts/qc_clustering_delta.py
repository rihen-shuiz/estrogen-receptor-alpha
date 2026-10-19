"""Tracker S13: clustering of the selected ESR1 set with and without the R-free failures.

Failing entries: r_free above qc_r_free (0.30, analysis_settings.json). The script checks that this rule
matches the qc_rfree flag in results/entries.csv before using it.
Method: average linkage on the submatrix of passing entries (the saved S04 matrix is only sliced, not
recomputed), same cuts and method as S08 via esr_analysis.cluster_structures.structural_groups.
Per cut, compares the full set with the passing set: entries, clusters, top-3 share, and full-set clusters
that have no passing member ("lost" clusters). Failing entries stay in the files; they are only left out of
the passing-set clustering.
Output: results/qc_clustering_delta.csv.
Run from the repo root: PYTHONPATH=. python scripts/qc_clustering_delta.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
CUTS = [1.0, 1.5, 2.0, 2.5]


def top3_fraction(labels):
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    return float(sizes.head(3).sum() / len(labels))


def main():
    settings = json.loads(SETTINGS.read_text())
    primary = settings["primary_structural_threshold_A"]
    index = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    entries = pd.read_csv(ROOT / "results/entries.csv", dtype=str, keep_default_na=False)
    merged = index.merge(entries[["entry_id", "r_free", "qc_rfree"]], on="entry_id", how="left", validate="one_to_one")
    r_free = merged.r_free.astype(float)
    fails = r_free > settings["qc_r_free"]
    if not (fails == merged.qc_rfree.eq("fail")).all():
        raise ValueError("R-free rule disagrees with the qc_rfree flag in results/entries.csv")
    passes = ~fails

    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    pass_positions = np.flatnonzero(passes.to_numpy())
    pass_keys = merged.chain_key[passes].tolist()
    sub = matrix[np.ix_(pass_positions, pass_positions)]
    pass_assign, _ = structural_groups(pass_keys, sub, CUTS, primary)

    full_labels = pd.read_csv(ROOT / "results/entry_clusters.csv")
    rows = []
    for cut in CUTS:
        full = full_labels[f"cut_{cut:.1f}"].to_numpy()
        reduced = pass_assign[f"cluster_{cut:.1f}A"].to_numpy()
        lost = len(set(full) - set(full[passes.to_numpy()]))
        rows.append({
            "cut": cut,
            "n_entries_full": len(full),
            "n_entries_pass": int(passes.sum()),
            "n_clusters_full": int(np.unique(full).size),
            "n_clusters_pass": int(np.unique(reduced).size),
            "top3_fraction_full": top3_fraction(full),
            "top3_fraction_pass": top3_fraction(reduced),
            "clusters_without_pass_members": lost,
        })
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/qc_clustering_delta.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
