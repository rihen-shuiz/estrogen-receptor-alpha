"""Tracker K09: alignment-core sensitivity, old core (84 residues) vs data-chosen core v3 (123 residues, primary).

Old-core matrix: pairwise C-alpha RMSD from the old-core fit stored by S03
(results/alignment/aligned_coordinates_old_core.npz), computed exactly like S04 (esr_analysis.calculate_rmsd.
compute_distances, same min_pairwise_common_ca) and saved as a NEW file results/rmsd_old_core.npy. The primary matrix
results/rmsd_entry.npy is only read. Row order must equal results/alignment/rmsd_entry_index.txt.
Clustering: S08 method (esr_analysis.cluster_structures.structural_groups) at every S08 cut.
Comparison with the primary result: clusters per cut, adjusted Rand index against the S08 labels, and Spearman
correlation of all pairwise RMSDs (upper triangle).
Output: results/core_sensitivity.csv (one row per core; the primary row is the reference) and results/rmsd_old_core.npy.
Run from the repo root: PYTHONPATH=. python scripts/core_sensitivity.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score

from esr_analysis._task_io import aligned_caches
from esr_analysis.calculate_rmsd import compute_distances
from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE_V3 = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
CUTS = [1.0, 1.5, 2.0, 2.5]


def core_size(ranges):
    return sum(end - start + 1 for start, end in ranges)


def main():
    keys = (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().split()
    primary = np.load(ROOT / "results/rmsd_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    if clusters.chain_key.tolist() != keys:
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")

    caches = aligned_caches(ROOT / "results/alignment/aligned_coordinates_old_core.npz", keys)
    kept, old, _, excluded = compute_distances(keys, caches, SETTINGS["min_pairwise_common_ca"])
    if kept != keys or np.isnan(old).any():
        raise ValueError(f"old-core matrix does not cover the primary chains ({len(kept)} kept, {len(excluded)} excluded)")
    np.save(ROOT / "results/rmsd_old_core.npy", old)

    upper = np.triu_indices(len(keys), k=1)
    rho = float(spearmanr(primary[upper], old[upper])[0])
    assignments, _ = structural_groups(keys, old, CUTS, SETTINGS["primary_structural_threshold_A"])
    rows = []
    for core, size, labels, spearman in [
        (CORE_V3["name"] + " (primary)", core_size(CORE_V3["core_ranges"]), {c: clusters[f"cut_{c:.1f}"] for c in CUTS}, 1.0),
        ("old_core", core_size(SETTINGS["esr1_core_ranges"]), {c: assignments[f"cluster_{c:.1f}A"] for c in CUTS}, rho),
    ]:
        row = {"core": core, "core_residues": size, "n_chains": len(keys), "spearman_pairwise": spearman,
               "median_abs_diff_A": float(np.median(np.abs(primary[upper] - old[upper]))) if core == "old_core" else 0.0,
               "max_abs_diff_A": float(np.max(np.abs(primary[upper] - old[upper]))) if core == "old_core" else 0.0}
        for cut in CUTS:
            tag = f"{cut:.1f}".replace(".", "_")
            row[f"n_clusters_{tag}"] = int(pd.Series(labels[cut]).nunique())
            row[f"ari_{tag}"] = float(adjusted_rand_score(clusters[f"cut_{cut:.1f}"], labels[cut]))
        row["ari_vs_primary"] = row["ari_1_5"]
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/core_sensitivity.csv", index=False)
    print(table.T.to_string())


if __name__ == "__main__":
    main()
