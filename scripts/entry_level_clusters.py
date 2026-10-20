"""Tracker K07: entry-level view next to the chain-level clustering. Both views are kept; neither replaces the other.

Chains of one entry are not independent crystals. One chain per entry, among the chains in the S04 matrix, by the
stated rule: most mapped LBD C-alpha (results/chain_status.csv, K06), then lowest mean C-alpha B (mean_b_chain,
results/helix_bfactor_by_chain.csv, K03), then chain_key order.
Those rows are sliced from the saved S04 matrix (not recomputed) and clustered by the S08 method at the four cuts.
Comparison per cut: clusters, top-3 share (entry view: share of entries; chain view: share of chains), and the ARI
between the entry-level labels and the chain-level S08 labels of the same picked chains.
Outputs: results/entry_level_clusters.csv (one row per entry) and results/entry_vs_chain_clusters.csv.
Run from the repo root: PYTHONPATH=. python scripts/entry_level_clusters.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from esr_analysis.cluster_structures import structural_groups

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CUTS = [1.0, 1.5, 2.0, 2.5]


def top3(labels):
    return float(pd.Series(labels).value_counts().head(3).sum() / len(labels))


def main():
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    chains = pd.read_csv(ROOT / "results/entry_clusters.csv")
    keys = (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().splitlines()
    if chains.chain_key.tolist() != keys:
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")
    status = pd.read_csv(ROOT / "results/chain_status.csv", usecols=["chain_key", "mapped_lbd_ca"])
    bfac = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv", usecols=["chain_key", "mean_b_chain"])
    ranked = chains.merge(status, on="chain_key", validate="one_to_one").merge(bfac, on="chain_key", validate="one_to_one")
    ranked["row"] = np.arange(len(ranked))
    ranked = ranked.sort_values(["entry_id", "mapped_lbd_ca", "mean_b_chain", "chain_key"],
                                ascending=[True, False, True, True], na_position="last")
    picked = ranked.groupby("entry_id", sort=True).head(1).sort_values("row")
    rows = picked.row.to_numpy()

    assignments, _ = structural_groups(picked.chain_key.tolist(), matrix[np.ix_(rows, rows)], CUTS,
                                       SETTINGS["primary_structural_threshold_A"])
    entry_table = picked[["entry_id", "chain_key", "mapped_lbd_ca", "mean_b_chain"]].reset_index(drop=True)
    entry_table["n_chains_in_entry"] = entry_table.entry_id.map(chains.entry_id.value_counts())
    compare = []
    for cut in CUTS:
        tag = f"{cut:.1f}"
        entry_labels = assignments[f"cluster_{tag}A"].to_numpy()
        entry_table[f"cut_{tag}"] = entry_labels
        entry_table[f"chain_level_cut_{tag}"] = picked[f"cut_{tag}"].to_numpy()
        chain_labels = chains[f"cut_{tag}"].to_numpy()
        compare.append({"cut": tag, "n_chains": len(chains), "n_entries": len(picked),
                        "n_clusters_chain": int(len(np.unique(chain_labels))),
                        "n_clusters_entry": int(len(np.unique(entry_labels))),
                        "top3_chain": top3(chain_labels), "top3_entry": top3(entry_labels),
                        "ari": float(adjusted_rand_score(picked[f"cut_{tag}"], entry_labels))})
    entry_table.to_csv(ROOT / "results/entry_level_clusters.csv", index=False)
    compare = pd.DataFrame(compare)
    compare.to_csv(ROOT / "results/entry_vs_chain_clusters.csv", index=False)
    print(f"K07: {len(picked)} entries (one chain each) from {len(chains)} chains")
    print(compare.to_string(index=False))


if __name__ == "__main__":
    main()
