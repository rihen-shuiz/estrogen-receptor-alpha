"""Tracker R06: entropy and 90% saturation N, full chain set vs chains passing R-free <= qc_r_free (0.30).

Subsets: all 788 chains (S08 labels) and the chains of entries with R-free <= qc_r_free (results/qc_tiers.csv, K05;
an entry with no R-free stays out of the passing set and is counted). The passing set is reclustered as in S13:
average linkage on the slice of the saved S04 matrix (not recomputed), same cuts. tier_A (K05 qc_label A) is added
as a stricter third row per cut for comparison with S13. Nothing is removed from any file.
Entropy: Shannon bits of the cluster sizes (as R01). n_90: smallest N whose mean clusters seen over 100 random
orderings reaches 90% of that subset's total (as R02/R04, random_seed from analysis_settings.json).
Output: results/entropy_qc_compare.csv.
Run from the repo root: PYTHONPATH=. python scripts/entropy_qc_compare.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cluster_entropy import entropy_bits  # noqa: E402
from provisional_common import N_ORDERINGS, clusters_seen  # noqa: E402

from esr_analysis.cluster_structures import structural_groups  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CUTS = [1.0, 1.5, 2.0, 2.5]


def n_90(labels, seed):
    rng = np.random.default_rng(seed)
    orders = np.array([rng.permutation(len(labels)) for _ in range(N_ORDERINGS)])
    mean = clusters_seen(labels, orders).mean(axis=0)
    return int(np.argmax(mean >= 0.9 * len(np.unique(labels)))) + 1


def main():
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    if clusters.chain_key.tolist() != (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().splitlines():
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")
    rfree = clusters.entry_id.map(pd.read_csv(ROOT / "results/qc_tiers.csv").set_index("entry_id").r_free)
    subsets = {
        "all": np.ones(len(clusters), bool),
        f"rfree<={SETTINGS['qc_r_free']:.2f}": (rfree <= SETTINGS["qc_r_free"]).to_numpy(),
        "tier_A": clusters.qc_label.eq("A").to_numpy(),
    }
    rows = []
    for name, mask in subsets.items():
        keys = clusters.chain_key[mask].tolist()
        if name == "all":
            labels = {cut: clusters[f"cut_{cut:.1f}"].to_numpy() for cut in CUTS}
        else:
            assignments, _ = structural_groups(keys, matrix[np.ix_(mask, mask)], CUTS,
                                               SETTINGS["primary_structural_threshold_A"])
            labels = {cut: assignments[f"cluster_{cut:.1f}A"].to_numpy() for cut in CUTS}
        for cut in CUTS:
            sizes = pd.Series(labels[cut]).value_counts()
            rows.append({"cut": f"{cut:.1f}", "subset": name, "n_chains": int(mask.sum()),
                         "n_entries": int(clusters.entry_id[mask].nunique()),
                         "n_chains_left_out": int((~mask).sum()), "n_clusters": int(sizes.size),
                         "entropy_bits": entropy_bits(sizes), "n_90": n_90(labels[cut], SETTINGS["random_seed"]),
                         "n_90_fraction": n_90(labels[cut], SETTINGS["random_seed"]) / int(mask.sum())})
    table = pd.DataFrame(rows)
    full = table[table.subset == "all"].set_index("cut")
    table["entropy_delta_vs_all"] = table.entropy_bits - table.cut.map(full.entropy_bits)
    table["n_90_fraction_delta_vs_all"] = table.n_90_fraction - table.cut.map(full.n_90_fraction)
    table.to_csv(ROOT / "results/entropy_qc_compare.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
