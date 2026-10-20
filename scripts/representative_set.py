"""Tracker S15: minimal representative set covering >= 90% of conformational diversity at the primary cut.

Diversity is counted as R02/R04 count it: distinct structural clusters at primary_structural_threshold_A (S08 labels).
Covering 90% therefore needs ceil(0.9 x n_clusters) clusters, one representative each, taken largest cluster first.
cum_cluster_fraction and cum_chain_fraction show both readings of coverage for every prefix of the list.
Representative per cluster: the S12 medoid when its entry is QC-pass (K05 qc_label A); otherwise the QC-pass chain
with the smallest mean RMSD to its cluster mates (saved S04 matrix, only read); if the cluster has no QC-pass chain,
the medoid is kept and qc_pass is False, so the shortfall is visible.
Outputs: results/representative_set.csv (one row per representative, with a one-line justification).
Run from the repo root: python scripts/representative_set.py
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CUT = f"cut_{SETTINGS['primary_structural_threshold_A']:.1f}"
TARGET = 0.9


def main():
    matrix = np.load(ROOT / "results/rmsd_entry.npy")
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    if clusters.chain_key.tolist() != (ROOT / "results/alignment/rmsd_entry_index.txt").read_text().splitlines():
        raise ValueError("entry_clusters.csv row order does not match the RMSD matrix")
    medoids = pd.read_csv(ROOT / "results/medoids.csv").set_index("cluster")
    sizes = clusters[CUT].value_counts()
    needed = math.ceil(TARGET * len(sizes))
    rows = []
    for rank, (cluster, size) in enumerate(sizes.items(), start=1):
        members = np.flatnonzero(clusters[CUT].to_numpy() == cluster)
        medoid = medoids.loc[cluster]
        medoid_row = int(np.flatnonzero(clusters.chain_key.to_numpy() == medoid.chain_key)[0])
        passing = members[clusters.qc_label.to_numpy()[members] == "A"]
        mean_to_mates = (matrix[np.ix_(members, members)].sum(axis=1) / max(len(members) - 1, 1))
        if clusters.qc_label[medoid_row] == "A":
            pick, why = medoid_row, "S12 medoid, QC tier A"
        elif len(passing):
            pick = int(passing[np.argmin(mean_to_mates[np.isin(members, passing)])])
            why = (f"medoid {medoid.chain_key} is tier {clusters.qc_label[medoid_row]}; most central tier-A chain "
                   f"of {len(passing)} in the cluster")
        else:
            pick, why = medoid_row, f"S12 medoid kept; no tier-A chain in this cluster (medoid tier {clusters.qc_label[medoid_row]})"
        mean_pick = float(mean_to_mates[list(members).index(pick)]) if len(members) > 1 else np.nan
        rows.append({
            "rank": rank, "cluster": int(cluster), "entry_id": clusters.entry_id[pick], "chain_key": clusters.chain_key[pick],
            "qc_label": clusters.qc_label[pick], "qc_pass": clusters.qc_label[pick] == "A",
            "cluster_size_chains": int(size), "cluster_entries": int(clusters.entry_id[members].nunique()),
            "mean_rmsd_to_mates_A": mean_pick, "is_medoid": pick == medoid_row,
            "in_minimal_set": rank <= needed,
            "justification": f"cluster {cluster} at {CUT[4:]} A ({size} chains, rank {rank} by size): {why}",
        })
    table = pd.DataFrame(rows)
    table["cum_cluster_fraction"] = table["rank"] / len(sizes)
    table["cum_chain_fraction"] = table.cluster_size_chains.cumsum() / len(clusters)
    table.to_csv(ROOT / "results/representative_set.csv", index=False)
    minimal = table[table.in_minimal_set]
    print(f"S15: {len(sizes)} clusters at {CUT[4:]} A; 90% needs {needed}; minimal set {len(minimal)} chains, "
          f"QC-pass {int(minimal.qc_pass.sum())}/{len(minimal)}; chain coverage {minimal.cum_chain_fraction.iloc[-1]:.3f}")
    print(table[["rank", "cluster", "chain_key", "qc_label", "cluster_size_chains", "is_medoid", "justification"]].to_string(index=False))


if __name__ == "__main__":
    main()
