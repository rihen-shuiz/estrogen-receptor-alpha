"""Tracker R02: cumulative diversity curve over all eligible chains (100 random orderings).

Input: results/entry_clusters.csv (S08; keep-all run, one row per eligible chain, 788 chains).
For each cut, 100 random orderings of the chains (one generator, random_seed from analysis_settings.json) give the
number of distinct clusters seen after N chains, N = 1..total. Mean and SD (ddof=1) over the orderings.
n_90 (smallest N whose mean reaches 90% of the cluster total) is repeated on every row of its cut for reading.
Helpers: scripts/provisional_common.py (clusters_seen; same counting as the provisional P02 curve).
Output: results/diversity_curve.csv.
Run from the repo root: python scripts/diversity_curve.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provisional_common import N_ORDERINGS, clusters_seen  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
SEED = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["random_seed"]
CUTS = ["1.0", "1.5", "2.0", "2.5"]


def main():
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    rng = np.random.default_rng(SEED)
    orders = np.array([rng.permutation(len(clusters)) for _ in range(N_ORDERINGS)])
    rows, summary = [], []
    for cut in CUTS:
        labels = clusters[f"cut_{cut}"].to_numpy()
        curves = clusters_seen(labels, orders)
        mean, sd = curves.mean(axis=0), curves.std(axis=0, ddof=1)
        total = len(np.unique(labels))
        n_90 = int(np.argmax(mean >= 0.9 * total)) + 1
        summary.append(f"{cut} A: {total} clusters, n_90 {n_90}")
        for i in range(len(labels)):
            rows.append({"cut": cut, "n": i + 1, "mean_clusters": mean[i], "sd_clusters": sd[i],
                         "n_clusters_total": total, "n_90": n_90, "seed": SEED, "n_orderings": N_ORDERINGS})
    pd.DataFrame(rows).to_csv(ROOT / "results/diversity_curve.csv", index=False)
    print(f"R02: {len(clusters)} chains, {len(rows)} rows; " + "; ".join(summary))


if __name__ == "__main__":
    main()
