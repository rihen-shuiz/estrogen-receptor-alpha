"""Tracker R04: 90% saturation N at each cut, from the R02 curve.

Input: results/diversity_curve.csv (R02; 788 chains, 100 orderings). n_90 = smallest N whose mean number of
clusters seen reaches 90% of the cluster total at that cut. Also the share of chains needed (n_90 / n_chains).
Output: results/saturation.csv.
Run from the repo root: python scripts/saturation.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    curve = pd.read_csv(ROOT / "results/diversity_curve.csv", dtype={"cut": str})
    rows = []
    for cut, part in curve.groupby("cut", sort=False):
        total = int(part.n_clusters_total.iloc[0])
        reached = part[part.mean_clusters >= 0.9 * total]
        n_90 = int(reached.n.iloc[0])
        if "n_90" in part and int(part.n_90.iloc[0]) != n_90:
            raise ValueError(f"cut {cut}: n_90 disagrees with the value written by R02")
        rows.append({"cut": cut, "n_90": n_90, "n_clusters_total": total, "threshold_90pct": 0.9 * total,
                     "mean_at_n90": float(reached.mean_clusters.iloc[0]), "n_chains": int(part.n.max()),
                     "fraction_of_chains": n_90 / int(part.n.max()),
                     "seed": int(part.seed.iloc[0]), "n_orderings": int(part.n_orderings.iloc[0])})
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/saturation.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
