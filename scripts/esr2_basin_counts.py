"""Tracker X05: ESR2 basin sizes with a confidence caveat.

Input: results/esr2_clusters.csv (X03). Basin = X03 cluster at the primary cut. Per basin: chains, entries, the
entries themselves, and interpretable = at least MIN_ENTRIES distinct entries (chains of one entry are one crystal,
so they do not add independent evidence). Basins below that are kept and marked, not removed.
Output: results/esr2_basin_counts.csv.
Run from the repo root: python scripts/esr2_basin_counts.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CUT = f'cut_{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
MIN_ENTRIES = 3


def main():
    clusters = pd.read_csv(ROOT / "results/esr2_clusters.csv")
    table = clusters.groupby(CUT).agg(n_chains=("chain_key", "size"), n_entries=("entry_id", "nunique"),
                                      entries=("entry_id", lambda s: ";".join(sorted(set(s)))))
    table.index.name = "basin"
    table = table.reset_index().sort_values("n_chains", ascending=False)
    table["interpretable"] = table.n_entries >= MIN_ENTRIES
    table["caveat"] = table.interpretable.map({True: "", False: f"fewer than {MIN_ENTRIES} entries: not interpretable"})
    table["cut"] = CUT.removeprefix("cut_")
    table.to_csv(ROOT / "results/esr2_basin_counts.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
