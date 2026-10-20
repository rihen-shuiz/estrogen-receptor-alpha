"""Tracker E10: ESR1 positions consistently displaced between matched ESR1 and ESR2 chains.

Input: results/displacement_matrix.csv (E07; rows = ESR1 positions, columns = 135 matched chain pairs, blank =
unobserved), results/numbering_map_esr1_esr2.csv (E03, residue identity). For every position: pairs observed,
pairs with displacement > THRESHOLD_A, fraction_pairs = that count / pairs observed, and median displacement.
consistent = fraction_pairs > 0.5 with at least MIN_OBSERVED observed pairs (positions observed in few pairs are kept
with consistent False and the count shown). All positions are written, sorted by position.
Output: results/consistent_positions.csv.
Run from the repo root: python scripts/consistent_positions.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
THRESHOLD_A = 1.5
MIN_OBSERVED = 20


def main():
    wide = pd.read_csv(ROOT / "results/displacement_matrix.csv", index_col="esr1_pos")
    numbering = pd.read_csv(ROOT / "results/numbering_map_esr1_esr2.csv").set_index("esr1_pos")
    core = {p for start, end in CORE["core_ranges"] for p in range(start, end + 1)}
    h12 = range(SETTINGS["h12_range"][0], SETTINGS["h12_range"][1] + 1)
    observed = wide.notna().sum(axis=1)
    above = (wide > THRESHOLD_A).sum(axis=1)
    table = pd.DataFrame({
        "position": wide.index, "n_pairs_observed": observed.to_numpy(), "n_pairs_above": above.to_numpy(),
        "fraction_pairs": (above / observed.where(observed > 0)).to_numpy(),
        "median_displacement_A": wide.median(axis=1).to_numpy(), "threshold_A": THRESHOLD_A,
    })
    table["esr2_pos"] = table.position.map(numbering.esr2_pos)
    table["esr1_aa"] = table.position.map(numbering.esr1_aa)
    table["esr2_aa"] = table.position.map(numbering.esr2_aa)
    table["same_residue"] = table.position.map(numbering.same_residue)
    table["in_core"] = table.position.isin(core)
    table["in_h12"] = table.position.isin(h12)
    table["consistent"] = (table.fraction_pairs > 0.5) & (table.n_pairs_observed >= MIN_OBSERVED)
    table.to_csv(ROOT / "results/consistent_positions.csv", index=False)
    hits = table[table.consistent]
    print(f"E10: {len(table)} positions; consistent (> {THRESHOLD_A} A in > 50% of >= {MIN_OBSERVED} observed pairs): "
          f"{len(hits)}; in core {int(hits.in_core.sum())}, in H12 {int(hits.in_h12.sum())}, "
          f"different residue {int((hits.same_residue == False).sum())}")  # noqa: E712
    print(hits[["position", "esr2_pos", "esr1_aa", "esr2_aa", "n_pairs_observed", "fraction_pairs",
                "median_displacement_A", "in_core", "in_h12"]].to_string(index=False))


if __name__ == "__main__":
    main()
