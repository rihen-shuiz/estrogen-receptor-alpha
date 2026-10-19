"""Tracker S16 (keep-all): H12 position per chain, measured after the S03 core fit, for every eligible chain.

For each chain in results/alignment/aligned_coordinates.npz (S03 primary core, 788 chains), take the aligned
C-alpha atoms of H12 (settings h12_range, residues 538-547) and compare them with the same residues of the
reference chain in the saved reference frame. The core fit is not repeated on H12.
Labels: n_h12_ca = shared H12 C-alpha. 0 shared -> h12_label "not_modelled" and h12_rmsd blank;
fewer than the full H12 range -> "partial"; otherwise "complete". Nothing is filtered: large values stay,
next to n_h12_ca. h12_mean_b and h12_b_group come from results/helix_bfactor_by_chain.csv (K03).
Output: results/h12_metric.csv
Run from the repo root: python scripts/h12_metric.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = ROOT / "esr_analysis/resources/config/analysis_settings.json"
REFERENCE = ROOT / "esr_analysis/resources/frozen/reference.json"
FITS = ROOT / "results/alignment/aligned_coordinates.npz"
BFACTORS = ROOT / "results/helix_bfactor_by_chain.csv"
OUT = ROOT / "results/h12_metric.csv"


def main():
    low, high = json.loads(SETTINGS.read_text())["h12_range"]
    h12 = list(range(low, high + 1))
    fits = np.load(FITS, allow_pickle=False)
    keys = fits["chain_keys"].astype(str).tolist()
    position_column = {int(p): i for i, p in enumerate(fits["positions"])}
    h12_cols = [position_column[p] for p in h12]
    coords = fits["coords"]
    row_of = {key: i for i, key in enumerate(keys)}
    ref_key = json.loads(REFERENCE.read_text())["chain_key"]
    ref = coords[row_of[ref_key]][h12_cols]

    bfactors = pd.read_csv(BFACTORS, keep_default_na=False, na_values=[""]).set_index("chain_key")
    rows = []
    for chain in keys:
        mobile = coords[row_of[chain]][h12_cols]
        shared = np.isfinite(mobile).all(axis=1) & np.isfinite(ref).all(axis=1)
        n_shared = int(shared.sum())
        if n_shared != int(bfactors.loc[chain, "h12_n_ca"] or 0):
            raise ValueError(f"H12 CA count disagrees with K03 for {chain}")
        value = float(np.sqrt(np.mean(np.sum((mobile[shared] - ref[shared]) ** 2, axis=1)))) if n_shared else np.nan
        label = "not_modelled" if n_shared == 0 else ("partial" if n_shared < len(h12) else "complete")
        rows.append({"entry_id": chain.split("|")[0], "chain_key": chain, "n_h12_ca": n_shared, "h12_rmsd": value,
                     "h12_label": label, "h12_mean_b": bfactors.loc[chain, "h12_mean_b"],
                     "h12_b_group": bfactors.loc[chain, "h12_b_group"]})
    table = pd.DataFrame(rows)
    table.to_csv(OUT, index=False)
    counts = table.h12_label.value_counts().to_dict()
    above_5 = int((table.h12_rmsd > 5).sum())
    print(f"S16: {len(table)} chains; labels {counts}; h12_rmsd median {table.h12_rmsd.median():.3f} "
          f"max {table.h12_rmsd.max():.3f} A; above 5 A: {above_5}")


if __name__ == "__main__":
    main()
