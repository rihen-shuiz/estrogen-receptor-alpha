"""Tracker K03: normalised C-alpha B-factor per helix, for every eligible ESR1 chain.

Chains: results/lbd_chains.csv, eligible_historical True (788). No chain is dropped, high-B chains included.
Per chain: each LBD C-alpha B is divided by that chain's own median C-alpha B (normalised B), then averaged per helix
(helix ranges from results/helix_table.csv, K02). Raw means are kept beside them. A helix with no observed
C-alpha in the chain is blank, not zero.
Per chain also: core mean B (current settings core ranges; K04 has not chosen the new core yet), H12 mean B
(settings h12_range 538-547), their raw and normalised ratios, b_group (tertile of mean_b_chain) and h12_b_group
(tertile of h12_mean_b). Tertiles use rank order so each group has about the same number of chains.
Residues with an invalid B (NaN or negative) are left blank in the B means, never dropped.
Output: results/helix_bfactor_by_chain.csv
Run from the repo root: python scripts/helix_bfactors.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
LBD = (SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"])
CORE_RANGES = SETTINGS["esr1_core_ranges"]
H12 = tuple(SETTINGS["h12_range"])
OUT = ROOT / "results/helix_bfactor_by_chain.csv"


def tertile(values):
    ranks = values.rank(method="first")
    labels = pd.qcut(ranks, 3, labels=["low", "mid", "high"])
    return labels.astype(object).where(values.notna(), "")


def main():
    chains = pd.read_csv(ROOT / "results/lbd_chains.csv", dtype={"entry_id": str}, keep_default_na=False, na_values=[""])
    eligible = chains[(chains.identity == "esr1") & chains.eligible_historical.astype(str).str.lower().eq("true")]
    keys = sorted(eligible.chain_key)
    helices = pd.read_csv(ROOT / "results/helix_table.csv", keep_default_na=False)

    residues = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", dtype={"chain_key": str}, keep_default_na=False,
                           na_values=[""], low_memory=False, usecols=["chain_key", "esr1_position", "has_ca", "B_factor_A2"])
    residues = residues[residues.chain_key.isin(set(keys)) & residues.has_ca.astype(str).str.lower().eq("true")]
    residues = residues.dropna(subset=["esr1_position"]).copy()
    residues["position"] = residues.esr1_position.astype(float).astype(int)
    residues = residues[residues.position.between(*LBD)]
    residues["b"] = residues.B_factor_A2.astype(float)
    residues.loc[~(np.isfinite(residues.b) & (residues.b >= 0)), "b"] = np.nan

    core_positions = [p for start, end in CORE_RANGES for p in range(start, end + 1)]
    h12_positions = list(range(H12[0], H12[1] + 1))
    rows = []
    for key in keys:
        block = residues[residues.chain_key == key].set_index("position")
        valid = block.b.dropna()
        median_b = float(valid.median()) if len(valid) else np.nan
        norm = block.b / median_b
        record = {"chain_key": key, "entry_id": key.split("|")[0], "n_lbd_ca_valid_b": int(len(valid)),
                  "median_ca_b_chain": median_b, "mean_b_chain": float(valid.mean()) if len(valid) else np.nan}

        def region(positions):
            raw = block.b.reindex(positions).dropna()
            normed = norm.reindex(positions).dropna()
            return (float(raw.mean()) if len(raw) else np.nan, float(normed.mean()) if len(normed) else np.nan, int(len(raw)))

        core_raw, core_norm, core_n = region(core_positions)
        h12_raw, h12_norm, h12_n = region(h12_positions)
        record.update({"core_mean_b": core_raw, "core_norm_b": core_norm, "core_n_ca": core_n,
                       "h12_mean_b": h12_raw, "h12_norm_b": h12_norm, "h12_n_ca": h12_n,
                       "h12_to_core_ratio": h12_raw / core_raw if core_n and h12_n else np.nan,
                       "h12_to_core_norm_ratio": h12_norm / core_norm if core_n and h12_n else np.nan})
        for helix in helices.itertuples():
            positions = list(range(int(helix.start), int(helix.end) + 1))
            raw, normed, n = region(positions)
            record[f"{helix.helix}_raw_mean_b"] = raw
            record[f"{helix.helix}_norm_mean_b"] = normed
            record[f"{helix.helix}_n_ca"] = n
        rows.append(record)
    table = pd.DataFrame(rows)
    table["b_group"] = tertile(table.mean_b_chain)
    table["h12_b_group"] = tertile(table.h12_mean_b)
    table.to_csv(OUT, index=False)
    cut = table.mean_b_chain.quantile([1/3, 2/3]).round(2).tolist()
    print(f"K03: {len(table)} eligible chains; b_group counts {table.b_group.value_counts().to_dict()}; "
          f"mean_b_chain tertile cut points near {cut}")
    print(f"K03: chains with no H12 CA: {int((table.h12_n_ca == 0).sum())}; "
          f"chains missing at least one helix: {int(table[[c for c in table.columns if c.endswith('_n_ca') and c not in ('core_n_ca','h12_n_ca')]].eq(0).any(axis=1).sum())}")


if __name__ == "__main__":
    main()
