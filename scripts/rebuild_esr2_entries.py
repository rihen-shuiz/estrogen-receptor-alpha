"""Rebuild tracker X01: ESR2 entries with LBD chain counts and QC flags (results/esr2_entries.csv).

Inputs: data/pdb_ids_esr2.txt (row order), results/rebuild/entries.csv (Q06 rebuild),
results/rebuild/lbd_chains.csv (Q04 rebuild).
Method: Q04 and Q06 logic already run on every ESR2 entry; this joins them. esr2_chains = number of
chain/model rows with identity esr2; esr2_eligible_chains = those with eligible_historical True.
Output: results/rebuild/esr2_entries.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_esr2_entries.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd

from rebuild_common import OUT, esr2_ids, read_csv


def main():
    entries = read_csv(OUT / "entries.csv").set_index("entry_id")
    chains = read_csv(OUT / "lbd_chains.csv")
    esr2 = chains[chains.identity == "esr2"]
    eligible = esr2.eligible_historical.astype(str).str.lower().eq("true")
    total = esr2.groupby("pdb_id").size()
    good = esr2[eligible].groupby("pdb_id").size()
    ids = esr2_ids()
    table = entries.loc[ids, ["resolution", "r_free", "mean_b", "qc_res", "qc_rfree", "qc_b", "qc_pass", "ligand_ids"]].copy()
    table["esr2_chains"] = total.reindex(ids).fillna(0).astype(int).values
    table["esr2_eligible_chains"] = good.reindex(ids).fillna(0).astype(int).values
    table.index.name = "entry_id"
    table.reset_index().to_csv(OUT / "esr2_entries.csv", index=False)
    print(f"X01 rebuild: {len(table)} ESR2 entries, {int(table.esr2_chains.sum())} ESR2 chains")


if __name__ == "__main__":
    main()
