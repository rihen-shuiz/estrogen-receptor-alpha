"""Tracker K06: master chain status table, one row per ESR1 chain-model, none dropped.

Rows: every ESR1 chain in results/chain_loss_reasons.csv (K01, 901 chains). Columns joined, all read, none recomputed:
  completeness_status / completeness_reason  K01 keep-all screen (complete, or the stated reason it has nothing to measure)
  qc_label, qc tiers, warning                K05 results/qc_tiers.csv (per entry)
  h12_b_group, h12_to_core_ratio, b_group    K03 results/helix_bfactor_by_chain.csv (eligible chains only)
  h12_label, h12_rmsd                        S16 results/h12_metric.csv
  core_rmsd                                  S03 results/alignment_core_rmsd.csv
  in_matrix, in_matrix_reason                S04 results/alignment/rmsd_entry_index.txt; false only with a stated reason
  cut_1.0 .. cut_2.5                         S08 results/entry_clusters.csv (blank when not in the matrix)
  old_stage_lost, old_reason, in_old_matrix  K01: why the chain was lost in the old one-chain-per-entry pipeline
Output: results/chain_status.csv.
Run from the repo root: python scripts/chain_status.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    loss = pd.read_csv(ROOT / "results/chain_loss_reasons.csv")
    table = loss[["chain_key", "entry_id", "chain_id", "auth_asym_id", "mapped_lbd_ca", "mapped_core_ca"]].copy()
    table["completeness_status"] = loss.stage_lost_keep_all.isna().map({True: "complete", False: "incomplete"})
    table["completeness_reason"] = loss.reason_keep_all.fillna("")

    qc = pd.read_csv(ROOT / "results/qc_tiers.csv")
    table = table.merge(qc[["entry_id", "resolution_A", "r_free", "mean_b_A2", "tier_resolution", "tier_rfree",
                            "tier_b", "qc_label", "warning"]], on="entry_id", how="left", validate="many_to_one")
    table["qc_label"] = table.qc_label.fillna("unknown")
    table["warning"] = table.warning.fillna("")

    bfac = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv",
                       usecols=["chain_key", "b_group", "h12_b_group", "h12_to_core_ratio"])
    h12 = pd.read_csv(ROOT / "results/h12_metric.csv", usecols=["chain_key", "h12_label", "h12_rmsd"])
    core = pd.read_csv(ROOT / "results/alignment_core_rmsd.csv", usecols=["chain_key", "core_rmsd", "core_rmsd_label"])
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv").drop(columns=["entry_id", "qc_label"])
    for extra in (bfac, h12, core, clusters):
        table = table.merge(extra, on="chain_key", how="left", validate="one_to_one")

    in_matrix = set((ROOT / "results/alignment/rmsd_entry_index.txt").read_text().split())
    excluded = json.loads((ROOT / "results/alignment/distance_exclusions.json").read_text())
    excluded = excluded if isinstance(excluded, dict) else {}
    table["in_matrix"] = table.chain_key.isin(in_matrix)
    table["in_matrix_reason"] = ""
    out = ~table.in_matrix
    table.loc[out, "in_matrix_reason"] = [
        excluded.get(key) or (f"completeness screen: {reason}" if status == "incomplete" else "")
        for key, status, reason in zip(table.chain_key[out], table.completeness_status[out], table.completeness_reason[out])]
    if (table.in_matrix_reason[out] == "").any():
        raise ValueError("a chain is outside the matrix without a stated reason")

    table["old_stage_lost"] = loss.stage_lost.fillna("kept")
    table["old_reason"] = loss.reason.fillna("")
    table["in_old_matrix"] = loss.in_old_matrix
    if len(table) != len(loss) or not table.chain_key.is_unique:
        raise ValueError("row count changed while joining")
    table.to_csv(ROOT / "results/chain_status.csv", index=False)
    print(f"K06: {len(table)} chains, {table.entry_id.nunique()} entries; in_matrix {int(table.in_matrix.sum())}, "
          f"out {int(out.sum())}; completeness {table.completeness_status.value_counts().to_dict()}; "
          f"qc_label {table.qc_label.value_counts().to_dict()}; old pipeline {table.old_stage_lost.value_counts().to_dict()}")


if __name__ == "__main__":
    main()
