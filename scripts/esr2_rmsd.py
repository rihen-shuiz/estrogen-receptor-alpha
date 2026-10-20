"""Tracker X02: ESR2 pairwise C-alpha RMSD matrix, aligned to an ESR2 reference on the core helices.

Chains: every eligible ESR2 chain in results/lbd_chains.csv (identity esr2, eligible_historical; keep-all, none dropped
for quality). ESR2 residues carry ESR1 numbering (esr1_position in results/lbd_extract/residue_mapping.csv, the E03
map), so the ESR1 core v3 helices and LBD range are used unchanged.
Reference: the ESR1 block-08 policy applied to ESR2 (esr_analysis/choose_reference.py): quality pass, LBD fraction
>= reference_min_lbd_fraction, H12 fraction >= its minimum, then best resolution, most mapped LBD C-alpha, fewest
mismatches, chain_key. Stated deviation: the core-fraction minimum (1.0, on the old ESR1 core) is not applied, because
no ESR2 chain can reach it (the best of 86 is 0.976; some old-core positions have no observed ESR2 equivalent in any
chain). The deviation and the best reachable core fraction are written to esr2_reference.json. If no chain meets the
remaining policy the script stops (no silent relaxing).
Fit: esr_analysis/align_structures.fit_to_reference on shared core v3 C-alpha (min_core_ca_for_alignment).
Distances: esr_analysis/calculate_rmsd.compute_distances over shared LBD C-alpha (min_pairwise_common_ca), as S04.
Outputs: results/rmsd_esr2.npy, results/rmsd_esr2_index.txt (row order), results/esr2_alignment.csv (per-chain fit),
results/esr2_reference.json (chosen reference and the candidate count).
Run from the repo root: PYTHONPATH=. python scripts/esr2_rmsd.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis._task_io import coordinate_mapping
from esr_analysis.align_structures import fit_to_reference
from esr_analysis.calculate_rmsd import compute_distances
from esr_analysis.run_io import positions

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
LBD = range(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)


def main():
    chains = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False)
    esr2 = chains[(chains.identity == "esr2") & chains.eligible_historical.astype(str).str.lower().eq("true")].copy()
    entries = pd.read_csv(ROOT / "results/esr2_entries.csv").rename(columns={"entry_id": "pdb_id"})
    esr2 = esr2.merge(entries[["pdb_id", "resolution", "r_free", "qc_pass"]], on="pdb_id", how="left", validate="many_to_one")
    keys = sorted(esr2.chain_key)

    candidates = esr2[(esr2.quality == "pass") & (esr2.lbd_fraction >= SETTINGS["reference_min_lbd_fraction"])
                      & (esr2.h12_fraction >= SETTINGS["reference_min_h12_fraction"])]
    if candidates.empty:
        raise RuntimeError("No ESR2 chain meets the reference policy; inspect lbd_chains.csv, do not relax silently")
    chosen = candidates.sort_values(["resolution", "mapped_lbd_ca", "mismatch_count", "chain_key"],
                                    ascending=[True, False, True, True]).iloc[0]

    mapping = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                          usecols=["chain_key", "has_ca", "esr1_position", "x_A", "y_A", "z_A"])
    mapping = mapping[mapping.chain_key.isin(keys)]
    coords = {key: coordinate_mapping(mapping, key) for key in keys}
    reference = coords[chosen.chain_key]
    core = positions(CORE["core_ranges"])
    aligned, summary = {}, []
    for key in keys:
        fitted, info = fit_to_reference(coords[key], reference, core, SETTINGS["min_core_ca_for_alignment"])
        aligned[key] = {p: xyz for p, xyz in fitted.items() if p in LBD}
        summary.append({"chain_key": key, "entry_id": key.split("|")[0], "reference_chain_key": chosen.chain_key,
                        "n_common_core_ca": info["n_common_core_ca"], "core_rmsd_A": info["core_alignment_rmsd_A"],
                        "n_lbd_ca": len(aligned[key]), "status": info["status"], "reason": info["reason"] or ""})
    summary = pd.DataFrame(summary)
    fitted_keys = summary.loc[summary.status == "ok", "chain_key"].tolist()
    kept, matrix, overlap, excluded = compute_distances(fitted_keys, aligned, SETTINGS["min_pairwise_common_ca"])
    if np.isnan(matrix).any():
        raise ValueError("ESR2 matrix has NaN pairs after pruning; do not save")
    summary["in_matrix"] = summary.chain_key.isin(kept)
    np.save(ROOT / "results/rmsd_esr2.npy", matrix)
    (ROOT / "results/rmsd_esr2_index.txt").write_text("\n".join(kept) + "\n")
    summary.to_csv(ROOT / "results/esr2_alignment.csv", index=False)
    record = {"chain_key": chosen.chain_key, "pdb_id": chosen.pdb_id, "resolution_A": float(chosen.resolution),
              "lbd_fraction": float(chosen.lbd_fraction), "mapped_lbd_ca": int(chosen.mapped_lbd_ca),
              "n_candidates_meeting_policy": len(candidates), "n_chains": len(keys), "n_in_matrix": len(kept),
              "core_fraction": float(chosen.core_fraction), "best_core_fraction_any_esr2": float(esr2.core_fraction.max()),
              "excluded": excluded, "core": CORE["name"],
              "policy": "ESR1 block-08 policy on ESR2: quality pass, LBD and H12 fraction minimums, best resolution",
              "deviation": "reference_min_core_fraction 1.0 not applied: unreachable for ESR2 (best 0.976)"}
    (ROOT / "results/esr2_reference.json").write_text(json.dumps(record, indent=2, default=str) + "\n")
    upper = np.triu_indices(len(kept), 1)
    print(f"X02: reference {chosen.chain_key} ({chosen.resolution} A, {len(candidates)} candidates); {len(keys)} chains, "
          f"{len(kept)} in matrix, {len(excluded)} excluded; core RMSD median {summary.core_rmsd_A.median():.3f} A; "
          f"pairwise RMSD median {np.median(matrix[upper]):.3f} A (max {matrix[upper].max():.3f}); "
          f"min overlap {int(overlap[upper].min())}")


if __name__ == "__main__":
    main()
