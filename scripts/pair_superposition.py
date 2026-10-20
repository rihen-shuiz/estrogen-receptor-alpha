"""Tracker E05: fit each ESR2 chain onto its ESR1 partner on the shared core helices.

Pairs: results/matched_pairs.csv (E04, ESR1/ESR2 entries with the same ligand). Every eligible ESR1 chain of the
ESR1 entry is paired with every eligible ESR2 chain of the ESR2 entry (keep-all: no chain is dropped);
primary_pair marks the historical ESR1 chain (esr_analysis/resources/frozen/rmsd_entry_index.txt) with the first
ESR2 chain in chain_key order, for a one-row-per-entry-pair reading.
Coordinates: deposited C-alpha from results/lbd_extract/residue_mapping.csv. ESR2 residues there already carry
ESR1 numbering (esr1_position, the E03 map), so both chains are compared position by position in ESR1 numbering.
Fit: esr_analysis/align_structures.fit_to_reference, Kabsch on the shared core v3 C-alpha
(esr_analysis/resources/config/alignment_core_v3.json), min_core_ca_for_alignment from analysis_settings.json.
After the fit: C-alpha RMSD over all shared LBD positions and over H12 (h12_range), no refit.
A pair that cannot be fitted keeps its row with blank values and a reason.
Output: results/pair_superposition.csv.
Run from the repo root: PYTHONPATH=. python scripts/pair_superposition.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis._task_io import coordinate_mapping
from esr_analysis.align_structures import fit_to_reference, rmsd
from esr_analysis.run_io import positions

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
LBD = range(SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"] + 1)
H12 = range(SETTINGS["h12_range"][0], SETTINGS["h12_range"][1] + 1)


def shared_rmsd(aligned, reference, window):
    common = sorted(set(aligned) & set(reference) & set(window))
    if not common:
        return np.nan, 0
    return rmsd(np.array([aligned[p] for p in common]), np.array([reference[p] for p in common])), len(common)


def main():
    pairs = pd.read_csv(ROOT / "results/matched_pairs.csv")
    chains = pd.read_csv(ROOT / "results/lbd_chains.csv", low_memory=False,
                         usecols=["chain_key", "pdb_id", "identity", "eligible_historical"])
    eligible = chains[chains.eligible_historical.astype(str).str.lower().eq("true")]
    by_entry = {key: sorted(group.chain_key) for key, group in eligible.groupby(["identity", "pdb_id"])}
    historical = set((ROOT / "esr_analysis/resources/frozen/rmsd_entry_index.txt").read_text().split())
    needed = {k for keys in by_entry.values() for k in keys}
    mapping = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                          usecols=["chain_key", "has_ca", "esr1_position", "x_A", "y_A", "z_A"])
    mapping = mapping[mapping.chain_key.isin(needed)]
    coords = {}
    core = positions(CORE["core_ranges"])

    rows = []
    for pair in pairs.itertuples():
        esr1_keys = by_entry.get(("esr1", pair.esr1_entry), [])
        esr2_keys = by_entry.get(("esr2", pair.esr2_entry), [])
        if not esr1_keys or not esr2_keys:
            rows.append({"esr1_entry": pair.esr1_entry, "esr2_entry": pair.esr2_entry, "ligand": pair.esr1_ligand,
                         "status": "no_chain", "reason": f"eligible chains esr1 {len(esr1_keys)}, esr2 {len(esr2_keys)}"})
            continue
        for k1 in esr1_keys:
            for k2 in esr2_keys:
                for key in (k1, k2):
                    if key not in coords:
                        coords[key] = coordinate_mapping(mapping, key)
                aligned, info = fit_to_reference(coords[k2], coords[k1], core, SETTINGS["min_core_ca_for_alignment"])
                lbd_rmsd, n_lbd = shared_rmsd(aligned, coords[k1], LBD)
                h12_rmsd, n_h12 = shared_rmsd(aligned, coords[k1], H12)
                rows.append({
                    "esr1_entry": pair.esr1_entry, "esr2_entry": pair.esr2_entry, "ligand": pair.esr1_ligand,
                    "esr1_chain": k1, "esr2_chain": k2,
                    "primary_pair": k1 in historical and k2 == esr2_keys[0],
                    "n_common_core_ca": info["n_common_core_ca"], "core_rmsd_A": info["core_alignment_rmsd_A"],
                    "n_common_lbd_ca": n_lbd if aligned else 0, "lbd_rmsd_A": lbd_rmsd if aligned else np.nan,
                    "n_common_h12_ca": n_h12 if aligned else 0, "h12_rmsd_A": h12_rmsd if aligned else np.nan,
                    "core_used": CORE["name"], "status": info["status"], "reason": info["reason"] or "",
                })
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/pair_superposition.csv", index=False)
    ok = table[table.status == "ok"]
    primary = ok[ok.primary_pair == True]  # noqa: E712
    print(f"E05: {len(pairs)} entry pairs -> {len(table)} chain pairs ({len(ok)} fitted, {len(primary)} primary); "
          f"core RMSD median {ok.core_rmsd_A.median():.3f} A (max {ok.core_rmsd_A.max():.3f}); "
          f"LBD median {ok.lbd_rmsd_A.median():.3f} A; H12 median {ok.h12_rmsd_A.median():.3f} A; "
          f"not fitted {int((table.status != 'ok').sum())}")


if __name__ == "__main__":
    main()
