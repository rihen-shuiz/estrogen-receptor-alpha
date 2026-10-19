"""Tracker S03 (keep-all): align EVERY eligible ESR1 chain on the new core and on the old core.

Chains: results/lbd_chains.csv rows with identity esr1 and eligible_historical (788). The completeness screen is unchanged.
Reference: the same reference as the last run (esr_analysis/resources/frozen/reference.json, chain 7NFB|B|1).
Primary core: core_ranges in esr_analysis/resources/config/alignment_core_v3.json (written by K04).
Old core: esr1_core_ranges in analysis_settings.json (340-360, 400-420, 450-470, 500-520), fitted and stored too.
Fit: esr_analysis/align_structures.fit_to_reference (Kabsch on shared core C-alpha, min_core from settings).
A chain whose fit cannot run keeps its row with a blank RMSD and a reason. A large core RMSD is labelled
(core_rmsd_label; LARGE_CORE_RMSD_A = 2.0 A, labelling only), never a reason to drop the chain.
Outputs: results/alignment/aligned_coordinates.npz, alignment_transforms.json, global_alignment.csv (primary core),
results/alignment/aligned_coordinates_old_core.npz, alignment_transforms_old_core.json, global_alignment_old_core.csv,
results/alignment/chain_keys.txt and results/alignment_core_rmsd.csv (one row per chain).
The old-core fit is checked against the archived strict run (results/archive_v2_strict/alignment) before saving.
Run from the repo root: PYTHONPATH=. /home/rihen/.venvs/esr1/bin/python scripts/align_selected_chains.py
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from esr_analysis._task_io import coordinate_mapping, json_write
from esr_analysis.align_structures import fit_to_reference
from esr_analysis.run_io import positions

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE_CONFIG = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
CHAINS = ROOT / "results/lbd_chains.csv"
MAPPING = ROOT / "results/lbd_extract/residue_mapping.csv"
REFERENCE = ROOT / "esr_analysis/resources/frozen/reference.json"
ALIGN_OUT = ROOT / "results/alignment"
ARCHIVE = ROOT / "results/archive_v2_strict/alignment"
TABLE = ROOT / "results/alignment_core_rmsd.csv"
LBD = (SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"])
MIN_CORE = SETTINGS["min_core_ca_for_alignment"]
LARGE_CORE_RMSD_A = 2.0  # labelling only, not a QC cut; a chain above it is still kept


def eligible_chain_keys():
    chains = pd.read_csv(CHAINS, keep_default_na=False, na_values=[""], low_memory=False)
    eligible = chains[(chains.identity == "esr1") & chains.eligible_historical.astype(str).str.lower().eq("true")]
    keys = sorted(eligible.chain_key)
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate chain keys among eligible ESR1 chains")
    return keys


def fit_all(mobiles, reference, reference_key, keys, core_ranges):
    canonical = np.arange(LBD[0], LBD[1] + 1)
    core = positions(core_ranges)
    coords, transforms, summary = [], {}, []
    for key in keys:
        aligned, info = fit_to_reference(mobiles[key], reference, core, MIN_CORE)
        transforms[key] = info
        coords.append([aligned.get(int(p), np.full(3, np.nan)) for p in canonical])
        summary.append({"chain_key": key, "pdb_id": key.split("|")[0], "reference_chain_key": reference_key,
                        "n_common_core_ca": info["n_common_core_ca"], "core_alignment_rmsd_A": info["core_alignment_rmsd_A"],
                        "status": info["status"], "reason": info["reason"] or ""})
    return np.asarray(coords), transforms, pd.DataFrame(summary)


def save_fit(tag, keys, coords, transforms, summary):
    canonical = np.arange(LBD[0], LBD[1] + 1)
    suffix = "" if tag == "primary" else "_old_core"
    np.savez_compressed(ALIGN_OUT / f"aligned_coordinates{suffix}.npz", chain_keys=np.asarray(keys),
                        positions=canonical, coords=coords)
    json_write(ALIGN_OUT / f"alignment_transforms{suffix}.json", transforms)
    summary.to_csv(ALIGN_OUT / f"global_alignment{suffix}.csv", index=False)


def main():
    keys = eligible_chain_keys()
    ALIGN_OUT.mkdir(parents=True, exist_ok=True)
    (ALIGN_OUT / "chain_keys.txt").write_text("\n".join(keys) + "\n")
    mapping = pd.read_csv(MAPPING, dtype={"chain_key": str}, keep_default_na=False, na_values=[""], low_memory=False)
    reference_key = json.loads(REFERENCE.read_text())["chain_key"]
    reference = coordinate_mapping(mapping, reference_key)
    mobiles = {key: coordinate_mapping(mapping, key) for key in keys}

    new_coords, new_transforms, new_summary = fit_all(mobiles, reference, reference_key, keys, CORE_CONFIG["core_ranges"])
    old_coords, old_transforms, old_summary = fit_all(mobiles, reference, reference_key, keys, SETTINGS["esr1_core_ranges"])

    archived = pd.read_csv(ARCHIVE / "global_alignment.csv").set_index("chain_key")
    diff = np.abs(old_summary.set_index("chain_key").core_alignment_rmsd_A - archived.loc[keys, "core_alignment_rmsd_A"])
    archive_npz = np.load(ARCHIVE / "aligned_coordinates.npz", allow_pickle=False)
    archived_keys = archive_npz["chain_keys"].astype(str).tolist()
    archived_coords = archive_npz["coords"][[archived_keys.index(k) for k in keys]]
    same_gaps = np.array_equal(np.isnan(old_coords), np.isnan(archived_coords))
    coord_diff = float(np.nanmax(np.abs(old_coords - archived_coords)))
    if float(diff.max()) > 1e-12 or coord_diff != 0.0 or not same_gaps:
        raise ValueError(f"Old-core fit differs from the archived strict run: RMSD {diff.max()} A, coords {coord_diff} A, gaps equal {same_gaps}")
    print(f"S03: old-core check vs archived strict run: RMSD max diff {float(diff.max()):.1e} A, coordinates max diff {coord_diff} A")

    save_fit("primary", keys, new_coords, new_transforms, new_summary)
    save_fit("old", keys, old_coords, old_transforms, old_summary)

    label_threshold = LARGE_CORE_RMSD_A
    table = pd.DataFrame({
        "entry_id": [key.split("|")[0] for key in keys],
        "chain_key": keys,
        "core_rmsd": new_summary.core_alignment_rmsd_A.to_numpy(),
        "core_rmsd_old_core": old_summary.core_alignment_rmsd_A.to_numpy(),
        "core_used": CORE_CONFIG["name"],
        "core_rmsd_label": np.where(new_summary.status.eq("ok"),
                                    np.where(new_summary.core_alignment_rmsd_A > label_threshold, "large", "ok"), "no_fit"),
        "n_common_core_ca": new_summary.n_common_core_ca.to_numpy(),
        "n_common_core_ca_old": old_summary.n_common_core_ca.to_numpy(),
        "status": new_summary.status.to_numpy(),
        "reason": new_summary.reason.to_numpy(),
    })
    table.to_csv(TABLE, index=False)
    new_ok = table.core_rmsd.dropna()
    print(f"S03: {len(table)} eligible chains; reference {reference_key}; primary core residues {sum(b - a + 1 for a, b in CORE_CONFIG['core_ranges'])}")
    print(f"S03: primary core RMSD min {new_ok.min():.3f} median {new_ok.median():.3f} max {new_ok.max():.3f} A; "
          f"labelled large (> {label_threshold} A): {int((table.core_rmsd_label == 'large').sum())}; no fit: {int((table.core_rmsd_label == 'no_fit').sum())}")


if __name__ == "__main__":
    main()
