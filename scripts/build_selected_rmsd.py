"""Tracker S04 (keep-all): pairwise C-alpha RMSD over ALL eligible ESR1 chains, from the S03 primary core fit.

Rows: every chain in results/alignment/aligned_coordinates.npz (S03, the 788 eligible chains, sorted chain_key order).
Distances: esr_analysis/calculate_rmsd.compute_distances over the common LBD C-alpha positions, H12 included, after the
new core fit. min_pairwise_common_ca comes from analysis_settings.json (150) and is never lowered. A chain is left
out only if it cannot meet that minimum; each one is listed with its reason in results/alignment/distance_exclusions.json.
Outputs: results/rmsd_entry.npy (no NaN), results/alignment/rmsd_entry_index.txt (row order, read by S05),
results/alignment/overlap_entry.npy, results/alignment/distance_exclusions.json.
The archived 351-chain matrix (results/archive_v2_strict) is only read, for a comparison on the shared chains.
Run from the repo root: PYTHONPATH=. /home/rihen/.venvs/esr1/bin/python scripts/build_selected_rmsd.py
"""
import json
from pathlib import Path

import numpy as np

from esr_analysis._task_io import aligned_caches
from esr_analysis.calculate_rmsd import compute_distances

ROOT = Path(__file__).resolve().parents[1]
ALIGN_OUT = ROOT / "results/alignment"
ARCHIVE = ROOT / "results/archive_v2_strict"
DIFFER_A = 1e-6  # a pair counts as different when its RMSD moves by more than this, in angstrom


def main():
    settings = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
    fits = ALIGN_OUT / "aligned_coordinates.npz"
    ids = np.load(fits, allow_pickle=False)["chain_keys"].astype(str).tolist()
    caches = aligned_caches(fits, ids)
    kept, matrix, overlap, excluded = compute_distances(ids, caches, settings["min_pairwise_common_ca"])
    if np.isnan(matrix).any():
        raise ValueError("Matrix still has NaN pairs after pruning; do not save")
    np.save(ROOT / "results/rmsd_entry.npy", matrix)
    np.save(ALIGN_OUT / "overlap_entry.npy", overlap)
    (ALIGN_OUT / "rmsd_entry_index.txt").write_text("\n".join(kept) + "\n")
    (ALIGN_OUT / "distance_exclusions.json").write_text(json.dumps(excluded, indent=2) + "\n")
    print(f"S04: {len(kept)} chains in the matrix; {len(excluded)} left out for overlap; "
          f"{len(kept) * (len(kept) - 1) // 2} pairs; min pairwise overlap {int(overlap[np.triu_indices(len(kept), 1)].min())}")

    archived = ARCHIVE / "rmsd_entry.npy"
    archived_index = (ARCHIVE / "alignment/rmsd_entry_index.txt").read_text().splitlines()
    old = np.load(archived)
    position = {key: i for i, key in enumerate(kept)}
    shared = [key for key in archived_index if key in position]
    new_rows = [position[key] for key in shared]
    old_rows = [archived_index.index(key) for key in shared]
    upper = np.triu_indices(len(shared), k=1)
    diff = np.abs(matrix[np.ix_(new_rows, new_rows)][upper] - old[np.ix_(old_rows, old_rows)][upper])
    print(f"S04 vs archived 351-chain matrix on {len(shared)} shared chains (new core vs old core): "
          f"{int((diff > DIFFER_A).sum())} of {diff.size} pairs differ by more than {DIFFER_A} A; max |diff| {diff.max():.3f} A")


if __name__ == "__main__":
    main()
