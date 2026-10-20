"""Tracker E06 / E07: per-residue C-alpha displacement between matched ESR1 and ESR2 chains after the E05 core fit.

Pairs: results/pair_superposition.csv (E05) rows with primary_pair True and status ok, in file order
(one chain pair per matched entry pair). Each ESR2 chain is fitted onto its ESR1 partner exactly as in E05
(scripts/pair_superposition.py: Kabsch on shared core v3 C-alpha), then the distance between equivalent C-alpha
atoms is measured at every shared LBD position in ESR1 numbering. The refit core RMSD must equal the E05 value.
ESR2 numbering and residue identities come from results/numbering_map_esr1_esr2.csv (E03).
  --first 5 (E06): long table, results/first_five_displacements.csv (one row per pair and position).
  --all     (E07): wide table, results/displacement_matrix.csv (rows = ESR1 positions, one column per pair;
            blank where a position is unobserved in either chain) plus results/displacement_pairs.csv (column key).
Run from the repo root: PYTHONPATH=. python scripts/receptor_displacements.py --first 5 | --all
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pair_superposition import CORE, H12, LBD, SETTINGS  # noqa: E402

from esr_analysis._task_io import coordinate_mapping  # noqa: E402
from esr_analysis.align_structures import fit_to_reference  # noqa: E402
from esr_analysis.run_io import positions  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CORE_POSITIONS = set(positions(CORE["core_ranges"]))
SAME_FIT_TOL_A = 1e-9


def displacements(pairs):
    keys = set(pairs.esr1_chain) | set(pairs.esr2_chain)
    mapping = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                          usecols=["chain_key", "has_ca", "esr1_position", "x_A", "y_A", "z_A"])
    mapping = mapping[mapping.chain_key.isin(keys)]
    coords = {key: coordinate_mapping(mapping, key) for key in keys}
    for i, pair in enumerate(pairs.itertuples()):
        aligned, info = fit_to_reference(coords[pair.esr2_chain], coords[pair.esr1_chain], sorted(CORE_POSITIONS),
                                         SETTINGS["min_core_ca_for_alignment"])
        if abs(info["core_alignment_rmsd_A"] - pair.core_rmsd_A) > SAME_FIT_TOL_A:
            raise ValueError(f"{pair.esr1_chain} / {pair.esr2_chain}: refit core RMSD differs from E05")
        reference = coords[pair.esr1_chain]
        for p in sorted(set(aligned) & set(reference) & set(LBD)):
            yield i, pair, p, float(np.linalg.norm(aligned[p] - reference[p]))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--first", type=int)
    mode.add_argument("--all", action="store_true")
    args = parser.parse_args()

    table = pd.read_csv(ROOT / "results/pair_superposition.csv")
    pairs = table[(table.primary_pair == True) & (table.status == "ok")].reset_index(drop=True)  # noqa: E712
    if args.first:
        pairs = pairs.head(args.first)
    numbering = pd.read_csv(ROOT / "results/numbering_map_esr1_esr2.csv").set_index("esr1_pos")

    rows = []
    for i, pair, p, d in displacements(pairs):
        mapped = numbering.loc[p] if p in numbering.index else None
        rows.append({"pair": i + 1, "esr1_entry": pair.esr1_entry, "esr2_entry": pair.esr2_entry,
                     "ligand": pair.ligand, "esr1_chain": pair.esr1_chain, "esr2_chain": pair.esr2_chain,
                     "esr1_pos": p, "esr2_pos": mapped.esr2_pos if mapped is not None else np.nan,
                     "esr1_aa": mapped.esr1_aa if mapped is not None else "",
                     "esr2_aa": mapped.esr2_aa if mapped is not None else "",
                     "same_residue": mapped.same_residue if mapped is not None else "",
                     "in_core": p in CORE_POSITIONS, "in_h12": p in H12, "displacement_A": d})
    long = pd.DataFrame(rows)

    if args.first:
        long.to_csv(ROOT / "results/first_five_displacements.csv", index=False)
        per_pair = long.groupby("pair").displacement_A.agg(["size", "median", "max"])
        print(f"E06: {len(pairs)} pairs, {len(long)} rows\n{per_pair.to_string()}")
        return

    long["column"] = long.esr1_chain + "~" + long.esr2_chain
    wide = long.pivot(index="esr1_pos", columns="column", values="displacement_A")
    wide = wide.reindex(index=list(LBD), columns=(pairs.esr1_chain + "~" + pairs.esr2_chain).tolist())
    wide.index.name = "esr1_pos"
    wide.to_csv(ROOT / "results/displacement_matrix.csv")
    key = pairs[["esr1_entry", "esr2_entry", "ligand", "esr1_chain", "esr2_chain", "core_rmsd_A"]].copy()
    key.insert(0, "column", pairs.esr1_chain + "~" + pairs.esr2_chain)
    key.to_csv(ROOT / "results/displacement_pairs.csv", index=False)
    median = wide.median(axis=1)
    print(f"E07: {wide.shape[0]} positions x {wide.shape[1]} pairs; observed cells {int(wide.notna().sum().sum())}; "
          f"largest median displacement at {median.idxmax()} ({median.max():.2f} A)")


if __name__ == "__main__":
    main()
