"""Chemistry cluster x structural cluster test for the September entries (tracker task P04).

Each September entry (one chain, from results/task_outputs/block_10 index) gets its ligand of interest
from results/chain_ligand_of_interest.csv (Q18). A chain with several ligands keeps the one with the most
pocket residues in contact; ties go to the smallest ligand_key. Entries with no ligand of interest are
excluded and counted. The chemistry cluster of the ligand comes from results/ligand_clusters_loi.csv (C11,
K=7), the structural cluster from the 1.5 Å cut of P01 (results/provisional/clusters.csv).

Test 1: Pearson chi-squared of the contingency table (no continuity correction), with the count of
cells whose expected count is below 5. Test 2: permutation p with 10000 permutations, shuffling the
chemistry cluster across ligands (all entries of one ligand move together), seed 42. p = (k + 1) / (N + 1).
Writes results/provisional/chem_struct_test.csv and results/provisional/chem_struct_contingency.csv.

Run from anywhere: python scripts/provisional_chem_struct.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provisional_common import PROVISIONAL, ROOT, SEED, load_september, entry_id  # noqa: E402

N_PERMUTATIONS = 10000
STRUCT_COLUMN = "cut_1.5"


def chi2_from_counts(counts):
    total = counts.sum()
    expected = np.outer(counts.sum(axis=1), counts.sum(axis=0)) / total
    kept = expected > 0
    return float(((counts - expected)[kept] ** 2 / expected[kept]).sum())


def read(path):
    return pd.read_csv(ROOT / path, dtype=str, keep_default_na=False)


def ligand_per_entry(chain_keys):
    """Ligand of interest per September chain: (ligand_id, n_multi_resolved) or None if apo."""
    rows = read("results/chain_ligand_of_interest.csv")
    chosen, multi = {}, 0
    for key in chain_keys:
        found = rows[(rows.chain_key == key) & (rows.ligand_id != "")].copy()
        if found.empty:
            chosen[key] = None
            continue
        if len(found) > 1:
            multi += 1
        found["pocket"] = found.n_pocket_residues.astype(int)
        best = found.sort_values(["pocket", "ligand_key"], ascending=[False, True]).iloc[0]
        chosen[key] = best.ligand_id
    return chosen, multi


def main():
    chain_keys, _ = load_september()
    chosen, multi = ligand_per_entry(chain_keys)
    chemistry = read("results/ligand_clusters_loi.csv")
    chem_of = dict(zip(chemistry.ligand_id, chemistry.chem_cluster.astype(int)))
    clusters = pd.read_csv(PROVISIONAL / "clusters.csv", dtype=str, keep_default_na=False)
    struct_of = dict(zip(clusters.chain_key, clusters[STRUCT_COLUMN].astype(int)))

    included = pd.DataFrame([{"entry_id": entry_id(key), "chain_key": key, "ligand_id": chosen[key],
                              "struct": struct_of[key]} for key in chain_keys if chosen[key] is not None])
    excluded = len(chain_keys) - len(included)
    missing = sorted(set(included.ligand_id) - set(chem_of))
    if missing:
        raise ValueError(f"ligands without a chemistry cluster: {missing}")
    included["chem"] = included.ligand_id.map(chem_of)

    chem_levels, struct_levels = sorted(included.chem.unique()), sorted(included.struct.unique())
    counts = pd.crosstab(included.chem, included.struct).reindex(index=chem_levels, columns=struct_levels,
                                                                  fill_value=0)
    observed = counts.to_numpy()
    stat, p_chi2, dof, expected = chi2_contingency(observed, correction=False)
    assert np.isclose(stat, chi2_from_counts(observed))
    cells_below_5 = int((expected < 5).sum())

    ligands = sorted(included.ligand_id.unique())
    ligand_index = {ligand: i for i, ligand in enumerate(ligands)}
    entry_ligand = included.ligand_id.map(ligand_index).to_numpy()
    ligand_chem = np.array([chem_levels.index(chem_of[ligand]) for ligand in ligands])
    struct_index = included.struct.map({level: i for i, level in enumerate(struct_levels)}).to_numpy()
    n_chem, n_struct = len(chem_levels), len(struct_levels)
    rng = np.random.default_rng(SEED)
    hits = 0
    for _ in range(N_PERMUTATIONS):
        entry_chem = rng.permutation(ligand_chem)[entry_ligand]
        table = np.bincount(entry_chem * n_struct + struct_index, minlength=n_chem * n_struct).reshape(n_chem, n_struct)
        if chi2_from_counts(table) >= stat - 1e-9:
            hits += 1
    p_perm = (hits + 1) / (N_PERMUTATIONS + 1)

    base = {"n": len(included), "df": int(dof), "cells_expected_below_5": cells_below_5,
            "n_cells": int(expected.size), "n_ligands": len(ligands), "n_excluded_no_ligand": excluded,
            "n_multi_resolved": multi}
    results = pd.DataFrame([
        {"test": "chi_squared_independence", **base, "statistic": stat, "p_value": p_chi2,
         "method": "pearson_chi2_no_correction", "n_permutations": 0, "seed": SEED},
        {"test": "permutation_ligand_level", **base, "statistic": stat, "p_value": p_perm,
         "method": "chemistry_label_shuffled_across_ligands", "n_permutations": N_PERMUTATIONS, "seed": SEED},
    ])
    PROVISIONAL.mkdir(parents=True, exist_ok=True)
    results.to_csv(PROVISIONAL / "chem_struct_test.csv", index=False)
    counts.index.name = "chem_cluster"
    counts.columns = [f"struct_{column}" for column in counts.columns]
    counts.to_csv(PROVISIONAL / "chem_struct_contingency.csv")
    print(results.to_string(index=False))
    print(counts.to_string())
    print(f"wrote {PROVISIONAL / 'chem_struct_test.csv'} and {PROVISIONAL / 'chem_struct_contingency.csv'}")


if __name__ == "__main__":
    main()
