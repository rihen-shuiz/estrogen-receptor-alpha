"""Rebuild tracker X14: match ChEMBL ESR1 compounds to PDB ligands (results/chembl_pdb_matches.csv).

Inputs: results/rebuild/chembl_esr1_activities.csv (X13 rebuild; falls back to data/chembl_esr1_activities.csv
if the rebuild download is absent), results/rebuild/fp_index.csv (C02 rebuild, PDB ligand SMILES).
Method: distinct ChEMBL molecules with a parseable canonical_smiles; ECFP4 (Morgan r=2, 2048 bits, same
generator as C02) on both sides; RDKit Tanimoto; keep every ligand/molecule pair with Tanimoto >= 0.9.
match_type 'exact' when the RDKit canonical isomeric SMILES are equal, else 'tanimoto>=0.9'.
Tanimoto rounded to 4 decimals. Sorted by ligand_id, then Tanimoto descending; ties keep the order in which
molecules first appear in the activity file (sorted by activity_id).
Output: results/rebuild/chembl_pdb_matches.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_chembl_pdb_matches.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from rebuild_common import OUT, ROOT, read_csv


def main():
    RDLogger.DisableLog("rdApp.*")
    source = OUT / "chembl_esr1_activities.csv"
    if not source.exists():
        source = ROOT / "data/chembl_esr1_activities.csv"
    chembl = read_csv(source)
    chembl = chembl[chembl.canonical_smiles.notna()].drop_duplicates("molecule_chembl_id")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    ids, fps, smiles = [], [], []
    for cid, smi in zip(chembl.molecule_chembl_id, chembl.canonical_smiles):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        ids.append(cid)
        fps.append(generator.GetFingerprint(mol))
        smiles.append(Chem.MolToSmiles(mol))
    index = read_csv(OUT / "fp_index.csv")
    rows = []
    for lid, smi in zip(index.ligand_id.astype(str), index.canonical_isomeric_smiles):
        mol = Chem.MolFromSmiles(smi)
        canonical = Chem.MolToSmiles(mol)
        sims = np.array(DataStructs.BulkTanimotoSimilarity(generator.GetFingerprint(mol), fps))
        for j in np.where(sims >= 0.9)[0]:
            exact = smiles[j] == canonical
            rows.append({"ligand_id": lid, "molecule_chembl_id": ids[j], "tanimoto": round(float(sims[j]), 4),
                         "match_type": "exact" if exact else "tanimoto>=0.9", "_exact": exact})
    table = pd.DataFrame(rows).sort_values(["ligand_id", "tanimoto"], ascending=[True, False],
                                           kind="stable").drop(columns="_exact")
    table.to_csv(OUT / "chembl_pdb_matches.csv", index=False)
    print(f"X14 rebuild: {len(table)} matches, {table.ligand_id.nunique()} PDB ligands, "
          f"{table.molecule_chembl_id.nunique()} ChEMBL molecules (from {source.relative_to(ROOT).as_posix()})")


if __name__ == "__main__":
    main()
