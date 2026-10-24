"""Rebuild tracker C02 (ECFP4), C03 (MACCS) and C04 (Tanimoto matrix).

Inputs: results/rebuild/ligands.csv (Q08 rebuild), read text-safe (ligand 'NA' is sodium).
Method: row order = ligands.csv order (sorted CCD code); molecule = RDKit parse of canonical_isomeric_smiles.
C02 Morgan radius 2, 2048 bits, no chirality (rdFingerprintGenerator), uint8 0/1.
C03 RDKit MACCS keys, bits 1..166 (bit 0 is always unset in RDKit and is dropped), uint8.
C04 Tanimoto on the ECFP4 bits, |A&B| / |A|B| (float64; no empty fingerprints occur).
Outputs: results/rebuild/fp_ecfp4.npy, fp_maccs.npy, tanimoto_ecfp4.npy, fp_index.csv (row, ligand_id,
canonical_isomeric_smiles)
Run from the repo root: PYTHONPATH=. python scripts/rebuild_fingerprints.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from rebuild_common import OUT, read_csv


def main():
    ligands = read_csv(OUT / "ligands.csv")
    ligands = ligands[ligands.chemistry_status == "valid"].reset_index(drop=True)
    mols = [Chem.MolFromSmiles(s) for s in ligands.canonical_isomeric_smiles]
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=False)
    ecfp = np.array([generator.GetFingerprintAsNumPy(m) for m in mols], dtype=np.uint8)
    maccs = []
    for m in mols:
        bits = np.zeros(167, dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(MACCSkeys.GenMACCSKeys(m), bits)
        maccs.append(bits[1:])
    maccs = np.array(maccs, dtype=np.uint8)
    b = ecfp.astype(np.int64)
    inter = b @ b.T
    counts = b.sum(1)
    union = counts[:, None] + counts[None, :] - inter
    tanimoto = np.where(union > 0, inter / np.maximum(union, 1), 0.0)
    np.save(OUT / "fp_ecfp4.npy", ecfp)
    np.save(OUT / "fp_maccs.npy", maccs)
    np.save(OUT / "tanimoto_ecfp4.npy", tanimoto)
    pd.DataFrame({"row": range(len(ligands)), "ligand_id": ligands.ligand_id,
                  "canonical_isomeric_smiles": ligands.canonical_isomeric_smiles}).to_csv(OUT / "fp_index.csv", index=False)
    print(f"C02-C04 rebuild: {ecfp.shape} ECFP4, {maccs.shape} MACCS, {tanimoto.shape} Tanimoto")


if __name__ == "__main__":
    main()
