"""Tracker X09: top P2Rank pocket per chain, its size proxy and score, and whether it overlaps the real ligand.

Inputs: results/p2rank/<entry>_<chain>.pdb_predictions.csv (X08), data/clean/<entry>_<chain>.pdb (ligand of interest
as HETATM; scripts/clean_chain_pdbs.py), results/chain_ligand_of_interest.csv (Q18), pocket_positions of Q18 from
scripts/ligand_of_interest.py (POCKET_POSITIONS).
Keep-all: one row per eligible chain (788), entry_id repeated; nothing is filtered.
Top pocket = P2Rank rank 1. Size proxy = sas_points (solvent-accessible points P2Rank assigns to the pocket; P2Rank
reports no volume). overlaps_ligand: DCA criterion used in P2Rank papers, distance from pocket centre to the nearest
ligand heavy atom <= DCA_A; blank for chains with no ligand of interest. first_overlapping_rank = rank of the first
pocket that meets DCA (blank if none). q18_residue_jaccard compares the top pocket residues with the Q18 pocket list.
Output: results/pockets.csv.
Run from the repo root: python scripts/pocket_table.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
DCA_A = 4.0


def q18_pocket():
    """POCKET_POSITIONS from scripts/ligand_of_interest.py, read without importing it (it needs Biopython)."""
    import ast
    tree = ast.parse((ROOT / "scripts/ligand_of_interest.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "POCKET_POSITIONS":
            return set(ast.literal_eval(node.value))
    raise ValueError("POCKET_POSITIONS not found")


def ligand_atoms(pdb):
    """Heavy atoms of the ligand section: HETATM records after the chain's TER (modified residues before TER are protein)."""
    lines = pdb.read_text().splitlines()
    after_ter = lines[next((i for i, l in enumerate(lines) if l.startswith("TER")), len(lines)):]
    rows = [l for l in after_ter if l.startswith("HETATM") and l[76:78].strip() not in ("H", "D")]
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows]) if rows else np.empty((0, 3))


def main():
    pocket_ref = q18_pocket()
    loi = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv", keep_default_na=False)
    ligands = loi[loi.ligand_id != ""].groupby("chain_key").ligand_id.agg(lambda s: ";".join(sorted(set(s))))
    chains = sorted(set(loi.chain_key))
    rows = []
    for key in chains:
        entry, chain, _ = key.split("|")
        stem = f"{entry}_{chain}"
        pred = pd.read_csv(ROOT / f"results/p2rank/{stem}.pdb_predictions.csv", skipinitialspace=True)
        pred.columns = [c.strip() for c in pred.columns]
        lig = ligand_atoms(ROOT / f"data/clean/{stem}.pdb")
        row = {"entry_id": entry, "chain_key": key, "ligand_id": ligands.get(key, ""), "n_pockets": len(pred)}
        dca = None
        if len(lig) and len(pred):
            centres = pred[["center_x", "center_y", "center_z"]].to_numpy(float)
            dca = np.sqrt(((centres[:, None, :] - lig[None, :, :]) ** 2).sum(-1)).min(axis=1)
        if len(pred):
            top = pred.iloc[0]
            residues = {int(re.match(r"-?\d+", r.split("_")[1]).group()) for r in str(top.residue_ids).split()}
            row.update(pocket_score=float(top.score), pocket_probability=float(top.probability),
                       sas_points=int(top.sas_points), surf_atoms=int(top.surf_atoms), n_residues=len(residues),
                       residue_ids=str(top.residue_ids).strip(),
                       q18_residue_jaccard=len(residues & pocket_ref) / len(residues | pocket_ref))
        if dca is not None:
            hits = np.flatnonzero(dca <= DCA_A)
            row.update(dca_top_A=float(dca[0]), overlaps_ligand=bool(dca[0] <= DCA_A),
                       first_overlapping_rank=int(pred["rank"].iloc[hits[0]]) if len(hits) else np.nan)
        else:
            row.update(dca_top_A=np.nan, overlaps_ligand=pd.NA, first_overlapping_rank=np.nan)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(ROOT / "results/pockets.csv", index=False)
    with_lig = table[table.ligand_id != ""]
    print(f"X09: {len(table)} chains ({table.entry_id.nunique()} entries); top pocket overlaps ligand in "
          f"{int(with_lig.overlaps_ligand.astype(bool).sum())}/{len(with_lig)} chains with a ligand; "
          f"median score {table.pocket_score.median():.2f}, median sas_points {table.sas_points.median():.0f}, "
          f"median Jaccard vs Q18 pocket {table.q18_residue_jaccard.median():.2f}")


if __name__ == "__main__":
    main()
