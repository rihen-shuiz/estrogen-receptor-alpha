"""SV14 diagnostic: where the residues that fail the apo check sit (descriptive; the apo rule itself is in
SV14_extract.py and is not changed by this).

For each blocking residue in data/sv/SV14/apo_check.csv (residue other than the chain's own polymer or water, >= 6
heavy atoms, within 4 A of the ligand-pocket residues): heavy atoms, distance of its centroid and of its nearest heavy
atom to the E2 reference centroid (data/sv/SV14/e2_centroid.json), core v3 frame.
Output: results/sv/SV14/apo_blocking_residues.csv
Run: .venv/bin/python -B scripts/sv/SV14_apo_blockers.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts/sv"))
from esr_analysis.read_structure import read_cif  # noqa: E402
from common.core_frame import load_transforms, to_core  # noqa: E402


def main():
    apo = pd.read_csv(ROOT / "data/sv/SV14/apo_check.csv", keep_default_na=False)
    e2c = np.array(json.loads((ROOT / "data/sv/SV14/e2_centroid.json").read_text())["centroid_A"])
    T = load_transforms()
    rows = []
    for r in apo[apo.residues_ge6_within_4A != ""].itertuples():
        _, atoms = read_cif(ROOT / "data/raw" / f"{r.entry_id}.cif")
        atoms = atoms[(atoms.pdbx_PDB_model_num == "1") & ~atoms.type_symbol.str.upper().isin(["H", "D"])]
        ins = atoms.pdbx_PDB_ins_code.where(~atoms.pdbx_PDB_ins_code.isin(["?", "."]), "")
        for res in r.residues_ge6_within_4A.split(";"):
            asym, comp, rk, seq = res.split("|")
            a = atoms[(atoms.label_asym_id == asym) & (atoms.label_comp_id == comp)
                      & ((atoms.auth_seq_id.astype(str) + ins) == rk) & (atoms.label_seq_id == seq)]
            x = to_core(a[["Cartn_x", "Cartn_y", "Cartn_z"]].to_numpy(float), T[r.chain_key])
            rows.append({"chain_key": r.chain_key, "entry_id": r.entry_id, "residue": res, "heavy_atoms": len(a),
                         "centroid_to_e2_centroid_A": float(np.linalg.norm(x.mean(0) - e2c)),
                         "nearest_atom_to_e2_centroid_A": float(np.linalg.norm(x - e2c, axis=1).min())})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/sv/SV14/apo_blocking_residues.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
