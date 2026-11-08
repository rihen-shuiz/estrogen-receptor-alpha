"""SV14 diagnostic: where the pocket component of the E2 reference chains extends beyond the ligand.

Descriptive only (no rule depends on it): the wild-type EST reference pockets fall in two groups of volume, so this
lists, for each reference chain, the pocket points farther than 4.0 A from every ligand heavy atom and the receptor
residues with a heavy atom within 4.5 A of those points (UniProt numbering). Same grid code as SV14_grid.py
(chain_grids, select_component; buried >= 10 rays).
Input: data/sv/SV14/atoms_core.npz, e2_centroid.json; results/sv/SV11/site_residues.csv.
Output: results/sv/SV14/e2_reference_pocket_lining.csv
Run: .venv/bin/python -B scripts/sv/SV14_e2_lining.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV14_grid as G  # noqa: E402

FAR_A = 4.0
LINING_A = 4.5


def main():
    a = np.load(G.DATA / "atoms_core.npz")
    keys = list(a["chain_keys"].astype(str))
    rs = np.concatenate([[0], np.cumsum(a["rec_count"])])
    ls = np.concatenate([[0], np.cumsum(a["lig_count"])])
    s = pd.read_csv(ROOT / "results/sv/SV11/site_residues.csv")
    site = np.array(sorted(s[(s.site == "ligand pocket") & (s.in_site.astype(str) == "True")].uniprot_position
                           .astype(int)))
    e2 = json.loads((G.DATA / "e2_centroid.json").read_text())
    e2c = np.array(e2["centroid_A"])
    six = ndimage.generate_binary_structure(3, 1)
    rows = []
    for k in e2["chains"]:
        i = keys.index(k)
        X, R, P = (a["rec_xyz"][rs[i]:rs[i + 1]].astype(float), a["rec_rad"][rs[i]:rs[i + 1]].astype(float),
                   a["rec_pos"][rs[i]:rs[i + 1]])
        Lx, Lr = a["lig_xyz"][ls[i]:ls[i + 1]].astype(float), a["lig_rad"][ls[i]:ls[i + 1]].astype(float)
        gr = G.chain_grids(X, R, P, Lx, Lr, site)
        cand = gr["cav"] & (gr["hits"] >= G.PRIMARY_T)
        lab, ncomp = ndimage.label(cand, structure=six)
        L, _, _ = G.select_component(lab, ncomp, gr["lig"], cand, True, gr["lo"], e2c)
        pts = gr["lo"] + np.argwhere(lab == L) * G.H
        dl, _ = cKDTree(Lx).query(pts)
        far = pts[dl > FAR_A]
        near = cKDTree(X).query_ball_point(far, LINING_A) if len(far) else []
        res = pd.Series([int(P[j]) for lst in near for j in lst], dtype=int)
        counts = res.value_counts()
        top = counts.head(12)
        rows.append({"chain_key": k, "pocket_points": int(len(pts)), "pocket_volume_A3": len(pts) * G.VOXEL,
                     "points_far_from_ligand": int(len(far)), "volume_far_from_ligand_A3": len(far) * G.VOXEL,
                     "lining_residues_far_part_top12": ";".join(f"{r}:{c}" for r, c in top.items()),
                     "n_lining_residues_far_part": int(len(counts))})
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "results/sv/SV14/e2_reference_pocket_lining.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
