"""Core v3 frame shared by SV03, SV04, SV07 and SV13 (RUN_SV1028, section 2: the superposition core stays v3).

Source: the SV00 reproduction-gate fit (S03 code path, core esr_analysis/resources/config/alignment_core_v3.json
core_ranges, reference chain 7NFB|B|1), saved by SV00 in data/sv/gate/:
  - alignment_transforms.json: per chain rotation R (3x3), mobile_centroid m and reference_centroid c;
  - aligned_coordinates.npz: C-alpha of every kept chain in the core v3 frame (chain_keys, positions 305-550, coords).
Convention (checked in SV03 against aligned_coordinates.npz, max |d| = 0 A):
  x_core = (x_deposited - m) @ R + c          (row vectors)
SV03 freezes a copy of the transforms it used in data/sv/SV03/core_v3_transforms.csv (one row per chain) so that SV04
and SV13 apply exactly the same frame to any deposited coordinates (backbone, side chains, ligands).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
GATE = ROOT / "data/sv/gate"
FROZEN = ROOT / "data/sv/SV03/core_v3_transforms.csv"
R_COLS = [f"r{i}{j}" for i in range(3) for j in range(3)]
M_COLS = ["mobile_cx", "mobile_cy", "mobile_cz"]
C_COLS = ["reference_cx", "reference_cy", "reference_cz"]


def transforms_from_gate() -> pd.DataFrame:
    """The SV00 gate transforms as a table (one row per chain)."""
    j = json.loads((GATE / "alignment_transforms.json").read_text())
    rows = []
    for key, t in j.items():
        R = np.asarray(t["rotation"], float)
        rows.append({"chain_key": key, "status": t["status"], "n_common_core_ca": t["n_common_core_ca"],
                     "core_alignment_rmsd_A": t["core_alignment_rmsd_A"],
                     **dict(zip(R_COLS, R.ravel())), **dict(zip(M_COLS, t["mobile_centroid"])),
                     **dict(zip(C_COLS, t["reference_centroid"]))})
    return pd.DataFrame(rows)


def load_transforms(path: Path | None = None) -> dict:
    """{chain_key: (R, m, c)} from the frozen SV03 copy (default) or another table with the same columns."""
    t = pd.read_csv(path or FROZEN)
    return {r.chain_key: (np.array([getattr(r, c) for c in R_COLS]).reshape(3, 3),
                          np.array([getattr(r, c) for c in M_COLS]), np.array([getattr(r, c) for c in C_COLS]))
            for r in t.itertuples()}


def to_core(xyz, transform) -> np.ndarray:
    """Deposited coordinates (..., 3) -> core v3 frame. NaN rows stay NaN."""
    R, m, c = transform
    return (np.asarray(xyz, float) - m) @ R + c


def load_ca_core():
    """C-alpha in the core v3 frame: (chain_keys [n], positions [246], coords [n, 246, 3], NaN = not observed)."""
    d = np.load(GATE / "aligned_coordinates.npz")
    return d["chain_keys"].astype(str), d["positions"].astype(int), d["coords"]
