"""Tracker S05: row order of results/rmsd_entry.npy as a CSV, so every analysis can load labels without guessing.

Input: results/alignment/rmsd_entry_index.txt (written by S04, same order as the matrix rows).
Output: results/rmsd_entry_index.csv with row, chain_key, entry_id.
Run from the repo root: python scripts/save_rmsd_index.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
INDEX = ROOT / "results/alignment/rmsd_entry_index.txt"
MATRIX = ROOT / "results/rmsd_entry.npy"
OUT = ROOT / "results/rmsd_entry_index.csv"


def main():
    keys = INDEX.read_text().splitlines()
    size = np.load(MATRIX).shape[0]
    if len(keys) != size or len(set(keys)) != size:
        raise ValueError(f"Index has {len(keys)} unique rows; matrix side is {size}")
    table = pd.DataFrame({
        "row": range(size),
        "chain_key": keys,
        "entry_id": [key.split("|")[0] for key in keys],
    })
    table.to_csv(OUT, index=False)
    print(f"S05: {size} rows, {table.entry_id.nunique()} entries written to {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
