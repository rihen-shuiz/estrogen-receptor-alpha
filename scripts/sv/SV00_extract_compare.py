"""SV00 step 4 check: regenerated Q04 extraction (data/sv) versus the frozen tables (results/), file by file.
Per pair: byte-identical (sha256), rows, columns, and for non-identical files the largest numeric difference and the
number of text cells that differ after reading every column as text-safe values.
Output: results/sv/SV00/extract_compare.csv.
Run: .venv/bin/python -I -B scripts/sv/SV00_extract_compare.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from esr_analysis.run_io import sha256  # noqa: E402

PAIRS = [("data/sv/lbd_extract/residue_mapping.csv", "results/lbd_extract/residue_mapping.csv"),
         ("data/sv/lbd_extract/receptor_chains.csv", "results/lbd_extract/receptor_chains.csv"),
         ("data/sv/lbd_extract/entries.csv", "results/lbd_extract/entries.csv"),
         ("data/sv/lbd_chains.csv", "results/lbd_chains.csv")]


def rd(p):
    return pd.read_csv(ROOT / p, keep_default_na=False, na_values=[""], low_memory=False)


rows = []
for new, old in PAIRS:
    same = sha256(ROOT / new) == sha256(ROOT / old)
    a, b = rd(new), rd(old)
    max_num, text_cells, num_cols = 0.0, 0, []
    if not same and a.shape == b.shape and list(a.columns) == list(b.columns):
        for c in a.columns:
            x, y = a[c], b[c]
            if x.dtype.kind in "fi" and y.dtype.kind in "fi":
                d = float(np.nanmax(np.abs(x.to_numpy(float) - y.to_numpy(float)))) if len(x) else 0.0
                if d > 0:
                    num_cols.append(c)
                max_num = max(max_num, d)
            else:
                both_na = x.isna() & y.isna()
                text_cells += int(((x.astype(str) != y.astype(str)) & ~both_na).sum())
    rows.append({"regenerated": new, "frozen": old, "byte_identical": same, "rows_regenerated": len(a),
                 "rows_frozen": len(b), "same_columns": list(a.columns) == list(b.columns),
                 "max_abs_numeric_diff": max_num, "numeric_columns_differing": ";".join(num_cols),
                 "text_cells_differing": text_cells})
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/sv/SV00/extract_compare.csv", index=False)
print(out.to_string(index=False))
