"""Shared helpers for the scripts/rebuild_*.py family (not a tracker task by itself).

The rebuild scripts recreate frozen upstream result files from raw inputs and write
only under results/rebuild/. This module holds paths, settings, input lookup (raw
mmCIF and CCD files, local copy first, then the rebuild download cache) and the
comparison helpers that fill results/rebuild/rebuild_report.csv.

Not run directly; imported by the rebuild scripts.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
SEED = 42
OUT = ROOT / "results/rebuild"
RAW_LOCAL = ROOT / "data/raw"
RAW_CACHE = OUT / "raw"            # mmCIFs missing from data/raw (git-ignored bulk data) are downloaded here
CCD_LOCAL = ROOT / "esr_analysis/resources/data/ccd"
CCD_CACHE = OUT / "ccd"            # CCD files missing locally are downloaded here
REPORT = OUT / "rebuild_report.csv"
REPORT_COLUMNS = ["target", "frozen_path", "rebuilt_path", "script", "status", "detail"]
TOL = 1e-6

OUT.mkdir(parents=True, exist_ok=True)


def read_ids():
    ids = []
    for name in ["data/pdb_ids_esr1.txt", "data/pdb_ids_esr2.txt"]:
        ids += [line.strip().upper() for line in (ROOT / name).read_text().splitlines() if line.strip()]
    return list(dict.fromkeys(ids))


def esr2_ids():
    return [line.strip().upper() for line in (ROOT / "data/pdb_ids_esr2.txt").read_text().splitlines() if line.strip()]


def raw_path(pdb_id):
    local = RAW_LOCAL / f"{pdb_id}.cif"
    return local if local.exists() else RAW_CACHE / f"{pdb_id}.cif"


def ccd_path(code):
    local = CCD_LOCAL / f"{code}.cif"
    return local if local.exists() else CCD_CACHE / f"{code}.cif"


def read_csv(path):
    """Read text-safe: 'NA' is the sodium CCD code, only empty fields are missing."""
    return pd.read_csv(path, keep_default_na=False, na_values=[""], low_memory=False)


def rel(path):
    return Path(path).resolve().relative_to(ROOT).as_posix()


# ---------------------------------------------------------------- comparison
def compare_frames(frozen, rebuilt, key, ignore=()):
    """Compare two tables keyed by `key` (str or list). Return (status, detail).

    `ignore` columns are reported but not counted (e.g. a local file path). Embedded CRLF in text
    is normalised to LF before comparing and the number of such cells is reported."""
    keys = [key] if isinstance(key, str) else list(key)
    notes = [f"rows frozen={len(frozen)} rebuilt={len(rebuilt)}"]
    fz = frozen.copy()
    rb = rebuilt.copy()
    for k in keys:
        fz[k] = fz[k].astype(str)
        rb[k] = rb[k].astype(str)
    fk = set(map(tuple, fz[keys].values))
    rk = set(map(tuple, rb[keys].values))
    only_f, only_r = fk - rk, rk - fk
    if only_f or only_r:
        notes.append(f"keys only_frozen={len(only_f)} only_rebuilt={len(only_r)}"
                     f" eg_frozen={sorted(only_f)[:3]} eg_rebuilt={sorted(only_r)[:3]}")
    missing_cols = [c for c in fz.columns if c not in rb.columns]
    extra_cols = [c for c in rb.columns if c not in fz.columns]
    if missing_cols:
        notes.append(f"columns missing in rebuild={missing_cols}")
    if extra_cols:
        notes.append(f"extra columns in rebuild={extra_cols}")
    merged = fz.merge(rb, on=keys, how="inner", suffixes=("__f", "__r"))
    exact_all, approx_only, diffs = True, [], []
    crlf = 0
    for col in [c for c in fz.columns if c in rb.columns and c not in keys]:
        a, b = merged[f"{col}__f"], merged[f"{col}__r"]
        is_bool = pd.api.types.is_bool_dtype(a) or pd.api.types.is_bool_dtype(b)
        an = pd.to_numeric(a, errors="coerce").astype(float) if not is_bool else a
        bn = pd.to_numeric(b, errors="coerce").astype(float) if not is_bool else b
        numeric = (not is_bool and (an.notna() == a.notna()).all() and (bn.notna() == b.notna()).all()
                   and a.notna().any())
        if numeric:
            both_nan = an.isna() & bn.isna()
            delta = (an - bn).abs()
            bad = ~both_nan & ~(delta <= TOL)
            if bad.any():
                exact_all = False
                diffs.append(f"{col}: {int(bad.sum())} differ (max |d|={np.nanmax(delta[~both_nan]) if (~both_nan).any() else 'nan'})")
            elif (delta[~both_nan] > 0).any():
                approx_only.append(f"{col} max|d|={delta[~both_nan].max():.2e}")
        else:
            sa = a.fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
            sb = b.fillna("").astype(str).str.replace(r"\.0$", "", regex=True)
            # booleans written as True/False vs true/false
            sa, sb = sa.str.lower().where(sa.str.lower().isin(["true", "false"]), sa), sb.str.lower().where(sb.str.lower().isin(["true", "false"]), sb)
            crlf_cells = (sa != sb) & (sa.str.replace("\r\n", "\n") == sb.str.replace("\r\n", "\n"))
            crlf += int(crlf_cells.sum())
            bad = (sa != sb) & ~crlf_cells
            if bad.any() and col in ignore:
                notes.append(f"ignored column {col}: {int(bad.sum())} differ")
            elif bad.any():
                exact_all = False
                ex = merged.loc[bad, keys + [f"{col}__f", f"{col}__r"]].head(1).values.tolist()
                diffs.append(f"{col}: {int(bad.sum())} differ eg={str(ex)[:160]}")
    if crlf:
        notes.append(f"{crlf} text cells differ only by embedded CRLF vs LF")
    if diffs:
        notes.append("value diffs: " + " | ".join(diffs))
    if approx_only:
        notes.append("within tol: " + ", ".join(approx_only))
    order_same = list(map(tuple, fz[keys].values)) == list(map(tuple, rb[keys].values))
    if not order_same and not (only_f or only_r):
        notes.append("row order differs")
    if exact_all and not only_f and not only_r and not missing_cols and len(fz) == len(rb):
        cosmetic = crlf or any(n.startswith("ignored column") for n in notes)
        status = "identical" if not approx_only and order_same and not cosmetic else "equivalent"
    else:
        status = "differs"
    return status, "; ".join(notes)


def compare_arrays(frozen, rebuilt):
    if frozen.shape != rebuilt.shape:
        return "differs", f"shape frozen={frozen.shape} rebuilt={rebuilt.shape}"
    d = np.abs(frozen.astype(float) - rebuilt.astype(float))
    m = float(np.nanmax(d)) if d.size else 0.0
    nan_mismatch = int((np.isnan(frozen.astype(float)) != np.isnan(rebuilt.astype(float))).sum())
    detail = f"shape={frozen.shape}; dtype frozen={frozen.dtype} rebuilt={rebuilt.dtype}; max|d|={m:.3g}; elements differing>{TOL}={int((d > TOL).sum())}; nan mismatch={nan_mismatch}"
    if m == 0 and nan_mismatch == 0:
        return ("identical" if frozen.dtype == rebuilt.dtype else "equivalent"), detail
    if m <= TOL and nan_mismatch == 0:
        return "equivalent", detail
    return "differs", detail


def report(target, frozen_path, rebuilt_path, script, status, detail):
    """Upsert one row (by target) in results/rebuild/rebuild_report.csv."""
    row = {"target": target, "frozen_path": frozen_path, "rebuilt_path": rebuilt_path,
           "script": script, "status": status, "detail": detail}
    table = pd.read_csv(REPORT, keep_default_na=False) if REPORT.exists() else pd.DataFrame(columns=REPORT_COLUMNS)
    table = table[table.target != target]
    table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
    table.to_csv(REPORT, index=False)
    print(f"[{status}] {target}: {detail[:300]}")
