"""Shared paths, readers and statistics for stage 12 (curated ligand action classes)."""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[2]))
DATA = ROOT / "data"
RESULTS_BASE = ROOT / "results"
OUT_ROOT = RESULTS_BASE / "12_curated_action"
TABLES_DIR = OUT_ROOT / "tables"

SEED = 42

CLASS_ORDER_CORE = ["agonist", "partial_agonist", "antagonist"]
# Ordinary variation thresholds from Stage 03 (dist_H12c_H11axis_A)
ORDINARY_THR_A = {"cl7": 0.91, "cl4": 1.08}


def ensure_dirs():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)


def read_csv_safe(path: Path | str, **kw) -> pd.DataFrame:
    """Read CSV preserving literal 'NA' (e.g., sodium CCD code)."""
    return pd.read_csv(path, keep_default_na=False, na_values=[""], **kw)


def group_of(cluster) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(cluster), "minor")


def bh_fdr(pvals) -> np.ndarray:
    pv = np.asarray(pvals, float)
    q = np.full_like(pv, np.nan)
    ok = np.isfinite(pv)
    x = pv[ok]
    if len(x):
        o = np.argsort(x)
        r = x[o] * len(x) / (np.arange(len(x)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty_like(x)
        qq[o] = np.minimum(r, 1.0)
        q[ok] = qq
    return q


def doi_block_id(df: pd.DataFrame) -> pd.Series:
    """Cluster ID for block bootstrap: DOI when present, otherwise 'NO_DOI:<pdb_id>'."""
    doi = df["doi"].fillna("").astype(str).str.strip()
    pdb = df["pdb_norm"].astype(str) if "pdb_norm" in df.columns else df["entry_id"].astype(str)
    return np.where(doi != "", doi, "NO_DOI:" + pdb)


def block_bootstrap_diff_medians(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    val_col: str,
    block_col: str = "doi_block",
    n_boot: int = 2000,
    seed: int = SEED,
    min_blocks: int = 4,
) -> dict:
    """Difference of medians (A - B) with DOI-block bootstrap 95% CI and Mann-Whitney p-value."""
    a = df_a.dropna(subset=[val_col]).copy()
    b = df_b.dropna(subset=[val_col]).copy()
    va = a[val_col].to_numpy(float)
    vb = b[val_col].to_numpy(float)
    if len(va) == 0 or len(vb) == 0:
        return dict(
            n_a=len(va), n_b=len(vb), blocks_a=0, blocks_b=0,
            med_a=np.nan, med_b=np.nan, diff=np.nan, ci_lo=np.nan, ci_hi=np.nan,
            p_mwu=np.nan, ci_valid=False,
        )
    med_a, med_b = float(np.median(va)), float(np.median(vb))
    diff = med_a - med_b
    p_mwu = float(mannwhitneyu(va, vb, alternative="two-sided").pvalue) if (len(va) >= 3 and len(vb) >= 3) else np.nan

    ba_groups = [g[val_col].to_numpy(float) for _, g in a.groupby(block_col)]
    bb_groups = [g[val_col].to_numpy(float) for _, g in b.groupby(block_col)]
    nb_a, nb_b = len(ba_groups), len(bb_groups)
    ci_lo, ci_hi, ci_valid = np.nan, np.nan, False
    if nb_a >= min_blocks and nb_b >= min_blocks:
        rng = np.random.default_rng(seed)
        boots = np.empty(n_boot, float)
        for i in range(n_boot):
            ia = rng.integers(0, nb_a, size=nb_a)
            ib = rng.integers(0, nb_b, size=nb_b)
            sa = np.concatenate([ba_groups[k] for k in ia])
            sb = np.concatenate([bb_groups[k] for k in ib])
            boots[i] = np.median(sa) - np.median(sb)
        ci_lo, ci_hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
        ci_valid = True
    return dict(
        n_a=int(len(va)), n_b=int(len(vb)),
        blocks_a=int(nb_a), blocks_b=int(nb_b),
        med_a=med_a, med_b=med_b, diff=diff,
        ci_lo=ci_lo, ci_hi=ci_hi, p_mwu=p_mwu, ci_valid=ci_valid,
    )


def hill_numbers(counts) -> dict:
    c = np.asarray(counts, float)
    c = c[c > 0]
    if len(c) == 0:
        return dict(q0=0, shannon_bits=np.nan, shannon_nats=np.nan, q1=np.nan, q2=np.nan)
    p = c / c.sum()
    h_nats = float(-(p * np.log(p)).sum())
    h_bits = float(-(p * np.log2(p)).sum())
    return dict(
        q0=int(len(p)),
        shannon_bits=h_bits,
        shannon_nats=h_nats,
        q1=float(np.exp(h_nats)),
        q2=float(1.0 / (p ** 2).sum()),
    )


def write_json(path: Path, data: dict):
    def conv(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating, float)):
            return None if not np.isfinite(v) else float(v)
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, dict):
            return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple, np.ndarray)):
            return [conv(x) for x in v]
        return v
    path.write_text(json.dumps(conv(data), indent=2, ensure_ascii=False))


