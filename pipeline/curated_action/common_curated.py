"""Shared helpers, styling, statistics and export utilities for Stage 12 (curated action analysis)."""
from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[2]))
DATA = ROOT / "data"
RESULTS_BASE = ROOT / "results"
OUT_ROOT = RESULTS_BASE / "12_curated_action"
TABLES_DIR = OUT_ROOT / "tables"
FIGS_DIR = OUT_ROOT / "figures"
PYMOL_DIR = OUT_ROOT / "pymol_ligands"
DOCS_DIR = ROOT / "docs"

SEED = 42

# Project structural group colors (from pipeline/figures/style.py)
CL7, CL4, MINOR = "#2a78d6", "#eb6834", "#a3a29c"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#f1f0ec"
GROUP_ORDER = ["cl7", "cl4", "minor"]
GROUP_COLOR = {"cl7": CL7, "cl4": CL4, "minor": MINOR, "multiple": "#4a3aa7"}
GROUP_MARKER = {"cl7": "o", "cl4": "s", "minor": "^", "multiple": "D"}
GROUP_LABEL = {
    "cl7": "cl7 (H12 over pocket)",
    "cl4": "cl4 (H12 in groove)",
    "minor": "minor (5 small clusters)",
    "multiple": "multiple groups",
}

# Distinct colorblind-safe palette for curated functional classes (does not collide with cl7/cl4)
CLASS_ORDER_CORE = ["agonist", "partial_agonist", "antagonist"]
CLASS_ORDER_ALL = ["agonist", "partial_agonist", "antagonist", "unknown", "mixed"]
CLASS_COLOR = {
    "agonist": "#1b9e77",          # teal green
    "partial_agonist": "#7b3294",  # royal purple
    "antagonist": "#c2185b",       # crimson magenta
    "unknown": "#8a8984",          # neutral grey
    "mixed": "#5d4037",            # brown
    "unlabelled": "#bdbcb6",
    "SERM": "#d95f02",
    "SERD": "#e7298a",
}
CLASS_MARKER = {
    "agonist": "o",
    "partial_agonist": "s",
    "antagonist": "D",
    "unknown": "^",
    "mixed": "P",
}
CLASS_LABEL = {
    "agonist": "agonist",
    "partial_agonist": "partial_agonist",
    "antagonist": "antagonist",
    "unknown": "unknown",
    "mixed": "mixed",
}

# Ordinary variation thresholds from Stage 03 (dist_H12c_H11axis_A)
ORDINARY_THR_A = {"cl7": 0.91, "cl4": 1.08}


def ensure_dirs():
    for d in (OUT_ROOT, TABLES_DIR, FIGS_DIR, PYMOL_DIR, DOCS_DIR):
        d.mkdir(parents=True, exist_ok=True)


def setup_mpl():
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 10.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": 9.5,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.8,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "lines.linewidth": 2.0,
        "svg.fonttype": "none",
        "axes.grid": False,
    })


def read_csv_safe(path: Path | str, **kw) -> pd.DataFrame:
    """Read CSV preserving literal 'NA' (e.g., sodium CCD code)."""
    return pd.read_csv(path, keep_default_na=False, na_values=[""], **kw)


def group_of(cluster) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(cluster), "minor")


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def wrap_title(fig_or_ax, text: str, width: int = 92, y: float = 0.99):
    t = "\n".join(textwrap.wrap(text, width))
    if hasattr(fig_or_ax, "suptitle"):
        fig_or_ax.suptitle(t, x=0.01, y=y, ha="left", va="top", fontsize=11, fontweight="bold", color=INK)
    else:
        fig_or_ax.set_title(t, color=INK)


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


def wilson_ci(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt((p * (1 - p) + z * z / (4 * n)) / n) / denom
    return (max(0.0, center - half), min(1.0, center + half))


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


def save_fig_bundle(
    fig,
    fig_id: str,
    data_df: pd.DataFrame | None,
    caption: str,
    manifest_rows: list[dict],
    question: str,
    sources: str,
    placement: str,
    unit: str,
):
    """Save PNG + SVG + CSV data + caption .txt and append to manifest."""
    ensure_dirs()
    png = FIGS_DIR / f"{fig_id}.png"
    svg = FIGS_DIR / f"{fig_id}.svg"
    fig.savefig(png, dpi=220, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    csv_name = ""
    if data_df is not None:
        csv_path = FIGS_DIR / f"{fig_id}_data.csv"
        data_df.to_csv(csv_path, index=False)
        csv_name = str(csv_path.relative_to(ROOT))
    cap_path = FIGS_DIR / f"{fig_id}_caption.txt"
    cap_path.write_text(caption.strip() + "\n", encoding="utf-8")
    manifest_rows.append({
        "fig_id": fig_id,
        "png": str(png.relative_to(ROOT)),
        "svg": str(svg.relative_to(ROOT)),
        "data_csv": csv_name,
        "caption_txt": str(cap_path.relative_to(ROOT)),
        "question": question,
        "sources": sources,
        "unit": unit,
        "placement": placement,
        "caption": caption.strip(),
    })
    print(f"[fig] {fig_id} -> {png.name}, {svg.name}")
