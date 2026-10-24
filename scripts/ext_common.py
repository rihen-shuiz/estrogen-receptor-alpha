"""Shared helpers for the extended analyses (tracker phase ext, tasks Y01...).

One place for: repo paths, loaders of the frozen inputs (788-chain RMSD matrix and its row order, S08 clusters,
QC tiers, H12 labels, B-factors, pockets, ligand of interest), the ONE structural-cluster palette used in every
figure, figure saving (200 dpi PNG + SVG) and the figure-card registry (results/ext/figure_cards.csv), from which
docs/figures/INDEX.md and docs/gallery.html are generated (scripts/ext_build_index.py).

Nothing here recomputes a frozen matrix; everything is loaded. Seeds: SEED = 42.
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[1]))
RES = ROOT / "results/ext"
FIG = ROOT / "figures/ext"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
SEED = 42
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
H12 = tuple(SETTINGS["h12_range"])

# ---------------------------------------------------------------- palette (validated: dataviz validate_palette.js)
# Structural clusters are the S08 average-linkage clusters at the primary 1.5 A cut. The two main clusters get
# slots 1-2; the five minor clusters (1, 2, 3, 5, 6; 14 chains in total) share slot 3 and are told apart by marker.
# Light-mode slots 1-3 pass the all-pairs CVD gate (worst dE 9.2), so they are safe on scatter plots.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
CATEGORICAL = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]
GREY = "#a3a29c"          # unknown / not_reported / other
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
SURFACE = "#fcfcfb"
CLUSTER_COLOR = {7: BLUE, 4: ORANGE, 1: AQUA, 2: AQUA, 3: AQUA, 5: AQUA, 6: AQUA}
CLUSTER_MARKER = {7: "o", 4: "o", 1: "^", 2: "s", 3: "D", 5: "v", 6: "P"}
MAIN = (7, 4)
SEQ_BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
DIVERGING = ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f2a9a8", "#e34948", "#a32a2a"]
TEXTURE_NOTE = "colour = S08 structural cluster at 1.5 A (blue cl7, orange cl4, aqua = minor clusters by marker)"


def cluster_name(c) -> str:
    c = int(c)
    return {7: "cl7 (main)", 4: "cl4 (main)"}.get(c, f"cl{c} (minor)")


def mpl_setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({
        "figure.dpi": 100, "savefig.dpi": 200, "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "legend.frameon": False,
        "figure.facecolor": "white", "axes.facecolor": "white", "svg.fonttype": "none",
        "axes.titleweight": "bold", "axes.titlelocation": "left", "lines.linewidth": 2,
    })
    for name, cols in (("esr_seq", SEQ_BLUE), ("esr_div", DIVERGING)):
        try:
            matplotlib.colormaps.register(LinearSegmentedColormap.from_list(name, cols))
        except ValueError:
            pass
    return plt


def save_fig(fig, name: str, *, section: str, task: str, shows: str, takeaway: str, source: str, script: str):
    """Write figures/ext/<name>.png (200 dpi) + .svg and register a figure card (replaces an older card of <name>)."""
    png, svg = FIG / f"{name}.png", FIG / f"{name}.svg"
    fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    import matplotlib.pyplot as plt
    plt.close(fig)
    register_card(name, section=section, task=task, shows=shows, takeaway=takeaway, source=source, script=script)
    return png


def register_card(name, *, section, task, shows, takeaway, source, script, file=None):
    path = RES / "figure_cards.csv"
    cols = ["figure", "file", "section", "task", "shows", "takeaway", "source", "script"]
    rows = []
    if path.exists():
        rows = [r for r in csv.DictReader(path.open(encoding="utf8")) if r["figure"] != name]
    rows.append(dict(figure=name, file=file or f"figures/ext/{name}.png", section=section, task=task, shows=shows,
                     takeaway=takeaway, source=source, script=script))
    with path.open("w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["task"], r["figure"])))


def write_numbers(task: str, values: dict):
    """Headline numbers per task -> results/ext/numbers_<task>.json (read by collect_numbers.py; never typed)."""
    def conv(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating,)):
            return float(v)
        if isinstance(v, dict):
            return {k: conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [conv(x) for x in v]
        return v
    (RES / f"numbers_{task}.json").write_text(json.dumps(conv(values), indent=1))


# ---------------------------------------------------------------- loaders
def rmsd():
    """S04 788 x 788 C-alpha RMSD (A) and its row order (row, chain_key, entry_id)."""
    m = np.load(ROOT / "results/rmsd_entry.npy")
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    assert m.shape == (len(idx), len(idx))
    return m, idx


def clusters():
    return pd.read_csv(ROOT / "results/entry_clusters.csv")


def chain_table() -> pd.DataFrame:
    """One row per matrix chain (788, matrix order) with the frozen labels joined. Nothing dropped."""
    _, idx = rmsd()
    t = idx.copy()
    cl = clusters()[["chain_key", "cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5", "qc_label"]]
    t = t.merge(cl, on="chain_key", how="left", validate="one_to_one")
    cs = pd.read_csv(ROOT / "results/chain_status.csv")
    keep = ["chain_key", "resolution_A", "r_free", "mean_b_A2", "tier_resolution", "tier_rfree", "tier_b",
            "h12_to_core_ratio", "b_group", "h12_b_group", "h12_rmsd", "h12_label", "core_rmsd", "core_rmsd_label",
            "mapped_lbd_ca", "mapped_core_ca"]
    t = t.merge(cs[keep], on="chain_key", how="left", validate="one_to_one")
    hb = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv")
    t = t.merge(hb[["chain_key", "mean_b_chain", "median_ca_b_chain", "core_mean_b", "h12_mean_b"]],
                on="chain_key", how="left", validate="one_to_one")
    loi = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv")
    loi = loi.drop_duplicates("chain_key")[["chain_key", "ligand_id", "status", "n_pocket_residues"]]
    t = t.merge(loi.rename(columns={"status": "ligand_status"}), on="chain_key", how="left", validate="one_to_one")
    pk = pd.read_csv(ROOT / "results/pockets.csv")
    t = t.merge(pk[["chain_key", "pocket_score", "sas_points", "n_residues"]].rename(
        columns={"n_residues": "pocket_n_residues"}), on="chain_key", how="left", validate="one_to_one")
    ent = pd.read_csv(ROOT / "results/entries.csv")[["entry_id", "space_group", "deposit_date"]]
    t = t.merge(ent, on="entry_id", how="left", validate="many_to_one")
    t["cluster"] = t["cut_1.5"].astype(int)
    t["tier"] = t["qc_label"].fillna("unknown")
    assert len(t) == len(idx) and (t.chain_key.values == idx.chain_key.values).all()
    return t


def coords():
    """S03 superposed C-alpha coordinates (788, 246, 3), NaN where unobserved, and the ESR1 positions (305-550)."""
    d = np.load(ROOT / "results/alignment/aligned_coordinates.npz", allow_pickle=True)
    keys = [str(k) for k in d["chain_keys"]]
    _, idx = rmsd()
    assert keys == idx.chain_key.tolist()
    return d["coords"].astype(float), d["positions"].astype(int)


def helices() -> pd.DataFrame:
    h = pd.read_csv(ROOT / "results/helix_table.csv")
    return h[h.helix_class == 1].reset_index(drop=True)


# Conventional ESR1 LBD helix names (Brzozowski et al. 1997 numbering) mapped onto the mmCIF HELX_P segments of 7NFB
# in results/helix_table.csv (helix_class 1). The table itself gives no names; these names are the textbook ones
# for the same residue ranges, used only as axis labels.
HELIX_NAMES = {"AB4": "H1", "AB5": "H3", "AB7": "H4/5", "AB8": "H6", "AB9": "H7", "AC1": "H8", "AC2": "H9",
               "AC3": "H10/11", "AC4": "H12"}


def helix_bands(ax, positions_axis="x", alpha=0.10, label=True, ymax_frac=1.0):
    """Shade the helices (results/helix_table.csv) along a residue axis and label them on top."""
    h = helices()
    for i, r in h.iterrows():
        col = RED if r.helix == "AC4" else BLUE
        if positions_axis == "x":
            ax.axvspan(r.start - 0.5, r.end + 0.5, color=col, alpha=alpha if r.helix != "AC4" else 0.18, lw=0)
            if label:
                ax.text((r.start + r.end) / 2, 1.0, HELIX_NAMES.get(r.helix, r.helix), transform=ax.get_xaxis_transform(),
                        ha="center", va="bottom", fontsize=7, color=INK2)
        else:
            ax.axhspan(r.start - 0.5, r.end + 0.5, color=col, alpha=alpha, lw=0)


def bh(p):
    """Benjamini-Hochberg q-values (NaN-safe)."""
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    q[ok] = out
    return q


def cramers_v(table) -> float:
    from scipy.stats import chi2_contingency
    t = np.asarray(table, float)
    t = t[t.sum(1) > 0][:, t.sum(0) > 0]
    if t.shape[0] < 2 or t.shape[1] < 2:
        return float("nan")
    chi2 = chi2_contingency(t, correction=False)[0]
    n = t.sum()
    return float(np.sqrt(chi2 / (n * (min(t.shape) - 1))))


def scatter_clusters(ax, x, y, cl, s=10, alpha=0.75, legend=True):
    cl = np.asarray(cl).astype(int)
    for c in sorted(set(cl), key=lambda c: (c not in MAIN, c)):
        m = cl == c
        ax.scatter(np.asarray(x)[m], np.asarray(y)[m], s=s if c in MAIN else s * 2.2, c=CLUSTER_COLOR[c],
                   marker=CLUSTER_MARKER[c], alpha=alpha if c in MAIN else 1.0, lw=0.3, edgecolor="white",
                   label=f"{cluster_name(c)} n={m.sum()}", zorder=2 if c in MAIN else 3)
    if legend:
        ax.legend(fontsize=7, markerscale=1.2, loc="best", handletextpad=0.3)
