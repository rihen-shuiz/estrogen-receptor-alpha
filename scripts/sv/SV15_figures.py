"""SV15 figures, drawn from saved CSVs only (rule 1.8).

SV15_dpx_by_group: dpX per ligand (primary set: full ChEMBL match; E2 = the 0 line, not a point), strip + box,
  (a) by SV10 ligand group, (b) by SV12 action class. Filled markers = 'relative' (median same-assay dpX, the only
  values in the boxes and tests); hollow markers = 'absolute' (no same-assay E2 pair: median binding pChEMBL minus the
  E2 median of the same type family, column dpx_plot; drawing only). Marker shape and colour = action class
  (scripts/sv/common/class_style.py, unlabelled grey); boxes of panel (a) in the cluster colours.
SV15_dpx_vs_filling: dpX vs ligand median filling % (SV14 primary-analysis chains), same marker rules; Spearman rho
  from tests.csv (relative ligands only).
Sources: results/sv/SV15/affinity_per_ligand.csv, tests.csv, numbers.json (class gate counts). Outputs: figures/sv/SV15/*.png + .svg; panels
figures/panels/sv/SV15/.
Run: .venv/bin/python -B scripts/sv/SV15_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.class_style import CLASS_COLOR, CLASS_LABEL, CLASS_MARKER, CLASS_ORDER  # noqa: E402
from common.sv_fig import CLUSTER_COLOR, GREY, INK2, MUTED, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV15"
YLABEL = "ΔpX = pX(ligand) − pX(E2), same binding assay\n(0 = E2, RBA 100; +1 = 10× tighter, RBA 1000)"


def trow(tests, comparison, level, set_name="primary"):
    r = tests[(tests.set == set_name) & (tests.comparison == comparison) & (tests.level == level)]
    return r.iloc[0] if len(r) else None


def draw_points(ax, x, s, rng, size=30):
    """s: rows of affinity_per_ligand at one x position; class style, filled = relative, hollow = absolute."""
    for r in s.itertuples():
        cls = r.action_class if r.action_class in CLASS_COLOR else "unlabelled"
        col, mk = CLASS_COLOR[cls], CLASS_MARKER[cls]
        xx = x + rng.uniform(-0.2, 0.2)
        if r.affinity_basis == "relative":
            ax.scatter(xx, r.dpx_plot, s=size, marker=mk, facecolors=col, edgecolors="white", linewidths=0.5,
                       zorder=4, alpha=0.95)
        else:
            ax.scatter(xx, r.dpx_plot, s=size * 0.8, marker=mk, facecolors="none", edgecolors=col, linewidths=0.9,
                       zorder=3, alpha=0.8)


def box(ax, pos, vals, color):
    if len(vals) < 3:
        return
    bp = ax.boxplot([vals], positions=[pos], widths=0.55, showfliers=False, patch_artist=True,
                    medianprops={"color": INK2, "lw": 1.3}, whiskerprops={"color": MUTED}, capprops={"color": MUTED},
                    boxprops={"edgecolor": MUTED})
    bp["boxes"][0].set_facecolor(color)
    bp["boxes"][0].set_alpha(0.22)


def legend_handles(present):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker=CLASS_MARKER[c], markersize=6, markerfacecolor=CLASS_COLOR[c],
                markeredgecolor="white", label=CLASS_LABEL[c]) for c in CLASS_ORDER if c in present]
    h += [Line2D([], [], ls="", marker="o", markersize=6, markerfacecolor=INK2, markeredgecolor="white",
                 label="filled: relative (same-assay E2 pair)"),
          Line2D([], [], ls="", marker="o", markersize=6, markerfacecolor="none", markeredgecolor=INK2,
                 label="hollow: absolute (no E2 pair;\nnot in boxes or tests)")]
    return h


def by_group_fig(plt, t, tests, gate):
    a = t[t.primary_set & ~t.is_e2_reference & t.dpx_plot.notna()].copy()
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.6, 5.2), gridspec_kw={"width_ratios": [3, 4]}, sharey=True)
    rng = np.random.default_rng(42)
    # (a) SV10 ligand group
    groups = ["cl7", "cl4", "other"]
    a["g"] = np.where(a.sv10_group.isin(["cl7", "cl4"]) & a.sv10_compared, a.sv10_group, "other")
    labels = []
    for i, g in enumerate(groups):
        s = a[a.g == g]
        rel = s[s.affinity_basis == "relative"]
        box(a1, i, rel.dpx_plot.to_numpy(), CLUSTER_COLOR.get(g, GREY))
        draw_points(a1, i, s, rng)
        name = {"cl7": "cl7", "cl4": "cl4", "other": "mixed / minor\n(not tested)"}[g]
        labels.append(f"{name}\n{len(rel)} rel. / {len(s) - len(rel)} abs.")
    a1.set_xticks(range(len(groups)), labels, fontsize=7.5)
    a1.axhline(0, color=INK2, lw=0.9, ls="--", zorder=1)
    a1.set_ylabel(YLABEL, fontsize=8)
    rl, re_ = trow(tests, "dpX by SV10 ligand group", "ligand"), trow(tests, "dpX by SV10 ligand group", "entry")
    a1.set_title(f"By SV10 ligand group, primary set\nmedian ΔpX cl7 {rl.median_1:.2f} (n {int(rl.n_1)}) vs cl4 "
                 f"{rl.median_2:.2f} (n {int(rl.n_2)})\nMann-Whitney ligands q {rl.q_bh:.2f}; entries "
                 f"({int(re_.n_1)} vs {int(re_.n_2)}) q {re_.q_bh:.2f}", fontsize=8.5)
    # (b) action class
    present = [c for c in CLASS_ORDER if (a.action_class == c).any()]
    cl = trow(tests, "dpX by action class", "ligand")
    tested = bool(cl is not None and cl.tested)
    labels = []
    for i, c in enumerate(present):
        s = a[a.action_class == c]
        rel = s[s.affinity_basis == "relative"]
        box(a2, i, rel.dpx_plot.to_numpy(), CLASS_COLOR[c])
        draw_points(a2, i, s, rng)
        labels.append(f"{CLASS_LABEL[c].replace(' (no class)', '')}\n{len(rel)} rel. / {len(s) - len(rel)} abs.")
    a2.set_xticks(range(len(present)), labels, fontsize=7.5)
    a2.axhline(0, color=INK2, lw=0.9, ls="--", zorder=1)
    a2.text(len(present) - 0.55, 0.05, "E2", fontsize=7, color=INK2, ha="right", va="bottom")
    if tested:
        title = f"By action class (SV12), primary set\n{cl.groups}, q {cl.q_bh:.2f}"
    else:
        ent = gate["entries_in_sv15_tested_set"]
        title = ("By action class (SV12), primary set\nnot tested: coverage gate (>= 5 entries per side)\n"
                 "entries with relative ΔpX: " + ", ".join(f"{c} {n}" for c, n in ent.items()))
    a2.set_title(title, fontsize=8.5)
    a2.legend(handles=legend_handles(set(a.action_class)), loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7,
              frameon=False, title="action class, SV12 (shape, colour)", title_fontsize=7, alignment="left")
    fig.tight_layout()
    return fig, a1, a2


def vs_filling_fig(plt, t, tests):
    a = t[t.primary_set & ~t.is_e2_reference & t.dpx_plot.notna() & t.filling_pct_median.notna()]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    for r in a.itertuples():
        cls = r.action_class if r.action_class in CLASS_COLOR else "unlabelled"
        col, mk = CLASS_COLOR[cls], CLASS_MARKER[cls]
        if r.affinity_basis == "relative":
            ax.scatter(r.filling_pct_median, r.dpx_plot, s=34, marker=mk, facecolors=col, edgecolors="white",
                       linewidths=0.5, zorder=4)
            ax.annotate(r.ccd, (r.filling_pct_median, r.dpx_plot), xytext=(3, 2), textcoords="offset points",
                        fontsize=6, color=INK2)
        else:
            ax.scatter(r.filling_pct_median, r.dpx_plot, s=24, marker=mk, facecolors="none", edgecolors=col,
                       linewidths=0.8, alpha=0.75, zorder=3)
    ax.axhline(0, color=INK2, lw=0.9, ls="--", zorder=1)
    e2 = t[t.is_e2_reference].iloc[0]
    ax.axvline(e2.filling_pct_median, color=MUTED, lw=0.8, ls=":", zorder=1)
    ax.text(e2.filling_pct_median, ax.get_ylim()[1], " E2 filling", fontsize=6.5, color=MUTED, va="top")
    ax.set_xlabel("pocket filling % (SV14; median over the ligand's primary-analysis chains)")
    ax.set_ylabel(YLABEL, fontsize=8)
    rl, re_ = trow(tests, "dpX vs filling % (SV14)", "ligand"), trow(tests, "dpX vs filling % (SV14)", "entry")
    n_rel = int((a.affinity_basis == "relative").sum())
    ax.set_title(f"ΔpX vs pocket filling: Spearman ρ ligands {rl.effect:.2f} (q {rl.q_bh:.2f}, n {int(rl.n_1)}), "
                 f"entries {re_.effect:.2f} (q {re_.q_bh:.2f}, n {int(re_.n_1)})\n{n_rel} relative (filled, "
                 f"labelled) and {len(a) - n_rel} absolute ligands (hollow, not in ρ)", fontsize=8.5)
    ax.legend(handles=legend_handles(set(a.action_class)), loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=7,
              frameon=False, title="action class, SV12 (shape, colour)", title_fontsize=7, alignment="left")
    fig.tight_layout()
    return fig, ax


def main():
    t = pd.read_csv(RES / "affinity_per_ligand.csv")
    for c in ("primary_set", "is_e2_reference", "sv10_compared"):
        t[c] = t[c].astype(bool)
    tests = pd.read_csv(RES / "tests.csv")
    gate = json.loads((RES / "numbers.json").read_text())["class_gate_primary"]
    plt = mpl_setup()
    fig, a1, a2 = by_group_fig(plt, t, tests, gate)
    paths = save_figure(fig, "SV15", "SV15_dpx_by_group", panels={"by_cluster": a1, "by_class": a2})
    fig, ax = vs_filling_fig(plt, t, tests)
    paths += save_figure(fig, "SV15", "SV15_dpx_vs_filling", panels={"scatter": ax})
    for p in paths:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
