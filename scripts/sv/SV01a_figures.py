"""SV01a figures, drawn only from the saved SV01a CSVs (rule 1.8).

1. SV01a_forest_completeness: odds ratios (95 % CI) of the full Q+X+B model for 'H12 complete' (538-547), predictors
   grouped and coloured by block, panels cl7 and cl4, chain level = block-coloured dots, entry level = grey dots;
   beside it, McFadden gain of each block entered first and entered last (block-coloured bars = chain level, grey
   bars = entry level; solid = entered first, hatched = entered last).
   Source: results/sv/SV01a/odds_ratios.csv, model_blocks.csv, verdict_table.csv.
2. SV01a_completeness_profile: fraction of chains with the C-alpha observed per residue 520-550, cl7 vs cl4
   (solid = chains, dashed = entry-weighted). Source: results/sv/SV01a/completeness_profile.csv.
Run: .venv/bin/python -I scripts/sv/SV01a_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.sv_fig import (BLUE, GREY, INK, INK2, MAGENTA, MUTED, ORANGE, VIOLET, YELLOW,  # noqa: E402
                           mpl_setup, save_figure)

RES = ROOT / "results/sv/SV01a"
STEP = "SV01a"
BLOCK_COLOR = {"Q": VIOLET, "X": YELLOW, "B": MAGENTA, "Q+X": INK2}
BLOCK_NAME = {"Q": "Q quality", "X": "X experiment", "B": "B biology"}
LIM = (1e-3, 1e3)
H12 = (538, 547)
SHORT = {"resolution_A": "resolution (A) *", "r_free": "R-free *", "mean_b": "mean B, entry *", "ph": "growth pH *",
         "growth_temp_K": "growth T *", "ambient_temp_K": "collection T *", "p21": "space group P21",
         "esr1_chains_asu": "ESR1 chains in ASU *", "dimer_symmetry": "partner by symmetry",
         "ph_missing": "pH not reported", "growth_temp_K_missing": "growth T not reported",
         "ambient_temp_K_missing": "collection T not reported", "ligand_heavy_atoms": "ligand heavy atoms *"}


def short_reason(status: str) -> str:
    s = status.replace("dropped: ", "")
    if s.startswith("constant"):
        return "constant"
    if s.startswith("collinear"):
        return "collinear"
    if s.startswith("joint"):
        return "joint separation"
    if s.startswith("separation"):
        return "separates"
    return s


def forest(ax, o, cluster, outcome="complete"):
    d = o[(o.outcome == outcome) & (o.cluster == cluster)]
    ch = d[d.level == "chain"].set_index("predictor")
    en = d[d.level == "entry"].set_index("predictor")
    preds = list(ch.index)  # block order Q, X, B as written by the analysis script
    n = len(preds)
    ys = np.arange(n)[::-1]
    for blk in ("Q", "X", "B"):
        idx = [i for i, p in enumerate(preds) if ch.loc[p, "block"] == blk]
        if idx:
            ax.axhspan(ys[idx[-1]] - 0.5, ys[idx[0]] + 0.5, color=BLOCK_COLOR[blk], alpha=0.07, lw=0)
    for y, p in zip(ys, preds):
        for lvl, tab, dy, col in (("chain", ch, 0.17, BLOCK_COLOR[ch.loc[p, "block"]]), ("entry", en, -0.17, GREY)):
            r = tab.loc[p]
            if r.status != "fitted":
                continue
            lo, hi, est = (np.clip(v, *LIM) for v in (r.ci_low, r.ci_high, r.odds_ratio))
            ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=1.5, solid_capstyle="butt", zorder=2)
            for edge, raw, mk in ((lo, r.ci_low, "<"), (hi, r.ci_high, ">")):
                if raw < LIM[0] or raw > LIM[1]:
                    ax.plot(edge, y + dy, marker=mk, color=col, ms=4, zorder=3)
            clipped = r.odds_ratio < LIM[0] or r.odds_ratio > LIM[1]
            ax.scatter(est, y + dy, s=30 if lvl == "chain" else 24, color="white" if clipped else col,
                       edgecolor=col if clipped else "white", lw=1.0 if clipped else 0.5, zorder=4)
        rc, re_ = ch.loc[p, "status"], en.loc[p, "status"]
        if rc != "fitted" or re_ != "fitted":
            a = short_reason(rc) if rc != "fitted" else None
            b = short_reason(re_) if re_ != "fitted" else None
            txt = (f"dropped: {a}" if a == b else
                   "; ".join(x for x in (f"chain: {a}" if a else "", f"entry: {b}" if b else "") if x))
            ax.text(LIM[1] * 0.8, y, txt, ha="right", va="center", fontsize=6.5, color=MUTED)
    ax.axvline(1, color=INK, lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(*LIM)
    ax.set_ylim(-0.6, n - 0.4)
    labels = [SHORT.get(p, ch.loc[p, "label"]) for p in preds]
    ax.set_yticks(ys, labels, fontsize=7.5)
    ax.set_xlabel("odds ratio for complete H12 (95% CI), full Q+X+B model;\ncontinuous predictors (*) per 1 SD")
    n_ch, inc_ch = int(ch.n.iloc[0]), int(ch.n_incomplete.iloc[0])
    n_en, inc_en = int(en.n.iloc[0]), int(en.n_incomplete.iloc[0])
    ax.set_title(f"cl{cluster}: {n_ch} chains ({inc_ch} incomplete H12), {n_en} entries ({inc_en})", fontsize=9.5)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=BLOCK_COLOR[b], label=f"{BLOCK_NAME[b]} (chain)") for b in "QXB"]
    h.append(Line2D([], [], marker="o", ls="", color=GREY, label="entry level"))
    h.append(Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, label=f"estimate beyond {LIM[0]:g}-{LIM[1]:g}"))
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, fontsize=6.5, handletextpad=0.3)


def gains(ax, mb, cluster, outcome="complete"):
    d = mb[(mb.outcome == outcome) & (mb.cluster == cluster)]
    blocks = ["Q", "X", "B", "Q+X"]
    offs = {("chain", "first"): 0.33, ("chain", "last"): 0.11, ("entry", "first"): -0.11, ("entry", "last"): -0.33}
    y0 = np.arange(len(blocks))[::-1] * 1.0
    for i, blk in enumerate(blocks):
        for (lvl, order), dy in offs.items():
            v = d[(d.level == lvl) & (d.block == blk) & (d.order == order)].mcfadden_gain.iloc[0]
            col = BLOCK_COLOR[blk] if lvl == "chain" else GREY
            ax.barh(y0[i] + dy, v, height=0.2, color=col, alpha=1.0 if order == "first" else 0.5,
                    hatch=None if order == "first" else "////", edgecolor="white", lw=0.4)
            ax.text(v + 0.008, y0[i] + dy, f"{v:.3f}", va="center", fontsize=6, color=INK2)
    ax.set_yticks(y0, blocks)
    ax.set_xlabel("McFadden R2 gain")
    ax.set_xlim(0, d.mcfadden_gain.max() * 1.2)
    n_ch = int(d[d.level == "chain"].n.iloc[0])
    n_en = int(d[d.level == "entry"].n.iloc[0])
    ax.set_title(f"cl{cluster}: block gains ({n_ch} chains, {n_en} entries)", fontsize=9.5)
    from matplotlib.patches import Patch
    h = [Patch(color=INK2, label="chain, entered first"),
         Patch(facecolor=INK2, alpha=0.5, hatch="////", edgecolor="white", label="chain, entered last"),
         Patch(color=GREY, label="entry, entered first"),
         Patch(facecolor=GREY, alpha=0.5, hatch="////", edgecolor="white", label="entry, entered last")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=6.5)


def forest_figure():
    o = pd.read_csv(RES / "odds_ratios.csv")
    mb = pd.read_csv(RES / "model_blocks.csv")
    vt = pd.read_csv(RES / "verdict_table.csv")
    plt = mpl_setup()
    fig = plt.figure(figsize=(17, 9.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 0.8], hspace=0.6, wspace=0.42)
    ax7 = fig.add_subplot(gs[:, 0])
    ax4 = fig.add_subplot(gs[:, 1])
    g7 = fig.add_subplot(gs[0, 2])
    g4 = fig.add_subplot(gs[1, 2])
    forest(ax7, o, 7)
    forest(ax4, o, 4)
    gains(g7, mb, 7)
    gains(g4, mb, 4)
    v = vt[(vt.outcome == "complete") & (vt.order == "first")].set_index(["level", "cluster"])
    lv = vt[vt.outcome == "complete"].groupby("level").level_verdict.first()
    step = vt[vt.outcome == "complete"].step_verdict.iloc[0]

    def pair(level, c):
        r = v.loc[(level, c)]
        return f"{r.gain_QX:.3f} vs {r.gain_B:.3f}"
    fig.suptitle(
        f"H12 complete (538-547): Q+X vs B gain - chain level cl7 {pair('chain', 7)}, cl4 {pair('chain', 4)} "
        f"({lv['chain']}); entry level cl7 {pair('entry', 7)}, cl4 {pair('entry', 4)} ({lv['entry']}); "
        f"[DEFAULT] verdict: {step}", x=0.01, ha="left", fontweight="bold", fontsize=10.5)
    return save_figure(fig, STEP, "SV01a_forest_completeness",
                       {"forest_cl7": ax7, "forest_cl4": ax4, "gains_cl7": g7, "gains_cl4": g4})


def profile_figure():
    p = pd.read_csv(RES / "completeness_profile.csv")
    p = p[(p.group_type == "cluster") & (p.residue >= 520) & (p.residue <= 550)]
    plt = mpl_setup()
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.axvspan(H12[0] - 0.5, H12[1] + 0.5, color=GREY, alpha=0.15, lw=0)
    ax.text((H12[0] + H12[1]) / 2, 0.03, f"H12 window {H12[0]}-{H12[1]}", ha="center", fontsize=7.5, color=INK2)
    summ = {}
    for g, col in (("cl7", BLUE), ("cl4", ORANGE)):
        d = p[p.group == g].sort_values("residue")
        ax.plot(d.residue, d.frac_observed, color=col, lw=2, label=f"{g} chains (n={int(d.n_chains.iloc[0])})")
        ax.plot(d.residue, d.frac_observed_entry, color=col, lw=1.2, ls="--",
                label=f"{g} entry-weighted (n={int(d.n_entries.iloc[0])})")
        summ[g] = d.set_index("residue").frac_observed
    ax.set_xlim(519.5, 550.5)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("ESR1 residue (UniProt)")
    ax.set_ylabel("fraction with C-alpha observed")
    ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title(f"C-alpha observed, residues 520-550: at {H12[1]} cl7 {summ['cl7'][H12[1]]:.2f} vs cl4 "
                 f"{summ['cl4'][H12[1]]:.2f}; at {H12[0]} cl7 {summ['cl7'][H12[0]]:.2f} vs cl4 {summ['cl4'][H12[0]]:.2f}",
                 fontsize=9.5)
    return save_figure(fig, STEP, "SV01a_completeness_profile", {"profile": ax})


def main():
    for pth in forest_figure() + profile_figure():
        print(pth.relative_to(ROOT), f"{pth.stat().st_size / 1e3:.1f} kB")


if __name__ == "__main__":
    main()
