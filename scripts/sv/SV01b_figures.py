"""SV01b figures, drawn only from the saved SV01b CSVs (rule 1.8). Plot functions of scripts/sv/SV01a_figures.py are reused.

1. SV01b_forest_completeness (redraw of the SV01a plot with the lattice terms in block X): odds ratios (95 % CI) of the
   full Q+X+B model for 'H12 complete' (538-547), predictors grouped and coloured by block, panels cl7 and cl4,
   chain level = block-coloured dots, entry level = grey dots; beside it the McFadden gain of each block entered first
   and entered last. Source: results/sv/SV01b/odds_ratios.csv, model_blocks.csv, verdict_table.csv.
2. SV01b_forest_mobility: the same layout for the mobility outcome (mean B z-score of the observed residues 528-547
   against the chain's core v3 C-alpha B): OLS coefficients (z units, 95 % CI) and OLS R2 gains.
Run: .venv/bin/python -B scripts/sv/SV01b_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV01a_figures as F  # noqa: E402
from common.sv_fig import GREY, INK, INK2, MUTED, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV01b"
STEP = "SV01b"
F.SHORT.update({"lattice_asu": "lattice, ASU partner *", "lattice_sym": "lattice, symmetry mate *"})


def load():
    o = pd.read_csv(RES / "odds_ratios.csv")
    mb = pd.read_csv(RES / "model_blocks.csv")
    mb["mcfadden_gain"] = mb["gain"]  # column name read by SV01a_figures.gains
    vt = pd.read_csv(RES / "verdict_table.csv")
    return o, mb, vt


def title(vt, outcome, what):
    v = vt[(vt.outcome == outcome) & (vt.order == "first")].set_index(["level", "cluster"])
    lv = vt[vt.outcome == outcome].groupby("level").level_verdict.first()
    ov = vt[vt.outcome == outcome].outcome_verdict.iloc[0]

    def pair(level, c):
        r = v.loc[(level, c)]
        return f"{r.gain_QX:.3f} vs {r.gain_B:.3f}"
    return (f"{what}: Q+X vs B gain (lattice terms in X) - chain level cl7 {pair('chain', 7)}, cl4 {pair('chain', 4)} "
            f"({lv['chain']}); entry level cl7 {pair('entry', 7)}, cl4 {pair('entry', 4)} ({lv['entry']}); "
            f"[DEFAULT] rule: {ov}")


def forest_lin(ax, o, cluster, lim, outcome="mobility"):
    d = o[(o.outcome == outcome) & (o.cluster == cluster)]
    ch = d[d.level == "chain"].set_index("predictor")
    en = d[d.level == "entry"].set_index("predictor")
    preds = list(ch.index)
    n = len(preds)
    ys = np.arange(n)[::-1]
    for blk in ("Q", "X", "B"):
        idx = [i for i, p in enumerate(preds) if ch.loc[p, "block"] == blk]
        if idx:
            ax.axhspan(ys[idx[-1]] - 0.5, ys[idx[0]] + 0.5, color=F.BLOCK_COLOR[blk], alpha=0.07, lw=0)
    for y, p in zip(ys, preds):
        for lvl, tab, dy, col in (("chain", ch, 0.17, F.BLOCK_COLOR[ch.loc[p, "block"]]), ("entry", en, -0.17, GREY)):
            r = tab.loc[p]
            if r.status != "fitted":
                continue
            lo, hi, est = (np.clip(v, -lim, lim) for v in (r.ci_low, r.ci_high, r.estimate))
            ax.plot([lo, hi], [y + dy, y + dy], color=col, lw=1.5, solid_capstyle="butt", zorder=2)
            for edge, raw, mk in ((lo, r.ci_low, "<"), (hi, r.ci_high, ">")):
                if abs(raw) > lim:
                    ax.plot(edge, y + dy, marker=mk, color=col, ms=4, zorder=3)
            clipped = abs(r.estimate) > lim
            ax.scatter(est, y + dy, s=30 if lvl == "chain" else 24, color="white" if clipped else col,
                       edgecolor=col if clipped else "white", lw=1.0 if clipped else 0.5, zorder=4)
        rc, re_ = ch.loc[p, "status"], en.loc[p, "status"]
        if rc != "fitted" or re_ != "fitted":
            a = F.short_reason(rc) if rc != "fitted" else None
            b = F.short_reason(re_) if re_ != "fitted" else None
            txt = (f"dropped: {a}" if a == b else
                   "; ".join(x for x in (f"chain: {a}" if a else "", f"entry: {b}" if b else "") if x))
            ax.text(lim * 0.97, y, txt, ha="right", va="center", fontsize=6.5, color=MUTED)
    ax.axvline(0, color=INK, lw=0.8, zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(ys, [F.SHORT.get(p, ch.loc[p, "label"]) for p in preds], fontsize=7.5)
    ax.set_xlabel("coefficient on mobility (mean B z-score, 528-547; 95% CI), full Q+X+B OLS;\n"
                  "continuous predictors (*) per 1 SD")
    ax.set_title(f"cl{cluster}: {int(ch.n.iloc[0])} chains, {int(en.n.iloc[0])} entries", fontsize=9.5)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], marker="o", ls="", color=F.BLOCK_COLOR[b], label=f"{F.BLOCK_NAME[b]} (chain)") for b in "QXB"]
    h.append(Line2D([], [], marker="o", ls="", color=GREY, label="entry level"))
    h.append(Line2D([], [], marker="o", ls="", mfc="white", mec=INK2, label=f"estimate beyond +-{lim:g}"))
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, fontsize=6.5, handletextpad=0.3)


def figure(outcome, name, what, linear):
    o, mb, vt = load()
    plt = mpl_setup()
    fig = plt.figure(figsize=(17, 9.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 0.8], hspace=0.6, wspace=0.42)
    ax7, ax4 = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[:, 1])
    g7, g4 = fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[1, 2])
    if linear:
        f = o[(o.outcome == outcome) & (o.status == "fitted")]
        lim = float(min(np.ceil(max(f.ci_low.abs().max(), f.ci_high.abs().max()) * 2) / 2, 3.0))
        forest_lin(ax7, o, 7, lim, outcome)
        forest_lin(ax4, o, 4, lim, outcome)
    else:
        F.forest(ax7, o, 7, outcome)
        F.forest(ax4, o, 4, outcome)
    F.gains(g7, mb, 7, outcome)
    F.gains(g4, mb, 4, outcome)
    if linear:
        for g in (g7, g4):
            g.set_xlabel("OLS R2 gain")
    fig.suptitle(title(vt, outcome, what), x=0.01, ha="left", fontweight="bold", fontsize=10.5)
    return save_figure(fig, STEP, name, {"forest_cl7": ax7, "forest_cl4": ax4, "gains_cl7": g7, "gains_cl4": g4})


def main():
    paths = figure("complete", "SV01b_forest_completeness", "H12 complete (538-547)", linear=False)
    paths += figure("mobility", "SV01b_forest_mobility", "Mobility (B z-score 528-547 vs core v3)", linear=True)
    for pth in paths:
        print(pth.relative_to(ROOT), f"{pth.stat().st_size / 1e3:.1f} kB")


if __name__ == "__main__":
    main()
