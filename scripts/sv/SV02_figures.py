"""SV02 figures, drawn only from the saved SV02 CSVs (rule 1.8).

1. SV02_helix_track (main): x = UniProt residue 305-550, y = fraction pydssp 'H' among the chains observing the
   position, lines for all chains, cl7 and cl4; dashed line = the 80 % rule; global consensus helices as a track;
   the three historical H12 windows (538-547, 537-546, 536-545) marked. Lower panel: fraction of chains observing each
   position (all, cl7, cl4), positions with < 50 % of all chains observing (low support) shaded.
   Source: results/sv/SV02/sse_per_position.csv, consensus_global.csv.
2. SV02_sse_heatmap (supplementary): pydssp state per chain (rows, sorted by cluster cl7, cl4, minor; matrix order
   inside) and position. Source: results/sv/SV02/sse_codes_per_chain.csv, results/sv/SV00/chain_table.csv.
3. rama/SV02_rama_<pos> (13 plots, positions fixed in SV02_phipsi.py): phi vs psi, one point per chain with both
   angles, colour = cluster. Source: results/sv/SV02/phipsi_per_chain.csv, rama_summary.csv (title context).
   Saved to figures/sv/SV02/rama/ (PNG + SVG) and each panel to figures/panels/sv/SV02/.
Run: .venv/bin/python -B -I scripts/sv/SV02_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.sv_fig import (AQUA, BLUE, CLUSTER_COLOR, GREY, INK, INK2, MAGENTA, MUTED, ORANGE,  # noqa: E402
                           mpl_setup, save_figure)

RES = ROOT / "results/sv/SV02"
STEP = "SV02"
WINDOWS = [("538-547 settings (primary)", 538, 547), ("537-546 helix table", 537, 546), ("536-545", 536, 545)]
SSE_COLORS = {"": "#ffffff", "-": "#e6e5e1", "E": MAGENTA, "H": INK2}


def helix_track(plt):
    pp = pd.read_csv(RES / "sse_per_position.csv")
    cg = pd.read_csv(RES / "consensus_global.csv")
    c7 = pd.read_csv(RES / "consensus_cl7.csv")
    c4 = pd.read_csv(RES / "consensus_cl4.csv")
    g = {k: pp[pp.group == k].set_index("position") for k in ("global", "cl7", "cl4")}

    def h12(c):
        x = c[c.helix_name == "H12"]
        return "; ".join(f"{int(r.start)}-{int(r.end)}" for r in x.itertuples()) if len(x) else "none"

    fig, (ax, bx) = plt.subplots(2, 1, figsize=(12, 7.0), sharex=True, gridspec_kw={"height_ratios": [3, 1.3],
                                                                                     "hspace": 0.08})
    style = {"global": (INK, 1.4, f"all chains (n={int(g['global'].n_chains.iloc[0])})"),
             "cl7": (BLUE, 1.2, f"cl7 (n={int(g['cl7'].n_chains.iloc[0])})"),
             "cl4": (ORANGE, 1.2, f"cl4 (n={int(g['cl4'].n_chains.iloc[0])})")}
    for k in ("cl7", "cl4", "global"):
        col, lw, lab = style[k]
        ax.plot(g[k].index, g[k].frac_alpha_H, color=col, lw=lw, label=lab, drawstyle="steps-mid",
                alpha=0.95 if k == "global" else 0.85)
    ax.axhline(0.8, color=INK2, ls="--", lw=0.9)
    ax.text(306, 0.81, "80 % rule", fontsize=7, color=INK2, va="bottom")
    # global consensus helices as a track
    for r in cg.itertuples():
        ax.add_patch(plt.Rectangle((r.start - 0.5, 1.035), r.end - r.start + 1, 0.04, color=INK2, lw=0, clip_on=False))
        ax.text((r.start + r.end) / 2, 1.085, r.helix_name, ha="center", va="bottom", fontsize=7, color=INK)
    ax.text(304, 1.055, "consensus\n(all chains)", ha="right", va="center", fontsize=6.5, color=INK2)
    # historical H12 windows
    for i, (lab, s, e) in enumerate(WINDOWS):
        y = 1.16 + 0.055 * i
        ax.plot([s - 0.5, e + 0.5], [y, y], color=MUTED, lw=2.2, solid_capstyle="butt", clip_on=False)
        ax.text(s - 2, y, lab, ha="right", va="center", fontsize=6.5, color=INK2)
        ax.axvspan(s - 0.5, e + 0.5, ymin=0, ymax=1, color=GREY, alpha=0.05, lw=0)
    ax.set_ylim(0, 1.0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("fraction 'H' among observing chains\n(pydssp: any helix)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.005, 1.0), ncol=1, fontsize=7.5)
    n = len(cg)
    ax.set_title(f"pydssp helix fraction per position: {n} global consensus helices (>= 4 positions with H in "
                 f">= 80 % of observing chains)\nconsensus H12 {h12(cg)} (all chains), {h12(c7)} (cl7), {h12(c4)} (cl4)",
                 fontsize=9, pad=88)
    # lower panel: observation
    low = g["global"][g["global"].low_support]
    for p in low.index:
        bx.axvspan(p - 0.5, p + 0.5, color=GREY, alpha=0.35, lw=0)
    for k in ("cl7", "cl4", "global"):
        col, lw, lab = style[k]
        bx.plot(g[k].index, g[k].frac_observing, color=col, lw=lw - 0.2, drawstyle="steps-mid")
    bx.axhline(0.5, color=INK2, ls="--", lw=0.8)
    bx.set_ylim(0, 1.02)
    bx.set_yticks([0, 0.5, 1.0])
    bx.set_ylabel("fraction of chains\nobserving")
    bx.set_xlabel("UniProt residue (P03372)")
    bx.set_xlim(304.5, 550.5)
    bx.set_xticks(np.arange(310, 551, 20))
    bx.text(345, 0.08, f"grey = low support (< 50 % of all chains observe it): {len(low)} positions",
            ha="left", va="bottom", fontsize=7, color=INK2)
    return save_figure(fig, STEP, "SV02_helix_track", panels={"helix_fraction": ax, "observing": bx})


def heatmap(plt):
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    m = pd.read_csv(RES / "sse_codes_per_chain.csv", index_col=0, keep_default_na=False, na_values=[""])
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv").set_index("chain_key")
    cl = ct.s08_cluster.reindex(m.index)
    rank = cl.map(lambda c: {7: 0, 4: 1}.get(int(c), 2 + int(c)))
    order = pd.DataFrame({"r": rank, "i": np.arange(len(m))}, index=m.index).sort_values(["r", "i"]).index
    m = m.loc[order]
    cl = cl.loc[order]
    codes = {"": 0, "-": 1, "E": 2, "H": 3}
    Z = m.fillna("").apply(lambda col: col.map(codes)).to_numpy(int)
    pos = np.array([int(c) for c in m.columns])
    fig = plt.figure(figsize=(12, 9))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.25, 12], wspace=0.01)
    sx = fig.add_subplot(gs[0])
    ax = fig.add_subplot(gs[1], sharey=sx)
    cmap = ListedColormap([SSE_COLORS[k] for k in ("", "-", "E", "H")])
    ax.imshow(Z, aspect="auto", interpolation="nearest", cmap=cmap, vmin=-0.5, vmax=3.5,
              extent=(pos[0] - 0.5, pos[-1] + 0.5, len(m) - 0.5, -0.5))
    cols = np.array([[int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)] for c in
                     [CLUSTER_COLOR["cl7"] if c == 7 else CLUSTER_COLOR["cl4"] if c == 4 else CLUSTER_COLOR["minor"]
                      for c in cl]]) / 255.0
    sx.imshow(cols[:, None, :], aspect="auto", interpolation="nearest", extent=(0, 1, len(m) - 0.5, -0.5))
    sx.set_xticks([])
    sx.set_ylabel("chains (sorted: cl7, cl4, minor clusters)")
    sx.set_yticks([])
    for spine in sx.spines.values():
        spine.set_visible(False)
    bounds = np.cumsum([int((cl == 7).sum()), int((cl == 4).sum())])
    for b in bounds:
        ax.axhline(b - 0.5, color=INK, lw=0.8)
    mids = {"cl7": bounds[0] / 2, "cl4": (bounds[0] + bounds[1]) / 2, "minor": (bounds[1] + len(m)) / 2}
    for k, y in mids.items():
        n = {"cl7": int((cl == 7).sum()), "cl4": int((cl == 4).sum()), "minor": int((~cl.isin([7, 4])).sum())}[k]
        sx.text(-0.3, y, f"{k}\n{n}", ha="right", va="center", fontsize=7, color=INK2, transform=sx.transData)
    ax.set_yticks([])
    ax.set_xlabel("UniProt residue (P03372)")
    ax.set_xticks(np.arange(310, 551, 20))
    ax.spines["left"].set_visible(False)
    h_frac = float((Z == 3).sum() / (Z > 0).sum())
    ax.set_title(f"pydssp state per chain and position: {len(m)} chains, {int((Z > 0).sum())} observed residues, "
                 f"{h_frac:.0%} of them H", fontsize=9)
    ax.legend(handles=[Patch(color=SSE_COLORS["H"], label="H (any helix)"), Patch(color=SSE_COLORS["E"], label="E"),
                       Patch(color=SSE_COLORS["-"], label="other ('-')"),
                       Patch(facecolor="white", edgecolor=MUTED, label="not observed")],
              loc="upper left", bbox_to_anchor=(0.0, -0.06), ncol=4, fontsize=7.5)
    return save_figure(fig, STEP, "SV02_sse_heatmap", panels={"heatmap": [sx, ax]})


def save_rama(fig, ax, name):
    from matplotlib.transforms import Bbox
    fdir = ROOT / "figures/sv" / STEP / "rama"
    pdir = ROOT / "figures/panels/sv" / STEP
    fdir.mkdir(parents=True, exist_ok=True)
    pdir.mkdir(parents=True, exist_ok=True)
    out = [fdir / f"{name}.png", fdir / f"{name}.svg", pdir / f"{name}__rama.png"]
    fig.savefig(out[0], dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(out[1], bbox_inches="tight", facecolor="white")
    fig.canvas.draw()
    bb = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    bb = Bbox.from_extents(bb.x0 - 0.05, bb.y0 - 0.05, bb.x1 + 0.05, bb.y1 + 0.05)
    fig.savefig(out[2], dpi=200, bbox_inches=bb, facecolor="white")
    return out


def rama(plt):
    df = pd.read_csv(RES / "phipsi_per_chain.csv")
    rs = pd.read_csv(RES / "rama_summary.csv")
    out = []
    for r in rs.itertuples():
        d = df[(df.position == r.position) & df.phi.notna() & df.psi.notna()]
        fig, ax = plt.subplots(figsize=(4.6, 4.9))
        for k, col, mk in (("cl7", BLUE, "o"), ("cl4", ORANGE, "o"), ("minor", AQUA, "^")):
            x = d[d.cluster == k]
            ax.scatter(x.phi, x.psi, s=9 if mk == "o" else 14, color=col, alpha=0.4, lw=0, marker=mk,
                       label=f"{k} (n={len(x)} chains, {x.entry_id.nunique()} entries)")
        ax.text(-60, -45, "α", fontsize=13, color=MUTED, alpha=0.45, ha="center", va="center")
        ax.text(-120, 130, "β", fontsize=13, color=MUTED, alpha=0.45, ha="center", va="center")
        ax.axhline(0, color=GREY, lw=0.5, alpha=0.5)
        ax.axvline(0, color=GREY, lw=0.5, alpha=0.5)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_yticks(np.arange(-180, 181, 90))
        ax.set_xlabel("φ (deg)")
        ax.set_ylabel("ψ (deg)")
        ax.set_title(f"Residue {r.position}: {r.context}", fontsize=9)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=7, handletextpad=0.2, markerscale=1.6)
        out += save_rama(fig, ax, f"SV02_rama_{r.position}")
        plt.close(fig)
    return out


def main():
    plt = mpl_setup()
    out = helix_track(plt) + heatmap(plt) + rama(plt)
    for p in out:
        print(p.relative_to(ROOT), f"{p.stat().st_size / 1e3:.1f} kB")


if __name__ == "__main__":
    main()
