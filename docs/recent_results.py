"""Slides: results added since the October report, in the order of the meeting requests (28 Oct, 4 Nov, 11 Nov 2026).

    python docs/recent_results.py        # writes docs/recent_results.pdf (16:9, vector, one finding per slide)

Every number on a slide is computed here from a saved table; each slide function names its source files in its
docstring. The charts are redrawn from the same tables as the figures in figures/ (same encodings, type sized for a
slide). Colours: group 7 (cl7, lid over the pocket) blue, group 4 (cl4, lid in the coactivator groove) orange, the
14 chains of the five small groups grey. 'Ordinary variation' = the larger of the within-group SD and the 95th
percentile of differences between copies of one crystal (common.ordinary_variation in pipeline/).
"""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.patches import Ellipse, FancyBboxPatch, RegularPolygon  # noqa: E402
from matplotlib.transforms import Affine2D  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "recent_results.pdf"
TITLE = "ESR1 ligand-binding domain: structural diversity of the crystal structures"
DATE = "13 November 2026"

G7, G4, GREY = "#2a78d6", "#eb6834", "#a3a29c"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#f1f0ec"
VIOLET = "#4a3aa7"
COLOR = {"cl7": G7, "cl4": G4, "minor": GREY}
MARK = {"cl7": "o", "cl4": "o", "minor": "^"}
NAME = {"cl7": "group 7", "cl4": "group 4", "minor": "grey chains"}

W, H = 13.333, 7.5                 # 16:9 page, inches
X0, X1 = 0.55, W - 0.55            # left / right margin
VB, VT = 1.18, H - 1.15            # visual box: bottom / top (inches from the page bottom)
VW, VH = X1 - X0, VT - VB

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 14, "axes.labelsize": 14.5, "axes.labelcolor": INK,
    "xtick.labelsize": 13.5, "ytick.labelsize": 13.5, "xtick.color": INK2, "ytick.color": INK2,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.9, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white", "legend.frameon": False,
    "lines.linewidth": 2.5, "pdf.fonttype": 42, "axes.titlesize": 14.5, "axes.titlelocation": "left",
    "axes.titlecolor": INK2, "axes.titleweight": "normal",
})


# ---------------------------------------------------------------------------------------------------------- helpers
def csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel, keep_default_na=False, na_values=[""])


def numbers(rel: str) -> dict:
    return json.loads((ROOT / rel).read_text())


def _width(fig, t) -> float:
    return t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi


def fit(fig, x, y, text, width, size, max_lines=1, min_size=12.0, **kw):
    """Text at (x, y) inches, wrapped to `width` inches in at most `max_lines` lines; the type shrinks if needed."""
    t = fig.text(x / W, y / H, "", fontsize=size, **kw)
    while size >= min_size:
        t.set_fontsize(size)
        for chars in range(max(len(text), 10), 9, -1):
            lines = textwrap.wrap(text, chars)
            t.set_text("\n".join(lines))
            if _width(fig, t) <= width:
                break
        if len(lines) <= max_lines and _width(fig, t) <= width:
            return t
        size -= 0.5
    raise ValueError(f"text does not fit: {text}")


def ax_at(fig, x, y, w, h, **kw):
    """Axes at (x, y, w, h) inches inside the visual box (origin at its lower left)."""
    return fig.add_axes([(X0 + x) / W, (VB + y) / H, w / W, h / H], **kw)


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.9)
    ax.set_axisbelow(True)


def notes(fig, *rows):
    """At most three short lines under the visual: (label, text)."""
    for i, (k, v) in enumerate(rows):
        y = 0.92 - 0.31 * i
        fig.text(X0 / W, y / H, k, fontsize=13.5, fontweight="bold", color=MUTED, va="center")
        fit(fig, X0 + 1.12, y, v, VW - 1.42, 14.5, va="center", color=INK)


def signed(v: float, fmt: str = "+.2f") -> str:
    """Number with a typographic minus sign."""
    return format(v, fmt).replace("-", "−")


def p_text(p: float) -> str:
    return "p < 0.001" if p < 0.001 else f"p = {p:.2f}"


def tiered_labels(ax, pts, base_y, min_gap, first, step, fontsize=12.5):
    """Label points (x, y, text) above / below base_y in tiers without overlaps (min_gap in x data units)."""
    last = {}
    for n, (x, y, txt) in enumerate(sorted(pts)):
        done = False
        for tier in range(6):
            for side in ((1, -1) if n % 2 == 0 else (-1, 1)):
                if (side, tier) not in last or x - last[(side, tier)] >= min_gap:
                    last[(side, tier)] = x
                    ax.annotate(txt, (x, y), (x, base_y + side * (first + tier * step)), fontsize=fontsize,
                                color=INK, ha="center", va="center",
                                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8, shrinkA=2, shrinkB=4))
                    done = True
                    break
            if done:
                break


class Deck:
    def __init__(self, path: Path):
        self.pdf = PdfPages(path, metadata={"Title": f"{TITLE}: results since the October report",
                                            "Subject": "Results for the requests of 28 October, 4 November and "
                                                       "11 November 2026", "CreationDate": None})
        self.n = 0

    def new(self, headline: str | None = None):
        self.n += 1
        fig = plt.figure(figsize=(W, H))
        if headline:
            fit(fig, X0, H - 0.30, headline, VW, 23, max_lines=2, min_size=20, va="top", fontweight="bold",
                color=INK, linespacing=1.15)
        if self.n > 1:
            fig.text((W - 0.32) / W, 0.30 / H, str(self.n), ha="right", va="center", fontsize=11, color=MUTED)
        return fig

    def add(self, fig):
        self.pdf.savefig(fig)
        plt.close(fig)

    def close(self):
        self.pdf.close()


# ----------------------------------------------------------------------------------------------------------- slides
def title_slide(deck):
    """Sources: results/01_structures/numbers.json (chains and entries analysed)."""
    n = numbers("results/01_structures/numbers.json")
    fig = deck.new()
    fit(fig, X0, 5.0, TITLE, VW, 34, max_lines=2, va="bottom", fontweight="bold", color=INK, linespacing=1.15)
    fig.text(X0 / W, 4.45 / H, "Results since the October report: requests of 28 October, 4 November and "
             "11 November", fontsize=19, color=INK2, va="top")
    fig.text(X0 / W, 3.85 / H, DATE, fontsize=19, color=INK2, va="top")
    fig.text(X0 / W, 1.0 / H, f"{n['chains_analysed']} chains of oestrogen receptor α from {n['entries_analysed']} "
             "PDB entries", fontsize=15, color=MUTED)
    deck.add(fig)


def asked_slide(deck):
    """Sources: docs/supervisor_feedback.md (requests); data/ligands/pocket_filling_tests.csv (pocket filling,
    cl7 vs cl4, entry-level Mann-Whitney, BH q)."""
    pf = csv("data/ligands/pocket_filling_tests.csv")
    q_fill = pf[(pf.comparison == "cluster") & (pf.level == "entry")].q_bh.iloc[0]
    rows = [
        ("28 October", "Does lid completeness follow the crystal or the biology?", "4, 15"),
        ("", "Helix ranges from consensus secondary structure", "4"),
        ("", "Compare incomplete helices by their midpoints", "5"),
        ("", "Colour the outliers in the clustering", "8, 10"),
        ("4 November", "Action classes against structure", "6"),
        ("", "Compare affinities measured in different assays", "7"),
        ("", "Pocket filling relative to the E2 pocket", f"no group difference (q {q_fill:.2f})"),
        ("", "Does an antagonist destabilise H12?", "open: class follows lid (6)"),
        ("11 November", "Outliers in lid-to-H11 distance: group 4, group 7, grey", "8, 9"),
        ("", "Sammon map, self-organising map, UMAP of PCA", "10"),
        ("", "Hide H12 and cluster again", "11"),
        ("", "Phi/psi, helix centres of mass, no superposition", "11, 12, 13"),
        ("", "Do crystallisation additives stabilise the lid?", "14, 15"),
        ("", "Biological variability or a real effect?", "rule on slide 3"),
    ]
    fig = deck.new("What was asked at the last three meetings, and where each answer is")
    xd, xr, xa = X0, X0 + 2.0, X0 + 8.75
    y = VT - 0.15
    for lab, x in (("Meeting", xd), ("Request", xr), ("Slide", xa)):
        fig.text(x / W, y / H, lab, fontsize=13, fontweight="bold", color=MUTED, va="center")
    y -= 0.22
    for date, req, ans in rows:
        if date:
            y -= 0.12
            fig.add_artist(plt.Line2D([X0 / W, X1 / W], [(y + 0.05) / H] * 2, color=GRID, lw=1))
            y -= 0.12
            fig.text(xd / W, y / H, date, fontsize=14.5, fontweight="bold", color=INK, va="center")
        fit(fig, xr, y, req, xa - xr - 0.3, 14.5, va="center", color=INK)
        fit(fig, xa, y, ans, X1 - xa, 14.5, va="center", color=INK, fontweight="bold" if ans[0].isdigit() else "normal")
        y -= 0.335
    deck.add(fig)


def capsule(ax, cx, cy, length, width, angle, color, z=2):
    p = FancyBboxPatch((cx - length / 2, cy - width / 2), length, width,
                       boxstyle=f"round,pad=0,rounding_size={width / 2}", fc=color, ec="none", zorder=z)
    p.set_transform(Affine2D().rotate_deg_around(cx, cy, angle) + ax.transData)
    ax.add_patch(p)


def receptor_sketch(ax, state):
    """Schematic ligand-binding domain: lid (H12) over the pocket (cl7) or in the coactivator groove (cl4)."""
    body, helix, helix_dark = "#ecebe7", "#dcdbd6", "#bdbcb6"
    ax.add_patch(FancyBboxPatch((0.8, 0.4), 7.6, 4.8, boxstyle="round,pad=0,rounding_size=1.5", fc=body, ec="none",
                                zorder=0))
    for cx, cy, length, a in ((3.7, 1.15, 2.6, 6), (6.0, 1.05, 2.0, -4), (3.9, 4.55, 1.6, 8)):
        capsule(ax, cx, cy, length, 0.55, a, helix, z=1)
    capsule(ax, 1.95, 2.8, 3.4, 0.62, 92, helix_dark)
    capsule(ax, 7.45, 2.7, 3.4, 0.62, 88, helix_dark)
    ax.text(1.95, 2.8, "H3", rotation=92, ha="center", va="center", fontsize=12.5, color=INK, zorder=3)
    ax.text(7.45, 2.7, "H11", rotation=88, ha="center", va="center", fontsize=12.5, color=INK, zorder=3)
    ax.add_patch(Ellipse((4.95, 2.85), 2.7, 1.65, fc="white", ec=MUTED, lw=1.2, ls=(0, (3, 2)), zorder=1))
    capsule(ax, 5.75, 4.2, 2.3, 0.75, 64, "white", z=1)        # opening of the pocket towards the lid site
    r = 0.30
    for x in (4.55, 4.55 + r * np.sqrt(3)):
        ax.add_patch(RegularPolygon((x, 2.75), 6, radius=r, fc=INK2, ec="white", lw=1, zorder=3))
    ax.text(4.95, 2.2, "ligand", ha="center", va="top", fontsize=12.5, color=INK, zorder=3)
    if state == "cl7":
        capsule(ax, 6.55, 5.3, 2.5, 0.64, -12, G7, z=4)
        ax.text(6.55, 5.95, "H12 (lid) closed over the pocket", ha="center", va="bottom", fontsize=12.5, color=INK)
        capsule(ax, 2.55, 5.25, 1.6, 0.42, 14, "#6f6e69", z=4)
        ax.text(2.25, 5.95, "coactivator peptide", ha="center", va="bottom", fontsize=12.5, color=INK)
    else:
        capsule(ax, 2.6, 5.3, 2.5, 0.64, 14, G4, z=4)
        ax.text(2.4, 5.95, "H12 (lid) in the coactivator groove", ha="center", va="bottom", fontsize=12.5, color=INK)
        ax.plot([5.3, 5.75, 5.95, 6.45], [2.85, 3.7, 4.55, 5.35], color=INK2, lw=2.2, zorder=3,
                solid_capstyle="round")
        ax.text(6.65, 5.45, "ligand side chain,\npocket left open", ha="left", va="bottom", fontsize=12.5, color=INK)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7)
    ax.set_aspect("equal")
    ax.axis("off")


def groups_slide(deck):
    """Sources: results/01_structures/numbers.json (chains, entries); results/02_lid_groups/numbers.json (group sizes,
    peptide-bound chains, number of groups at 1.5 A); results/03_lid_outliers/lid_distance_by_group.csv (median
    lid-to-H11 distance and copy 95th percentile per group)."""
    n1 = numbers("results/01_structures/numbers.json")
    n2 = numbers("results/02_lid_groups/numbers.json")
    lid = csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group")
    sizes = n2["group_sizes_1p5"]
    n_small = n2["n_groups_by_cut"]["1.5"] - 2
    fig = deck.new("Two lid positions split the chains into the two groups used on every slide: "
                   "group 7 and group 4")
    for i, (g, ttl) in enumerate((("cl7", "Group 7: lid closed over the pocket"),
                                  ("cl4", "Group 4: lid in the coactivator groove"))):
        x = 0.15 + i * 6.25
        receptor_sketch(ax_at(fig, x, 1.15, 5.7, 3.99), g)
        fig.text((X0 + x + 0.25) / W, (VB + 0.80) / H, ttl, fontsize=17, fontweight="bold", color=INK, va="center")
        fig.text((X0 + x + 0.25) / W, (VB + 0.35) / H,
                 f"{sizes[g]} chains · lid {lid.loc[g, 'median']:.1f} Å from H11 · peptide bound in "
                 f"{n2[f'peptide_bound_{g}']}", fontsize=14, color=INK2, va="center")
    notes(fig,
          ("Groups", f"{n1['chains_analysed']} chains from {n1['entries_analysed']} structures, clustered by C-alpha "
                     f"RMSD at 1.5 Å; {sizes['minor']} grey chains form {n_small} small groups"),
          ("Rule", "a difference counts only if it is larger than the spread within a group and between copies of "
                   "one crystal"),
          ("Example", f"lid-to-H11 distance: copies of one crystal differ by up to "
                      f"{lid.loc['cl7', 'copy_p95']:.2f} Å (group 7) and {lid.loc['cl4', 'copy_p95']:.2f} Å (group 4)"))
    deck.add(fig)


def h12_window_slide(deck):
    """Source: data/lid/h12_window_status.csv (chain level: lid status on window 538-547 and on 538-544)."""
    st = csv("data/lid/h12_window_status.csv")
    st = st[st.level == "chain"]
    rows = []
    for g in ("cl4", "cl7", "minor"):
        d = st[st.cluster_group == g]
        rows.append((g, int(d[d.h12_status_538_547 != "complete"].n.sum()),
                     int(d[d.h12_status_538_544 != "complete"].n.sum()), int(d.n.sum())))
    d4 = st[st.cluster_group == "cl4"]
    ends_only = int(d4[(d4.h12_status_538_547 == "partial") & (d4.h12_status_538_544 == "complete")].n.sum())
    fig = deck.new(f"Ending the lid window at residue 544 instead of 547 makes an incomplete lid rare: "
                   f"in group 4 it drops from {rows[0][1]} to {rows[0][2]} chains")
    ax = ax_at(fig, 1.9, 0.75, 9.9, 4.15)
    for y, (g, a, b, n) in zip((2, 1, 0), rows):
        pa, pb = 100 * a / n, 100 * b / n
        ax.annotate("", (pb, y), (pa, y), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6, mutation_scale=20,
                                                          shrinkA=9, shrinkB=9), zorder=1)
        ax.scatter(pa, y, s=190, facecolor="white", edgecolor=COLOR[g], lw=2.6, marker=MARK[g], zorder=3)
        ax.scatter(pb, y, s=190, color=COLOR[g], marker=MARK[g], zorder=3)
        if abs(pa - pb) > 8:
            ax.text(pa, y + 0.2, f"{a} of {n}", ha="center", va="bottom", fontsize=14, color=INK)
            ax.text(pb, y + 0.2, f"{b}", ha="center", va="bottom", fontsize=14, color=INK)
        else:
            ax.text(max(pa, pb) + 2.5, y, f"{a} of {n} → {b}", ha="left", va="center", fontsize=14, color=INK)
        if y == 2:
            ax.text(pa, y - 0.22, "window ends at 547", ha="center", va="top", fontsize=13, color=INK2)
            ax.text(pb, y - 0.22, "window ends at 544", ha="center", va="top", fontsize=13, color=INK2)
    ax.set_yticks([2, 1, 0], [f"{NAME[g]}\n{n} chains" for g, _, _, n in rows])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.5, 2.6)
    ax.set_xlim(-2, 100)
    ax.set_xlabel("Chains with an incomplete lid (% of the group)")
    grid(ax, "x")
    notes(fig,
          ("Shows", "chains whose lid (H12) is partly or not modelled, for the window 538–547 and for 538–544"),
          ("Means", f"{ends_only} of the {rows[0][1]} incomplete group 4 lids lack only residues 545–547: "
                    "the lid is there, its last turn is not"),
          ("Caveat", "the shorter window comes from consensus secondary structure and is not yet agreed"))
    deck.add(fig)


def midpoint_slide(deck):
    """Sources: data/lid/midpoint_error_by_k.csv (median midpoint shift of complete helices truncated by k residues);
    data/lid/midpoint_clustering_ari.csv (midpoint RMS distance, average linkage, K = 2, chain level, ARI vs groups)."""
    v = csv("data/lid/midpoint_error_by_k.csv")
    ari = csv("data/lid/midpoint_clustering_ari.csv")
    a = ari[(ari.distance == "midpoint_rms") & (ari.level == "chain") & (ari.K == 2)].ari.iloc[0]
    k4 = v[v.k_removed == 4].median_midpoint_error_A.max()
    fig = deck.new(f"A helix midpoint moves {k4:.2f} Å or less (median) when 4 end residues are missing, so helices "
                   "with missing ends can still be compared")
    ax = ax_at(fig, 0.9, 0.75, 9.6, 4.25)
    style = {"H12": (INK, 3.2), "H4/5": (VIOLET, 3.2)}
    ends = []
    for h in ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]:
        d = v[v.helix == h].set_index("k_removed").median_midpoint_error_A.dropna()
        c, lw = style.get(h, (MUTED, 1.6))
        ax.plot(d.index, d.values, color=c, lw=lw, marker="o", ms=6 if h in style else 4, zorder=3 if h in style else 2)
        ends.append([h, d.index[-1], d.values[-1]])
    for x in {e[1] for e in ends}:
        grp = sorted([e for e in ends if e[1] == x], key=lambda e: e[2])
        for k in range(1, len(grp)):
            grp[k][2] = max(grp[k][2], grp[k - 1][2] + 0.075)
        for h, _, y in grp:
            ax.text(x + 0.12, y, h, fontsize=13, va="center", color=INK if h in style else INK2,
                    fontweight="bold" if h in style else "normal")
    ax.axhline(1.0, color=INK, lw=1.2)
    ax.text(1.0, 1.03, "1 Å", fontsize=13, color=INK, va="bottom")
    ax.set_xlabel("Residues missing at one end of a complete helix")
    ax.set_ylabel("Median shift of the helix midpoint (Å)")
    ax.set_xlim(0.8, 6.7)
    ax.set_xticks(range(1, 7))
    grid(ax, "y")
    first_fail = int(v[(v.helix == "H4/5") & (v.median_midpoint_error_A > 1)].k_removed.min())
    notes(fig,
          ("Shows", "complete helices cut short by 1–6 residues; how far their midpoint moves"),
          ("Means", f"distances between helix midpoints reproduce the two lid groups (adjusted Rand index {a:.3f})"),
          ("Caveat", f"H4/5 crosses 1 Å at {first_fail} missing residues; the agreement mostly restates the H12 split"))
    deck.add(fig)


def class_slide(deck):
    """Source: data/ligands/chain_class.csv (uncurated ChEMBL mechanism class per chain; an entry counts in a group if
    at least one of its chains has that class and group)."""
    cc = csv("data/ligands/chain_class.csv")
    cc["action_class"] = cc.action_class.fillna("unlabelled")
    classes = ["agonist", "SERM", "SERD", "unlabelled"]
    groups = ["cl7", "cl4", "minor"]
    m = {(k, g): cc[(cc.action_class == k) & (cc.cluster == g)].entry_id.nunique() for k in classes for g in groups}
    tot = {k: cc[cc.action_class == k].entry_id.nunique() for k in classes}
    n_lab = cc[cc.action_class != "unlabelled"].entry_id.nunique()
    fig = deck.new(f"Action classes almost coincide with the lid groups: agonists {m[('agonist', 'cl7')]} of "
                   f"{tot['agonist']} entries in group 7, SERMs {m[('SERM', 'cl4')]} of {tot['SERM']} in group 4")
    ax = ax_at(fig, 2.6, 0.75, 8.6, 4.3)
    bh = 0.26
    for i, k in enumerate(classes):
        for j, g in enumerate(groups):
            y = i + (j - 1) * (bh + 0.04)
            share = m[(k, g)] / tot[k]
            ax.barh(y, 100 * share, height=bh, color=COLOR[g], edgecolor="white", linewidth=0)
            lab = f"{m[(k, g)]} of {tot[k]}" + (f"  {NAME[g]}" if i == 0 else "")
            ax.text(100 * share + 1.0, y, lab, va="center", fontsize=13, color=INK)
    ax.set_yticks(range(len(classes)), [f"{k}\n{tot[k]} entries" for k in classes])
    ax.set_ylim(len(classes) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 112)
    ax.set_xticks(range(0, 101, 25))
    ax.set_xlabel("Entries of the class with a chain in the group (%)")
    grid(ax, "x")
    notes(fig,
          ("Shows", "for each action class, the share of its crystal structures that hold a chain in each lid group"),
          ("Means", "class cannot yet be tested apart from lid position, nor its link to H12 stability or the outliers"),
          ("Caveat", f"uncurated ChEMBL labels: only {n_lab} of {cc.entry_id.nunique()} entries carry a class"))
    deck.add(fig)


def assay_slide(deck):
    """Sources: data/ligands/assay_spread.csv (SD across assays of absolute pChEMBL and of pX - pX(E2) in the same
    assay, ligands with >= 2 such assays); data/ligands/ligand_ids.csv (ligand names)."""
    a = csv("data/ligands/assay_spread.csv").sort_values("sd_px_abs", ascending=False).reset_index(drop=True)
    names = csv("data/ligands/ligand_ids.csv").drop_duplicates("ccd").set_index("ccd").parent_pref_name
    better = int((a.sd_dpx < a.sd_px_abs).sum())
    fig = deck.new(f"Referencing each affinity to E2 in the same assay narrows the spread between assays for {better} "
                   f"of {len(a)} ligands (median SD {a.sd_px_abs.median():.2f} → {a.sd_dpx.median():.2f})")
    ax = ax_at(fig, 3.7, 0.75, 8.3, 4.3)
    for y, r in a.iterrows():
        c = COLOR.get(r.ligand_group, MUTED)
        ax.annotate("", (r.sd_dpx, y), (r.sd_px_abs, y), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.5,
                                                                          mutation_scale=18, shrinkA=8, shrinkB=8))
        ax.scatter(r.sd_px_abs, y, s=170, facecolor="white", edgecolor=c, lw=2.4, zorder=3)
        ax.scatter(r.sd_dpx, y, s=170, color=c, zorder=3)
        ax.text(max(r.sd_px_abs, r.sd_dpx) + 0.1, y, f"{int(r.n_paired_assays)} assays", va="center", fontsize=13,
                color=INK2)
        if y == 0:
            ax.text(r.sd_px_abs, y - 0.3, "absolute", ha="center", va="bottom", fontsize=13, color=INK2)
            ax.text(r.sd_dpx, y - 0.3, "relative to E2", ha="center", va="bottom", fontsize=13, color=INK2)
    ax.set_yticks(range(len(a)), [f"{names[r.ccd].capitalize()}\n{r.action_class}, {NAME[r.ligand_group]}"
                                  for _, r in a.iterrows()])
    ax.set_ylim(len(a) - 0.5, -0.75)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, a[["sd_px_abs", "sd_dpx"]].max().max() + 0.75)
    ax.set_xlabel("Spread of the affinity across assays (SD, log units)")
    grid(ax, "x")
    notes(fig,
          ("Shows", "open: absolute affinities from different assays; filled: the same affinities minus E2 in each assay"),
          ("Means", "affinities from different systems become more comparable when each is referenced to E2"),
          ("Caveat", f"only {len(a)} ligands have two or more assays that also measured E2"))
    deck.add(fig)


def lid_outliers_slide(deck):
    """Sources: results/03_lid_outliers/lid_distance_per_chain.csv (lid-to-H11 distance per chain),
    results/03_lid_outliers/lid_distance_by_group.csv (medians, copy 95th percentiles),
    results/03_lid_outliers/lid_outliers.csv (examined chains, beyond_ordinary, lid_state, deviation_A)."""
    g = csv("results/03_lid_outliers/lid_distance_per_chain.csv")
    summ = csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group")
    fl = csv("results/03_lid_outliers/lid_outliers.csv")
    fl["beyond"] = fl.beyond_ordinary.astype(str).eq("True")
    col = "lid_to_h11_A"
    dev = fl[fl.group == "cl7"].set_index("chain_key").deviation_A
    top2 = sorted(dev.sort_values(ascending=False).index[:2])
    lab = {k: k.rsplit("|", 1)[0].replace("|", " ") for k in fl.chain_key}
    between = fl[(fl.group == "minor") & (fl.lid_state == "between cl7 and cl4")]
    n_ext4 = int(fl[(fl.group == "cl4") & fl.beyond].shape[0])
    assert int(summ.loc["cl4", "n_tukey"]) == 0     # the cl4 chains beyond copy variation are not Tukey outliers
    fig = deck.new(f"Two group 7 chains, {lab[top2[0]]} and {lab[top2[1]]}, sit furthest out; group 4 has only "
                   f"{n_ext4} mild extremes; {len(between)} grey chains hold the lid in between")
    ax = ax_at(fig, 1.55, 0.75, 10.6, 4.35)
    rows = {"cl7": 3.0, "minor": 1.5, "cl4": 0.0}
    rng = np.random.default_rng(42)
    for grp in ("cl7", "cl4"):
        m, thr = summ.loc[grp, "median"], summ.loc[grp, "ordinary_threshold"]
        ax.axvspan(m - thr, m + thr, color=COLOR[grp], alpha=0.10, lw=0)
        ax.text(m, -0.62, f"{NAME[grp]}: {m:.1f} ± {thr:.2f} Å", ha="center", va="top", fontsize=13, color=INK2)
    beyond = set(fl[fl.beyond].chain_key)
    for grp, y in rows.items():
        d = g[(g.group == grp) & g[col].notna()]
        jit = rng.uniform(-0.22, 0.22, len(d))
        hit = d.chain_key.isin(beyond).to_numpy()
        ax.scatter(d[col][~hit], y + jit[~hit], s=26 if grp != "minor" else 70, color=COLOR[grp], alpha=0.5, lw=0,
                   marker=MARK[grp])
        ax.scatter(d[col][hit], y + jit[hit], s=95, color=COLOR[grp], marker=MARK[grp], edgecolor=INK, lw=1.3,
                   zorder=3)
        pts = sorted(zip(d[col][hit], y + jit[hit], [lab[k] for k in d.chain_key[hit]]))
        if grp == "minor":
            tiered_labels(ax, pts, y, min_gap=1.05, first=0.55, step=0.25)
            continue
        runs = [[pts[0]]]                       # chains within 0.5 A of each other share one label
        for p in pts[1:]:
            if p[0] - runs[-1][-1][0] > 0.5:
                runs.append([p])
            else:
                runs[-1].append(p)
        for k, run in enumerate(runs):
            x = float(np.mean([p[0] for p in run]))
            side = 1 if k % 2 == 0 else -1
            bold = any(p[2] in (lab[top2[0]], lab[top2[1]]) for p in run)
            ax.annotate(", ".join(p[2] for p in run), (x, y + side * 0.22), (x, y + side * 0.55), ha="center",
                        va="center", fontsize=12.5, color=INK, fontweight="bold" if bold else "normal",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    ax.set_yticks(list(rows.values()), ["group 7", "grey chains", "group 4"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-1.0, 3.75)
    ax.set_xlim(10.8, 25.0)
    ax.set_xlabel("Distance from the lid (H12) centre to the H11 axis (Å)")
    grid(ax, "x")
    n_na = int(g[col].isna().sum())
    notes(fig,
          ("Shows", "every chain; shaded: group median ± the largest normal difference between copies of one crystal"),
          ("Means", f"{len(beyond)} chains lie beyond ordinary variation (ringed): these are the structures to "
                    "inspect one by one"),
          ("Caveat", f"{n_na} chains without a helical H12 have no distance"))
    deck.add(fig)


def copies_slide(deck):
    """Sources: results/03_lid_outliers/lid_distance_per_chain.csv, results/03_lid_outliers/lid_outliers.csv,
    results/03_lid_outliers/lid_distance_by_group.csv."""
    g = csv("results/03_lid_outliers/lid_distance_per_chain.csv")
    fl = csv("results/03_lid_outliers/lid_outliers.csv")
    summ = csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group")
    col = "lid_to_h11_A"
    order = {"cl7": 0, "minor": 1, "cl4": 2}
    ents = fl.assign(o=fl.group.map(order)).sort_values(["o", col]).entry_id.drop_duplicates().tolist()
    ents = [e for e in ents if (g.entry_id == e).sum() > 1]          # crystals with more than one copy
    flagged = set(fl.chain_key)
    n_alone = 0
    for e in ents:
        d = g[g.entry_id == e]
        others = d[~d.chain_key.isin(flagged)][col].dropna()
        inside = others.apply(lambda x: any(abs(x - summ.loc[k, "median"]) <= summ.loc[k, "ordinary_threshold"]
                                            for k in ("cl7", "cl4")))
        n_alone += int(len(others) > 0 and inside.all())
    sp = g[g.entry_id == "9B25"][col]
    fig = deck.new(f"An unusual lid usually belongs to one copy, not to the whole crystal: in {n_alone} of "
                   f"{len(ents)} crystals the other copies are ordinary")
    ax = ax_at(fig, 1.0, 0.75, 11.0, 4.4)
    for grp in ("cl7", "cl4"):
        m, thr = summ.loc[grp, "median"], summ.loc[grp, "ordinary_threshold"]
        ax.axvspan(m - thr, m + thr, color=COLOR[grp], alpha=0.10, lw=0)
    for y, e in enumerate(ents):
        d = g[g.entry_id == e].dropna(subset=[col]).sort_values(col)
        if d.empty:
            ax.text(12.0, y, "no copy has a helical H12", va="center", fontsize=12.5, color=INK2)
            continue
        if len(d) > 1:
            ax.plot([d[col].min(), d[col].max()], [y, y], color=GRID, lw=3, zorder=1, solid_capstyle="round")
        for _, r in d.iterrows():
            f = r.chain_key in flagged
            ax.scatter(r[col], y, s=95 if f else 60, color=COLOR[r.group] if f else "white",
                       edgecolor=INK if f else COLOR[r.group], lw=1.3 if f else 1.8, marker=MARK[r.group], zorder=3)
        lettered = d[d.chain_key.isin(flagged)]
        runs = []                               # copy letters; copies within 0.3 A share one label
        for _, r in lettered.iterrows():
            if runs and r[col] - runs[-1][-1][0] <= 0.3:
                runs[-1].append((r[col], r.chain_key.split("|")[1]))
            else:
                runs.append([(r[col], r.chain_key.split("|")[1])])
        for run in runs:
            x = run[-1][0]
            right_free = not ((d[col] > x) & (d[col] - x < 0.6)).any()
            ax.annotate(", ".join(c for _, c in run), (x, y), (10 if right_free else -10, 0),
                        textcoords="offset points", fontsize=11.5, color=INK2, va="center",
                        ha="left" if right_free else "right")
    ax.set_yticks(range(len(ents)), ents)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(len(ents) - 0.4, -0.8)
    ax.set_xlim(10.8, 25.0)
    ax.set_xlabel("Distance from the lid (H12) centre to the H11 axis (Å), every copy of each crystal")
    grid(ax, "x")
    ax.text(summ.loc["cl7", "median"], -0.75, "group 7 range", ha="center", va="bottom", fontsize=13, color=INK2)
    ax.text(summ.loc["cl4", "median"], -0.75, "group 4 range", ha="center", va="bottom", fontsize=13, color=INK2)
    notes(fig,
          ("Shows", "crystals holding an examined chain (filled, letter = copy) next to their other copies (open)"),
          ("Means", f"what copies share (ligand, drop) does not fix the lid: 9B25 spans {sp.min():.1f}–{sp.max():.1f} Å "
                    "in one crystal"),
          ("Caveat", "copies of one crystal sit in different lattice contacts (slide 15)"))
    deck.add(fig)


def maps_slide(deck):
    """Sources: results/05_maps/maps_coords.csv (Sammon, SOM, UMAP-on-PCA coordinates), results/05_maps/numbers.json
    (silhouettes, PCA components, SOM grid), results/05_maps/maps_quality.csv (Spearman of map vs RMSD distance),
    results/05_maps/maps_islands.csv (HDBSCAN islands: gap vs copy RMSD 95th percentile),
    results/03_lid_outliers/lid_outliers.csv (examined chains)."""
    m = csv("results/05_maps/maps_coords.csv")
    fl = csv("results/03_lid_outliers/lid_outliers.csv")
    n = numbers("results/05_maps/numbers.json")
    q = csv("results/05_maps/maps_quality.csv").set_index("map")
    isl = csv("results/05_maps/maps_islands.csv")
    rng = np.random.default_rng(42)
    m["som_x"] = m.som_col + rng.uniform(-0.32, 0.32, len(m))
    m["som_y"] = m.som_row + rng.uniform(-0.32, 0.32, len(m))
    sil = n["silhouette_cl7_cl4"]
    gap = isl.gap_A.max()
    copy95 = isl.copy_rmsd_p95.iloc[0]
    fig = deck.new(f"All three maps show the same two groups; no smaller island stands out more than copies of one "
                   f"crystal differ ({gap:.2f} vs {copy95:.2f} Å)")
    panels = [("sammon1", "sammon2", "Sammon map", sil["sammon"]),
              ("som_x", "som_y", f"Self-organising map ({n['som_grid'][0]} × {n['som_grid'][1]})", sil["som"]),
              ("umap1", "umap2", f"UMAP of {n['pca_components_95pct']} PCA scores", sil["umap_on_pca"])]
    flagged = m.chain_key.isin(fl.chain_key)
    halo = [patheffects.withStroke(linewidth=5, foreground="white")]
    for i, (x, y, ttl, s_) in enumerate(panels):
        ax = ax_at(fig, i * 4.14, 0.05, 3.95, 4.3)
        for grp in ("cl7", "cl4", "minor"):
            d = m[(m.group == grp) & ~flagged]
            ax.scatter(d[x], d[y], s=16 if grp != "minor" else 40, color=COLOR[grp], alpha=0.55, lw=0,
                       marker=MARK[grp])
        d = m[flagged]
        ax.scatter(d[x], d[y], s=55, facecolor=[COLOR[k] for k in d.group], edgecolor=INK, lw=1.1, marker="D",
                   zorder=3)
        ax.set_title(f"{ttl}\ngroup separation {s_:.2f} (silhouette)", fontsize=13.5, color=INK2, loc="left")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp_ in ("left", "bottom"):
            ax.spines[sp_].set_color(GRID)
        for grp in ("cl7", "cl4"):
            c = m[m.group == grp][[x, y]].median()
            ax.text(c[x], c[y], NAME[grp], ha="center", va="center", fontsize=15, fontweight="bold", color=INK,
                    zorder=5, path_effects=halo)
        ax.margins(0.08)
    notes(fig,
          ("Shows", "each chain placed by its RMSD to all others; diamonds: chains examined as lid outliers or grey"),
          ("Means", f"the Sammon map keeps RMSD distances best (rank correlation "
                    f"{q.loc['sammon', 'spearman_dist']:.2f}); the smaller islands follow crystal form"))
    deck.add(fig)


def lid_hidden_slide(deck):
    """Sources: results/06_lid_hidden/no_lid_k_scan.csv (average linkage on RMSD without residues 528-550, ARI vs
    groups), results/07_geometry_free/geometry_free_kscan.csv (Ward on superposition-free helix features with and
    without H12), results/06_lid_hidden/no_lid_separation.csv (gap between split groups vs copy 95th percentile)."""
    a = csv("results/06_lid_hidden/no_lid_k_scan.csv")
    b = csv("results/07_geometry_free/geometry_free_kscan.csv")
    sep = csv("results/06_lid_hidden/no_lid_separation.csv")
    with_h12 = b[b.features == "all_features"]
    no_h12 = b[b.features == "without_H12"]
    top = with_h12[with_h12.k == 2].ari_vs_rmsd_groups.iloc[0]
    best_without = max(no_h12.ari_vs_rmsd_groups.max(), a.ari_vs_lid_groups.max())
    fig = deck.new(f"Hide H12 and the groups disappear: helix geometry with H12 recovers them (agreement {top:.2f}), "
                   f"without H12 none exceeds {best_without:.2f}")
    ax = ax_at(fig, 1.0, 0.75, 7.3, 4.3)
    series = [(with_h12, "ari_vs_rmsd_groups", "helix geometry with H12", INK),
              (no_h12, "ari_vs_rmsd_groups", "helix geometry without H12", MUTED),
              (a, "ari_vs_lid_groups", "RMSD without residues 528–550", VIOLET)]
    ends = []
    for d, c_, lab, col in series:
        ax.plot(d.k, d[c_], color=col, lw=3, marker="o", ms=7)
        ends.append([d[c_].iloc[-1], lab])
    ends.sort()
    for k in range(1, len(ends)):
        ends[k][0] = max(ends[k][0], ends[k - 1][0] + 0.07)
    for y, lab in ends:
        ax.text(10.35, y, lab, fontsize=13.5, va="center", color=INK, clip_on=False)
    ax.set_xlabel("Number of clusters cut")
    ax.set_ylabel("Agreement with the lid groups\n(adjusted Rand index; 1 = identical)")
    ax.set_ylim(-0.03, 1.0)
    ax.set_xticks(range(2, 11))
    ax.set_xlim(1.7, 10.2)
    grid(ax, "y")
    best = sep.sort_values("gap").iloc[-1]
    k10 = a[a.k == a.k.max()].ari_vs_lid_groups.iloc[0]
    notes(fig,
          ("Shows", "how well each clustering finds group 7 and group 4, for 2 to 10 clusters"),
          ("Means", f"without the lid the chains form one ensemble: no split exceeds copy-to-copy differences "
                    f"({best.gap:.2f} vs {best.copy_p95:.2f} Å)"),
          ("Caveat", f"the {int(a.k.max())}-cluster RMSD split partly follows the lid groups ({k10:.2f}) but is smaller "
                     "than copy differences"))
    deck.add(fig)


def helix_pairs_slide(deck):
    """Source: results/07_geometry_free/geometry_free_features.csv (cl4 minus cl7 median of each helix-helix and
    helix-chain centre-of-mass distance; real = beyond max(pooled SD, copy 95th percentile))."""
    ft = csv("results/07_geometry_free/geometry_free_features.csv")
    ft = ft[ft.kind == "distance"].copy()
    ft["real"] = ft.real.astype(str).eq("True")
    ft["h12"] = (ft.part_1 == "H12") | (ft.part_2 == "H12")
    other = ft[~ft.h12]
    fig = deck.new(f"Without superposition only distances to H12 change: up to {ft.diff_cl4_minus_cl7.abs().max():.1f} "
                   f"Å, while all {len(other)} other helix pairs move {other.diff_cl4_minus_cl7.abs().max():.2f} Å "
                   "or less")
    ax = ax_at(fig, 2.6, 0.75, 9.4, 4.3)
    ax.axvline(0, color=MUTED, lw=1)
    h = ft[ft.h12].sort_values("diff_cl4_minus_cl7")
    partner = [r.part_1 if r.part_2 == "H12" else r.part_2 for r in h.itertuples()]
    partner = ["chain centre" if p == "COM" else p for p in partner]
    ax.scatter(h.diff_cl4_minus_cl7, np.ones(len(h)), s=150, color=np.where(h.real, INK, "white"), edgecolor=INK,
               lw=1.6, zorder=3)
    for k, (p, v) in enumerate(zip(partner, h.diff_cl4_minus_cl7)):
        side = 1 if k % 2 == 0 else -1
        ax.text(v, 1 + side * 0.17, f"{p}\n{signed(v, '+.1f')}", ha="center", va="bottom" if side > 0 else "top",
                fontsize=13, color=INK)
    rng = np.random.default_rng(42)
    ax.scatter(other.diff_cl4_minus_cl7, rng.uniform(-0.1, 0.1, len(other)), s=80,
               color=np.where(other.real, INK, "white"), edgecolor=MUTED, lw=1.4, zorder=3)
    ax.text(other.diff_cl4_minus_cl7.max() + 0.35, 0, f"all within ±{other.diff_cl4_minus_cl7.abs().max():.2f} Å",
            va="center", fontsize=13, color=INK2)
    ax.set_yticks([1, 0], [f"distances to H12\n({len(h)})", f"all other\nhelix pairs ({len(other)})"])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(-0.55, 1.6)
    lim = np.ceil(ft.diff_cl4_minus_cl7.abs().max() + 0.6)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("Group 4 minus group 7 (Å, medians): > 0 farther apart in group 4")
    grid(ax, "x")
    notes(fig,
          ("Shows", "distances between helix centres of mass from the deposited coordinates; filled = beyond ordinary "
                    "variation"),
          ("Means", f"{int(ft.real.sum())} of {len(ft)} distances change, all involving H12: the rest of the domain "
                    "keeps its shape"),
          ("Caveat", "a helix centre needs at least 70% of the helix modelled"))
    deck.add(fig)


def h12_mobility_slide(deck):
    """Source: results/07_geometry_free/geometry_free_per_chain.csv (len_H12 from phi/psi, mobility_h12_z = H12
    C-alpha B z-score vs core); Spearman on entry means within each group."""
    per = csv("results/07_geometry_free/geometry_free_per_chain.csv")
    res = {}
    for grp in ("cl7", "cl4"):
        e = per[per.group == grp].groupby("entry_id")[["len_H12", "mobility_h12_z"]].mean().dropna()
        r, p = spearmanr(e.len_H12, e.mobility_h12_z)
        res[grp] = (e, r, p)
    fig = deck.new(f"In group 4 a longer H12 is a less mobile H12 (ρ {signed(res['cl4'][1])}, {len(res['cl4'][0])} "
                   f"crystals); in group 7 there is no link (ρ {signed(res['cl7'][1])})")
    rng = np.random.default_rng(42)
    lo = min(e.mobility_h12_z.min() for e, _, _ in res.values())
    hi = max(e.mobility_h12_z.max() for e, _, _ in res.values())
    for i, grp in enumerate(("cl7", "cl4")):
        e, r, p = res[grp]
        ax = ax_at(fig, 1.0 + i * 5.95, 0.75, 5.0, 3.95)
        ax.scatter(e.len_H12 + rng.uniform(-0.15, 0.15, len(e)), e.mobility_h12_z, s=34, color=COLOR[grp],
                   alpha=0.6, lw=0)
        gb = e.groupby(e.len_H12.round()).mobility_h12_z
        med = gb.median()[gb.size() >= 5]       # median only where a length holds at least 5 crystals
        ax.plot(med.index, med.values, color=INK, lw=2.4)
        ax.set_title(f"{NAME[grp]}: ρ {signed(r)} ({p_text(p)}), {len(e)} crystals", fontsize=14, color=INK)
        ax.set_xlabel("H12 length from phi/psi (residues)")
        ax.set_ylim(lo - 0.3, hi + 0.3)
        if i == 0:
            ax.set_ylabel("H12 mobility\n(B-factor z-score vs core)")
        grid(ax, "y")
        if grp == "cl4":
            x0 = med.index[-1]
            ax.annotate("median per length\n(≥ 5 crystals)", (x0, med.loc[x0]), (x0, med.loc[x0] + 1.1), fontsize=13,
                        color=INK, ha="center", va="bottom", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8))
    notes(fig,
          ("Shows", "one point per crystal: helical length of H12 against how mobile it is"),
          ("Means", "in the groove, a partly unwound lid is also a looser lid; in group 7 H12 is almost always complete"),
          ("Caveat", "B-factors also depend on refinement (TLS) and resolution"))
    deck.add(fig)


def crystal_slide(deck):
    """Sources: results/08_crystal_environment/additives_lattice_tests.csv (entry-level effect, ordinary-variation
    threshold and BH q per predictor, outcome and group), results/08_crystal_environment/numbers.json (number of tests,
    real effects, median additives per chain, entries with an additive)."""
    t = csv("results/08_crystal_environment/additives_lattice_tests.csv").dropna(subset=["effect"])
    t["ratio"] = t.effect / t.threshold
    n = numbers("results/08_crystal_environment/numbers.json")
    labels = {"additives_per_chain": "bound additives", "other_het_per_chain": "other small molecules",
              "recipe_classes": "recipe components", "solvent_pct": "solvent content",
              "sym_atoms": "lattice contacts at lid", "n_residues_sym_contact": "lattice contacts, chain"}
    order = [p for p in labels if p in set(t.predictor)]       # predictors with an effect in at least one test
    fig = deck.new(f"Additives, recipe, solvent content and lattice contacts never move lid order or mobility beyond "
                   f"ordinary variation ({n['n_real']} of {n['n_tests']} tests)")
    for i, (o, ttl) in enumerate((("lid_order", "lid order (residues 538–547 modelled)"),
                                  ("mobility_h12_z", "lid mobility (B-factor z-score)"))):
        ax = ax_at(fig, 2.95 + i * 4.75, 0.75, 4.3, 3.95)
        ax.axvspan(-1, 1, color=BAND, lw=0)
        ax.axvline(0, color=MUTED, lw=1)
        for grp, off in (("cl7", -0.14), ("cl4", 0.14)):
            d = t[(t.outcome == o) & (t.group == grp)].set_index("predictor").reindex(order).dropna(subset=["ratio"])
            y = np.array([order.index(p) for p in d.index]) + off
            sig = (d.q < 0.05).to_numpy()
            ax.scatter(d.ratio[sig], y[sig], s=110, color=COLOR[grp], zorder=3)
            ax.scatter(d.ratio[~sig], y[~sig], s=110, facecolor="white", edgecolor=COLOR[grp], lw=2, zorder=3)
            if i == 0:
                ax.text(d.ratio.iloc[0] + 0.12, order.index(d.index[0]) + off, NAME[grp], va="center", ha="left",
                        fontsize=12.5, color=INK)
        ax.set_yticks(range(len(order)), [labels[p] for p in order] if i == 0 else [])
        ax.set_ylim(len(order) - 0.5, -0.5)
        ax.tick_params(axis="y", length=0)
        ax.spines["left"].set_visible(False)
        ax.set_xlim(-1.3, 1.3)
        ax.set_title(ttl, fontsize=14, color=INK)
        ax.set_xlabel("Effect ÷ ordinary variation")
        grid(ax, "x")
    mx = t.loc[t.ratio.abs().idxmax()]
    notes(fig,
          ("Shows", "shaded ±1 = ordinary variation; filled = detectable (q < 0.05) yet smaller than ordinary variation"),
          ("Means", f"more components in the drop do not make the lid more ordered or less mobile (largest "
                    f"{abs(mx.ratio):.2f} of ordinary)"),
          ("Caveat", f"bound additives are rare: median {n['median_additives_per_chain']:.0f} per chain, "
                     f"{n['entries_with_additive']} of {n['n_entries']} entries have any"))
    deck.add(fig)


def lattice_slide(deck):
    """Source: results/08_crystal_environment/lattice_within_crystal.csv (contact = symmetry-mate atoms on 528-550;
    per crystal, does the copy with more contacts have the better lid; sign test, BH q)."""
    w = csv("results/08_crystal_environment/lattice_within_crystal.csv")
    w = w[w.contact == "sym_atoms"]
    rows = [("lid_order", "cl4", "lid order"), ("lid_order", "cl7", "lid order"),
            ("mobility_h12_z", "cl4", "lid mobility"), ("mobility_h12_z", "cl7", "lid mobility")]
    r4 = w[(w.outcome == "lid_order") & (w.group == "cl4")].iloc[0]
    fig = deck.new(f"Within one crystal, the group 4 copy with more lattice contacts at the lid has more of the lid "
                   f"modelled: {int(r4.more_contacts_better)} crystals against {int(r4.more_contacts_worse)}")
    ax = ax_at(fig, 2.1, 0.75, 6.6, 4.3)
    for y, (o, g, lab) in enumerate(rows):
        r = w[(w.outcome == o) & (w.group == g)].iloc[0]
        tot = r.more_contacts_better + r.more_contacts_worse + r.ties
        left = 0
        for v, c, a in ((r.more_contacts_better, COLOR[g], 1.0), (r.ties, BAND, 1.0),
                        (r.more_contacts_worse, COLOR[g], 0.35)):
            ax.barh(y, v / tot, left=left, height=0.58, color=c, alpha=a, edgecolor="white", linewidth=2)
            left += v / tot
        ax.text(1.04, y, f"{int(r.more_contacts_better)} · {int(r.ties)} · {int(r.more_contacts_worse)}", va="center",
                fontsize=14, color=INK, clip_on=False)
        ax.text(1.42, y, f"q {r.q:.2g}", va="center", fontsize=14, color=INK, clip_on=False)
    ax.text(1.04, -0.5, "better · equal · worse", va="bottom", fontsize=13, color=INK2, clip_on=False)
    tot = r4.more_contacts_better + r4.more_contacts_worse + r4.ties
    for x, txt in ((r4.more_contacts_better / tot / 2, "better"), ((r4.more_contacts_better + r4.ties / 2) / tot, "equal"),
                   (1 - r4.more_contacts_worse / tot / 2, "worse")):
        ax.text(x, -0.42, txt, ha="center", va="bottom", fontsize=13, color=INK2)
    ax.set_yticks(range(len(rows)), [f"{lab}, {NAME[g]}" for o, g, lab in rows])
    ax.set_ylim(len(rows) - 0.5, -0.65)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xticks(np.linspace(0, 1, 5), [f"{v:.0%}" for v in np.linspace(0, 1, 5)])
    ax.set_xlabel("Share of crystals whose copies differ in lattice contacts at the lid")
    notes(fig,
          ("Shows", "copies share the drop, additives and solvent: does the copy with more contacts have the better lid?"),
          ("Means", "the crystal lattice, not the recipe, goes with how much of the lid is seen; mobility shows no rule"),
          ("Caveat", "contacts are counted on modelled atoms, so part of the link is built in"))
    deck.add(fig)


def main():
    deck = Deck(OUT)
    for make in (title_slide, asked_slide, groups_slide,
                 h12_window_slide, midpoint_slide,                       # 28 Oct
                 class_slide, assay_slide,                                # 4 Nov
                 lid_outliers_slide, copies_slide, maps_slide,            # 11 Nov
                 lid_hidden_slide, helix_pairs_slide, h12_mobility_slide,
                 crystal_slide, lattice_slide):
        make(deck)
    deck.close()
    print("wrote", OUT, f"({deck.n} slides)")


if __name__ == "__main__":
    main()
