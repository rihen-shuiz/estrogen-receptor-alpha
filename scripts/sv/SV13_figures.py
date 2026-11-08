"""SV13 figures, drawn only from saved CSV/JSON files (rule 1.8).

1. SV13_mds_map: the shared MDS map (Y03 metric MDS, results/ext/embedding_rmsd.csv mds1/mds2), all 788 chains;
   fill = S08 cluster (cl7 blue, cl4 orange, minor aqua), marker shape and edge = SV12 action class
   (scripts/sv/common/class_style.py; unlabelled = small grey-edged circles). Drawn with the SV03 plotting function
   scripts/sv/common/h12_map.plot_h12_map (x/y = mds1/mds2). Sources: embedding_rmsd.csv, results/sv/SV13/chain_data.csv.
2. SV13_h12_map: the SV03 H12 map (distance H12 centroid to H11 axis vs angle H12-H3; chains with an H12 axis),
   same styling, same function. Sources: results/sv/SV03/geometry_per_chain.csv, results/sv/SV13/chain_data.csv.
3. SV13_forest: class effects of C1 and C2a (chain level in the test-class colour, entry level grey; OR on a log
   axis for the logit outcomes, coefficients otherwise) and the C2b Spearman rho beside them; rows that fail the
   coverage gate carry the gate reason instead of an estimate. Source: results/sv/SV13/tests.csv.
4. SV13_h12_class_panels: one panel per SV12 class (6 classes; empty classes say so) + one overlay of the class
   representatives, all in the one saved view (results/sv/SV13/view.txt): thin semi-transparent C-alpha traces of
   residues 528-550 of every class chain (colour = S08 cluster), the representative's full observed trace (grey,
   matplotlib stand-in for the cartoon: PyMOL not available), its 528-550 trace thick with H12 538-547 coloured by
   the SV01b per-residue B z, and its ligand of interest as sticks. Unobserved residues are absent (line breaks).
   Sources: data/sv/SV13/panel_traces.csv, panel_ligand_atoms.csv, panel_ligand_bonds.csv, panel_meta.json,
   results/sv/SV13/representatives.csv.
Each graph is also saved as its own cropped PNG in figures/panels/sv/SV13/.
Run: .venv/bin/python -B scripts/sv/SV13_figures.py
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
from common.h12_map import GROUP_COLOR, GROUP_LABEL, GROUP_ORDER, plot_h12_map  # noqa: E402
from common.sv_fig import GREY, INK, INK2, MUTED, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV13"
DAT = ROOT / "data/sv/SV13"
STEP = "SV13"
CLASSES = [c for c in CLASS_ORDER if c != "unlabelled"]
H12 = (538, 547)


# ------------------------------------------------------------------------------------------------ maps
def class_layers(ax, df, x, y):
    """Unlabelled first (behind, small), then each labelled class on top; returns legend counts."""
    counts = {}
    u = df[df.cls == "unlabelled"]
    if len(u):
        plot_h12_map(ax, u, x=x, y=y, group_col="group", marker_col="cls", marker_map=CLASS_MARKER, size=9,
                     alpha=0.45, lw=0.3, legend=False, label_axes=False, zorder=2,
                     edge_map={"unlabelled": CLASS_COLOR["unlabelled"]})
        counts["unlabelled"] = len(u)
    for c in CLASSES:
        s = df[df.cls == c]
        if not len(s):
            continue
        plot_h12_map(ax, s, x=x, y=y, group_col="group", marker_col="cls", marker_map=CLASS_MARKER, size=36,
                     alpha=0.95, lw=1.3, legend=False, label_axes=False, zorder=3, edge_map=CLASS_COLOR)
        counts[c] = len(s)
    return counts


def map_legends(ax, df, counts):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    gh = [Patch(color=GROUP_COLOR[g], label=f"{GROUP_LABEL[g]} (n={int((df.group == g).sum())})")
          for g in GROUP_ORDER if (df.group == g).any()]
    leg1 = ax.legend(handles=gh, title="S08 cluster (fill)", loc="upper left", bbox_to_anchor=(1.02, 1.0),
                     fontsize=7, title_fontsize=7, frameon=False, alignment="left")
    ax.add_artist(leg1)
    ch = [Line2D([], [], ls="", marker=CLASS_MARKER[c], markersize=6 if c != "unlabelled" else 4.5,
                 markerfacecolor="white", markeredgecolor=CLASS_COLOR[c], markeredgewidth=1.3,
                 label=f"{CLASS_LABEL[c]} (n={counts.get(c, 0)})") for c in CLASS_ORDER]
    ax.legend(handles=ch, title="action class, SV12 (shape, edge);\nprimary set, chains", loc="upper left",
              bbox_to_anchor=(1.02, 0.72), fontsize=7, title_fontsize=7, frameon=False, alignment="left")


def chain_classes():
    d = pd.read_csv(RES / "chain_data.csv")
    d["cls"] = d.action_class.where(d.primary_set, "unlabelled")
    d["group"] = d.cluster
    return d


def mds_map(plt, d):
    emb = pd.read_csv(ROOT / "results/ext/embedding_rmsd.csv")
    assert (emb.chain_key.values == d.chain_key.values).all()
    df = d.assign(mds1=emb.mds1.values, mds2=emb.mds2.values)
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    fig.subplots_adjust(right=0.70)
    counts = class_layers(ax, df, "mds1", "mds2")
    ax.set_xlabel("MDS 1 (Å)")
    ax.set_ylabel("MDS 2 (Å)")
    n_lab = sum(v for k, v in counts.items() if k != "unlabelled")
    ax.set_title(f"Shared MDS map (Y03 metric MDS of the C-alpha RMSD matrix), {len(df)} chains;\n"
                 f"{n_lab} chains carry an SV12 class: agonist {counts.get('agonist', 0)}, SERM "
                 f"{counts.get('SERM', 0)}, SERD {counts.get('SERD', 0)}", fontsize=9)
    ax.grid(alpha=0.25, lw=0.5)
    map_legends(ax, df, counts)
    return save_figure(fig, STEP, "SV13_mds_map", panels={"mds_map": ax})


def h12_map(plt, d):
    g = pd.read_csv(ROOT / "results/sv/SV03/geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    assert (g.chain_key.values == d.chain_key.values).all()
    df = d[["chain_key", "cls", "group"]].assign(dist_H12c_H11axis_A=g.dist_H12c_H11axis_A.astype(float).values,
                                                 angle_H12_H3_deg=g.angle_H12_H3_deg.astype(float).values)
    df = df[df.dist_H12c_H11axis_A.notna() & df.angle_H12_H3_deg.notna()]
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    fig.subplots_adjust(right=0.70)
    counts = class_layers(ax, df, "dist_H12c_H11axis_A", "angle_H12_H3_deg")
    ax.set_xlabel("distance H12 centroid to H11 axis (Å)")
    ax.set_ylabel("angle H12–H3 (°)")
    tot = chain_classes().groupby("cls").size()
    ax.set_title(f"H12 map (SV03), {len(df)} of {len(d)} chains with an H12 axis (538-547, ≥ 6 helical Cα);\n"
                 + ", ".join(f"{c} {counts.get(c, 0)} of {int(tot.get(c, 0))}" for c in ("agonist", "SERM", "SERD")),
                 fontsize=9)
    ax.grid(alpha=0.25, lw=0.5)
    map_legends(ax, df, counts)
    return save_figure(fig, STEP, "SV13_h12_map", panels={"h12_map": ax})


# ------------------------------------------------------------------------------------------------ forest
COMP_LABEL = {"C1_SERD_vs_SERM": "C1 SERD vs SERM", "C1_antagonist_vs_SERM": "C1 antagonist vs SERM",
              "C2a_partial_agonist_vs_agonist": "C2a partial agonist vs agonist",
              "C2b_agonist_rho_dpX": "C2b agonists: ρ with ΔpX"}
STRATUM_LABEL = {"all": "all chains", "cl4": "within cl4", "cl7": "within cl7"}


def forest_axis(ax, t, null, logx, xlabel, title, short=False, wrap=90):
    import textwrap
    rows = t.drop_duplicates(["comparison", "stratum", "outcome"])[["comparison", "stratum", "outcome",
                                                                     "outcome_label", "test_class"]]
    n = len(rows)
    ax.set_ylim(n - 0.4, -0.6)
    if logx:
        ax.set_xscale("log")
        ax.set_xlim(0.05, 20)
    else:
        ax.set_xlim(-1.05, 1.05) if null == 0 and "ρ" in xlabel else ax.set_xlim(-3, 3)
    ax.axvline(null, color=MUTED, lw=0.8, zorder=1)
    labels = []
    for k, r in enumerate(rows.itertuples()):
        sub = t[(t.comparison == r.comparison) & (t.stratum == r.stratum) & (t.outcome == r.outcome)]
        labels.append(r.outcome_label if short else
                      f"{COMP_LABEL[r.comparison]} | {STRATUM_LABEL[r.stratum]} | {r.outcome_label}")
        if k % 2 == 0:
            ax.axhspan(k - 0.5, k + 0.5, color=GREY, alpha=0.08, lw=0, zorder=0)
        if not sub.gate_pass.iloc[0]:
            ax.text(0.5, k, textwrap.fill("not testable: " + str(sub.gate_reason.iloc[0]), wrap),
                    transform=ax.get_yaxis_transform(),
                    ha="center", va="center", fontsize=6.6, color=INK2)
            continue
        for lv, off, col in (("chain", -0.15, CLASS_COLOR.get(r.test_class, INK)), ("entry", 0.15, GREY)):
            s = sub[sub.level == lv]
            if not len(s) or not np.isfinite(s.estimate.iloc[0]):
                ax.text(0.98, k + off, f"{lv}: {s.status.iloc[0] if len(s) else 'missing'}",
                        transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6, color=INK2)
                continue
            e, lo, hi, q = s.estimate.iloc[0], s.ci_low.iloc[0], s.ci_high.iloc[0], s.q.iloc[0]
            if np.isfinite(lo):
                ax.plot([lo, hi], [k + off] * 2, color=col, lw=1.4, zorder=3)
            ax.scatter([e], [k + off], color=col, s=22, zorder=4)
            ax.text(0.99, k + off, f"{lv} {e:.2f}, q {q:.2g}", transform=ax.get_yaxis_transform(), ha="right",
                    va="center", fontsize=6, color=col)
    ax.set_yticks(range(n), labels, fontsize=6.8)
    ax.set_xlabel(xlabel, fontsize=7.5)
    ax.set_title(title, fontsize=8.2)
    ax.tick_params(axis="x", labelsize=7)


def forest(plt, t):
    from matplotlib.lines import Line2D
    p = t[t.set == "primary"]
    a = p[p.part.isin(["C1"]) & (p.family == "logit")]
    b = p[(p.part.isin(["C1"]) & (p.family == "linear")) | (p.part == "C2a")]
    c = p[p.part == "C2b"]
    n_rows = p.drop_duplicates(["comparison", "stratum", "outcome"])
    n_test = int(n_rows.gate_pass.sum())
    fig = plt.figure(figsize=(15.5, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0], height_ratios=[len(a.drop_duplicates(
        ["comparison", "stratum", "outcome"])), len(b.drop_duplicates(["comparison", "stratum", "outcome"]))],
        wspace=0.42, hspace=0.45)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[0, 1])
    forest_axis(ax1, a, 1.0, True, "odds ratio, class vs reference (log scale; Q + X + mutations + peptide)",
                "C1: H12 complete and fraction observed")
    forest_axis(ax2, b, 0.0, False, "coefficient, class vs reference (outcome units; same covariates)",
                "C1: H12 mobility z; C2a (within cl7): mobility, H12 geometry, ψ 537")
    forest_axis(ax3, c, 0.0, False, "Spearman ρ (outcome vs SV15 ΔpX)",
                "C2b (within cl7): agonists, Spearman ρ of each outcome with ΔpX", short=True, wrap=62)
    h = [Line2D([], [], color=INK2, marker="o", lw=1.4, label="chain level (entry-clustered errors), test-class "
                "colour"), Line2D([], [], color=GREY, marker="o", lw=1.4, label="entry level (one chain per entry)")]
    ax3.legend(handles=h, loc="upper left", bbox_to_anchor=(0.0, -0.30), fontsize=7, frameon=False)
    fig.suptitle(f"SV13 Part B, class effects on H12 (primary set): {n_test} of {len(n_rows)} comparison x stratum x "
                 f"outcome rows pass the coverage gate (≥ 5 entries per side; SV12 classes: agonist 27, SERM 20, "
                 f"SERD 4, partial agonist 0, antagonist 0 entries)", fontsize=9.5, x=0.01, ha="left", y=0.995)
    return save_figure(fig, STEP, "SV13_forest", panels={"C1_logit": ax1, "C1_mobility_C2a": ax2, "C2b_rho": ax3})


# ------------------------------------------------------------------------------------------------ H12 panels
def line_nan(df):
    """x/y arrays with NaN at unobserved residues so that matplotlib breaks the line there."""
    df = df.sort_values("position")
    x = np.where(df.observed, df.sx, np.nan)
    y = np.where(df.observed, df.sy, np.nan)
    return x, y


def draw_ligand(ax, la, lb, chain, color, lw=1.7, z=8):
    a = la[la.chain_key == chain].set_index("atom_index")
    b = lb[lb.chain_key == chain]
    for r in b.itertuples():
        ax.plot([a.sx[r.i], a.sx[r.j]], [a.sy[r.i], a.sy[r.j]], color=color, lw=lw, solid_capstyle="round",
                zorder=z)
    ax.scatter(a.sx, a.sy, s=4, color=color, zorder=z + 0.1, lw=0)


def scale_bar(ax, lim):
    x1 = lim["xlim"][1] - 3
    y0 = lim["ylim"][0] + 3
    ax.plot([x1 - 10, x1], [y0, y0], color=INK, lw=1.5)
    ax.text(x1 - 5, y0 + 1.0, "10 Å", ha="center", va="bottom", fontsize=6.5, color=INK)


def base_axes(ax, lim):
    ax.set_xlim(*lim["xlim"])
    ax.set_ylim(*lim["ylim"])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color(MUTED)
        s.set_linewidth(0.6)
    scale_bar(ax, lim)


def class_panel(ax, cls, tr, la, lb, rep, lim, cmap, norm):
    import matplotlib.patheffects as pe
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D
    base_axes(ax, lim)
    r = rep.set_index("action_class").loc[cls]
    lab = CLASS_LABEL[cls]
    if int(r.n_chains) == 0:
        ax.text(0.5, 0.5, f"no chain carries this class\nin SV12 (0 chains, 0 entries)", transform=ax.transAxes,
                ha="center", va="center", fontsize=8, color=INK2)
        ax.set_title(f"{lab}: 0 chains", fontsize=8.5, color=CLASS_COLOR[cls])
        return None
    t = tr[tr.action_class == cls]
    rk = r.representative
    full = t[(t.role == "rep_full")]
    x, y = line_nan(full)
    ax.plot(x, y, color=GREY, lw=2.0, alpha=0.45, zorder=2, solid_capstyle="round")
    ct = t[t.role == "class_trace"]
    ncl = {}
    for (ck, cl), g in ct.groupby(["chain_key", "cluster"]):
        if ck == rk:
            continue
        x, y = line_nan(g)
        ax.plot(x, y, color=GROUP_COLOR[cl], lw=0.7, alpha=0.4, zorder=3)
    for cl, g in ct.groupby("cluster"):
        ncl[cl] = g.chain_key.nunique()
    rw = ct[ct.chain_key == rk].sort_values("position")
    x, y = line_nan(rw)
    ax.plot(x, y, color=INK, lw=2.4, zorder=5, solid_capstyle="round")
    h = rw[rw.position.between(*H12)]
    segs, vals = [], []
    hp = h.set_index("position")
    for p in range(H12[0], H12[1]):
        if p in hp.index and p + 1 in hp.index and hp.observed[p] and hp.observed[p + 1]:
            segs.append([(hp.sx[p], hp.sy[p]), (hp.sx[p + 1], hp.sy[p + 1])])
            vals.append(np.nanmean([hp.b_z[p], hp.b_z[p + 1]]))
    if segs:
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidths=3.6, zorder=6, capstyle="round",
                            path_effects=[pe.Stroke(linewidth=5.2, foreground=INK), pe.Normal()])
        lc.set_array(np.asarray(vals))
        ax.add_collection(lc)
    single = [p for p in hp.index if hp.observed[p] and not ((p - 1 in hp.index and hp.observed.get(p - 1, False))
                                                             or (p + 1 in hp.index and hp.observed.get(p + 1, False)))]
    for p in single:
        ax.scatter([hp.sx[p]], [hp.sy[p]], c=[hp.b_z[p]], cmap=cmap, norm=norm, s=30, edgecolors=INK, zorder=6)
    draw_ligand(ax, la, lb, rk, CLASS_COLOR[cls])
    n_h12 = int(h.observed.sum())
    ax.set_title(f"{lab}: {int(r.n_chains)} chains, {int(r.n_entries)} entries, {int(r.n_complete_h12)} with "
                 f"complete H12 (538-547)", fontsize=8.2, color=CLASS_COLOR[cls])
    ax.text(0.01, 0.99, f"representative {rk} ({r.rep_ligand}, {r.rep_cluster}): median Cα RMSD "
            f"{r.rep_median_rmsd_A:.2f} Å\nto the other {int(r.n_chains) - 1} chains; H12 538-547: {n_h12} of 10 "
            f"residues modelled", transform=ax.transAxes, ha="left", va="top", fontsize=6.4, color=INK2)
    hl = [Line2D([], [], color=GROUP_COLOR[c], lw=1.0, label=f"{GROUP_LABEL[c]}: {ncl[c]} chain{'s' if ncl[c] != 1 else ''}, 528-550")
          for c in GROUP_ORDER if c in ncl]
    hl += [Line2D([], [], color=GREY, lw=2.0, alpha=0.6, label="representative, full Cα trace"),
           Line2D([], [], color=INK, lw=2.4, label="representative 528-550 (H12 by B z)"),
           Line2D([], [], color=CLASS_COLOR[cls], lw=1.7, label=f"ligand {r.rep_ligand} (sticks)")]
    ax.legend(handles=hl, loc="lower left", fontsize=5.8, frameon=False, handlelength=1.6)
    cax = ax.inset_axes([0.60, 0.05, 0.22, 0.025])
    import matplotlib.cm as cm
    cb = ax.figure.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_label("H12 Cα B z (SV01b)", fontsize=5.8, labelpad=1)
    cb.ax.tick_params(labelsize=5.5, length=2, pad=1)
    return cax


def overlay_panel(ax, tr, la, lb, rep, lim):
    from matplotlib.lines import Line2D
    base_axes(ax, lim)
    hl = []
    for r in rep[rep.n_chains > 0].itertuples():
        t = tr[(tr.action_class == r.action_class) & (tr.chain_key == r.representative)]
        x, y = line_nan(t[t.role == "rep_full"])
        ax.plot(x, y, color=CLASS_COLOR[r.action_class], lw=1.0, alpha=0.35, zorder=2)
        w = t[t.role == "class_trace"]
        x, y = line_nan(w)
        ax.plot(x, y, color=CLASS_COLOR[r.action_class], lw=2.8, zorder=5, solid_capstyle="round")
        draw_ligand(ax, la, lb, r.representative, CLASS_COLOR[r.action_class], lw=1.4)
        n_h12 = int(w[w.position.between(*H12)].observed.sum())
        hl.append(Line2D([], [], color=CLASS_COLOR[r.action_class], lw=2.8,
                         label=f"{CLASS_LABEL[r.action_class]}: {r.representative} ({r.rep_ligand}, {r.rep_cluster}; "
                               f"H12 {n_h12}/10)"))
    hl.append(Line2D([], [], color=INK2, lw=1.0, alpha=0.5, label="thin: full Cα trace; thick: 528-550; sticks: ligand"))
    ax.legend(handles=hl, loc="lower left", fontsize=5.8, frameon=False, handlelength=1.6)
    ax.set_title(f"Overlay of the {int((rep.n_chains > 0).sum())} class representatives (same view)", fontsize=8.2)


def class_panels(plt):
    from matplotlib.colors import Normalize
    tr = pd.read_csv(DAT / "panel_traces.csv")
    la = pd.read_csv(DAT / "panel_ligand_atoms.csv")
    lb = pd.read_csv(DAT / "panel_ligand_bonds.csv")
    meta = json.loads((DAT / "panel_meta.json").read_text())
    rep = pd.read_csv(RES / "representatives.csv")
    lim = meta["limits"]
    cmap = plt.get_cmap("viridis")
    norm = Normalize(*meta["colour_range_b_z"])
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.6))
    order = ["agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other"]
    panels = {}
    for ax, cls in zip(axes.ravel()[:6], order):
        cax = class_panel(ax, cls, tr, la, lb, rep, lim, cmap, norm)
        panels[cls] = [ax, cax] if cax is not None else ax
    overlay_panel(axes.ravel()[6], tr, la, lb, rep, lim)
    panels["overlay"] = axes.ravel()[6]
    lx = axes.ravel()[7]
    lx.axis("off")
    v = meta["view"]
    lx.text(0.0, 0.98, "One view for every panel (results/sv/SV13/view.txt):\n"
            f"principal axes of {v['view_chain']} (cl7 medoid) Cα at the\ncore v3 positions + H12 538-547, "
            "core v3 frame;\northographic projection, same limits everywhere.\n\n"
            "Class chains = SV12 primary-set chains (all clusters).\nRepresentative = lowest median Cα RMSD "
            "(existing\nmatrix results/rmsd_entry.npy) to the other chains of\nits class; ties by resolution.\n\n"
            "Unmodelled residues are absent (line breaks):\nno rebuilding, no interpolation.\n\n"
            "Fallback without PyMOL (step 12): Cα traces drawn\nwith matplotlib; the representative's cartoon is a\n"
            "thick Cα trace, its ligand of interest sticks.\n\n"
            "Crystal structures show H12 order and position;\nthey say nothing about receptor degradation or\n"
            "cellular activity.", transform=lx.transAxes, ha="left", va="top", fontsize=7.2, color=INK2)
    fig.suptitle("H12 (528-550) and the ligand of interest per SV12 action class, superposed on core v3 "
                 "(H12 window 538-547)", fontsize=10.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout()
    return save_figure(fig, STEP, "SV13_h12_class_panels", panels=panels)


def main():
    plt = mpl_setup()
    d = chain_classes()
    t = pd.read_csv(RES / "tests.csv", keep_default_na=False, na_values=[""])
    t["gate_pass"] = t.gate_pass.astype(str).eq("True")
    out = mds_map(plt, d) + h12_map(plt, d) + forest(plt, t) + class_panels(plt)
    for p in out:
        print(p.relative_to(ROOT), f"{p.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main()
