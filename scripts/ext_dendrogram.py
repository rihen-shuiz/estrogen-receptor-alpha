"""Tracker Y07: S08 dendrogram with metadata colour bars and the leaf-ordered RMSD heatmap with the same sidebars.

Inputs: results/linkage_entry.npy (S08 average linkage of results/rmsd_entry.npy), results/ext/chain_metadata.csv
(Y02), chain_table() (tiers, B-group, H12 label, R-free).
Bars (top to bottom): S08 cluster (palette), QC tier, R-free (sequential), B-group, H12 label, peptide class,
construct mutation, space group (5 most common + other), crystallisation pH (sequential).
Colours: the cluster bar uses the shared cluster palette; every other categorical bar uses palette slots 3-8 and grey
so blue/orange always mean cl7/cl4.
Outputs: results/ext/leaf_order.csv, figures/ext/Y07_dendrogram_metadata.png, figures/ext/Y07_rmsd_heatmap_annotated.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_dendrogram.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, to_rgb
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import dendrogram

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, CLUSTER_COLOR, GREEN, GREY, INK, MAGENTA, RED, RES, ROOT, SEQ_BLUE, VIOLET, YELLOW,
                        chain_table, mpl_setup, rmsd, save_fig, write_numbers)

OTHER = [AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]


def cat_colors(values, order=None, fixed=None):
    vals = pd.Series(values).fillna("not_reported").astype(str)
    order = order or list(vals.value_counts().index)
    fixed = fixed or {}
    cmap, k = {}, 0
    for v in order:
        if v in fixed:
            cmap[v] = fixed[v]
        elif v in ("not_reported", "unknown", "nan", "other", "none"):
            cmap[v] = GREY if v != "none" else "#e4e3df"
        else:
            cmap[v] = OTHER[k % len(OTHER)]
            k += 1
    return [to_rgb(cmap[v]) for v in vals], cmap


def seq_colors(values, vmin, vmax):
    import matplotlib
    cm = matplotlib.colormaps["esr_seq"]
    v = np.asarray(values, float)
    out = []
    for x in v:
        out.append(to_rgb(GREY) if np.isnan(x) else cm((np.clip(x, vmin, vmax) - vmin) / (vmax - vmin))[:3])
    return out


def bars(t):
    sg = t.space_group.fillna("not_reported")
    top = sg.value_counts().index[:5]
    sg = sg.where(sg.isin(top), "other")
    defs = []
    c = [to_rgb(CLUSTER_COLOR[int(x)]) for x in t.cluster]
    defs.append(("S08 cluster", c, {"cl7": CLUSTER_COLOR[7], "cl4": CLUSTER_COLOR[4], "minor": AQUA}))
    col, m = cat_colors(t.tier, ["A", "B", "C"], {"A": GREEN, "B": YELLOW, "C": RED})
    defs.append(("QC tier", col, m))
    defs.append(("R-free (0.18-0.32)", seq_colors(t.r_free, 0.18, 0.32), None))
    col, m = cat_colors(t.b_group, ["low", "mid", "high"], {"low": GREEN, "mid": YELLOW, "high": RED})
    defs.append(("B-group", col, m))
    col, m = cat_colors(t.h12_label, ["complete", "partial", "not_modelled"],
                        {"complete": VIOLET, "partial": YELLOW, "not_modelled": RED})
    defs.append(("H12 label", col, m))
    col, m = cat_colors(t.peptide_class)
    defs.append(("peptide", col, m))
    col, m = cat_colors(t.mutation_status)
    defs.append(("construct mutation", col, m))
    col, m = cat_colors(sg)
    defs.append(("space group", col, m))
    defs.append(("pH (5.5-9)", seq_colors(t.ph, 5.5, 9.0), None))
    return defs


def table():
    t = chain_table()
    cm = pd.read_csv(RES / "chain_metadata.csv")
    t = t.merge(cm[["chain_key", "mutation_status", "peptide_class", "ph"]], on="chain_key", how="left",
                validate="one_to_one")
    return t


def main():
    d, idx = rmsd()
    t = table()
    z = np.load(ROOT / "results/linkage_entry.npy")
    plt = mpl_setup()
    fig = plt.figure(figsize=(18, 9.5))
    defs = bars(t)
    nb = len(defs)
    gs = fig.add_gridspec(2 + nb, 1, height_ratios=[5] + [0.32] * nb + [0.05], hspace=0.08, bottom=0.17, top=0.95)
    ax = fig.add_subplot(gs[0])
    with plt.rc_context({"lines.linewidth": 0.6}):
        dn = dendrogram(z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color="#6b6a65")
    ax.axhline(1.5, color=RED, ls="--", lw=1)
    ax.text(ax.get_xlim()[1], 1.5, " 1.5 A cut", color=RED, va="center", fontsize=8)
    ax.set_ylabel("average-linkage height (RMSD, A)")
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    leaves = dn["leaves"]
    pd.DataFrame({"leaf_rank": np.arange(len(leaves)), "row": leaves, "chain_key": idx.chain_key.values[leaves],
                  "cluster": t.cluster.values[leaves]}).to_csv(RES / "leaf_order.csv", index=False)
    handles = []
    for k, (name, col, cmap) in enumerate(defs):
        a = fig.add_subplot(gs[1 + k])
        a.imshow(np.array(col)[leaves][None, :, :], aspect="auto", interpolation="nearest")
        a.set_yticks([0], [name], fontsize=8)
        a.set_xticks([])
        for s in a.spines.values():
            s.set_visible(False)
        if cmap:
            handles.append((name, cmap))
    # legends below
    x0 = 0.01
    for name, cmap in handles:
        hs = [Patch(color=v, label=k) for k, v in list(cmap.items())[:7]]
        leg = fig.legend(handles=hs, title=name, loc="lower left", bbox_to_anchor=(x0, 0.0), fontsize=7,
                         title_fontsize=7.5, ncol=1, frameon=False)
        x0 += 0.125
    fig.suptitle("Dendrogram of 788 chains: the two main branches line up with peptide, construct mutation and "
                 "H12 label, not with tier, R-free or B-group", x=0.01, ha="left", fontsize=11, fontweight="bold")
    save_fig(fig, "Y07_dendrogram_metadata", section="A. Structural clustering zoo", task="Y07",
             shows="S08 average-linkage dendrogram (788 leaves) with nine metadata colour bars",
             takeaway="The cl7/cl4 branches coincide with coactivator peptide (present/absent), Y537S vs L536S "
                      "constructs and complete vs partial H12, while quality bars (tier, R-free, B-group) are mixed "
                      "along both branches.",
             source="results/linkage_entry.npy, results/ext/chain_metadata.csv, results/ext/leaf_order.csv",
             script="scripts/ext_dendrogram.py")

    # ---------------- heatmap
    fig = plt.figure(figsize=(13, 12))
    sel = [defs[i] for i in (0, 1, 4, 5, 6)]
    gs = fig.add_gridspec(1 + len(sel) + 1, 1 + len(sel) + 1,
                          height_ratios=[0.35] * len(sel) + [10, 0.4], width_ratios=[0.35] * len(sel) + [10, 0.35],
                          hspace=0.03, wspace=0.03)
    ax = fig.add_subplot(gs[len(sel), len(sel)])
    dm = d[np.ix_(leaves, leaves)]
    im = ax.imshow(dm, cmap="esr_seq", vmin=0, vmax=6, interpolation="nearest", aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("788 chains in dendrogram leaf order")
    for k, (name, col, cmap) in enumerate(sel):
        a = fig.add_subplot(gs[k, len(sel)])
        a.imshow(np.array(col)[leaves][None, :, :], aspect="auto", interpolation="nearest")
        a.set_xticks([]); a.set_yticks([0], [name], fontsize=7.5)
        a.yaxis.tick_right()
        b = fig.add_subplot(gs[len(sel), k])
        b.imshow(np.array(col)[leaves][:, None, :], aspect="auto", interpolation="nearest")
        b.set_yticks([]); b.set_xticks([0], [name], rotation=90, fontsize=7.5)
        for s in list(a.spines.values()) + list(b.spines.values()):
            s.set_visible(False)
    cax = fig.add_subplot(gs[len(sel) + 1, len(sel)])
    fig.colorbar(im, cax=cax, orientation="horizontal", label="C-alpha RMSD (A), clipped at 6")
    hs = []
    for name, col, cmap in sel:
        if cmap:
            hs += [Patch(color=v, label=f"{name}: {k}") for k, v in list(cmap.items())[:5]]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.04, 0.965), fontsize=6.5, ncol=1)
    cl = t.cluster.values
    fig.suptitle(f"RMSD matrix in leaf order: two blocks (cl7 {np.sum(cl == 7)}, cl4 {np.sum(cl == 4)} chains) with a "
                 f"median of {np.median(d[np.ix_(cl == 7, cl == 4)]):.1f} A between them", x=0.01, ha="left",
                 fontsize=11, fontweight="bold", y=0.995)
    save_fig(fig, "Y07_rmsd_heatmap_annotated", section="A. Structural clustering zoo", task="Y07",
             shows="788 x 788 RMSD heatmap ordered by dendrogram leaves with cluster, tier, H12, peptide and mutation sidebars",
             takeaway="Two dark diagonal blocks; the cl4 block has visible sub-blocks (more internal heterogeneity) "
                      "and the off-diagonal is uniformly ~3.5 A.",
             source="results/rmsd_entry.npy, results/ext/leaf_order.csv", script="scripts/ext_dendrogram.py")
    write_numbers("Y07", dict(n_leaves=len(leaves)))


if __name__ == "__main__":
    main()
