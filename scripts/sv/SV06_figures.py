"""SV06 figures, drawn from saved CSVs only (rule 1.8).

1. SV06_outliers_map (two graphs):
   A  the shared MDS map (Y03 metric MDS, results/ext/embedding_rmsd.csv mds1/mds2): all 788 chains in the cluster
      colours (cl7 blue, cl4 orange, minor clusters aqua), the 22 Y16 consensus outliers on top, coloured by the Y16
      reason (quality red, conformation violet, the Y16 colours), marker = cluster group, labelled PDB|author chain;
      the two complete-H12 main medoids (Y10) as stars.
   B  per outlier, RMSD to its own main-cluster medoid (x) vs RMSD to the other main-cluster medoid (y), y = x drawn.
   Source: results/sv/SV06/outliers_on_map.csv, results/sv/SV06/numbers.json, results/ext/embedding_rmsd.csv,
   results/sv/SV00/chain_table.csv (cluster colour of the background chains).
2. SV06_rmsd_heatmap_outlier_bar: the Y07 annotated RMSD heatmap redrawn (same leaf order, colour scale and the
   same five side bars, built by scripts/ext_dendrogram.py bars()) with one extra side bar, 'Y16 outlier', that
   marks the 22 consensus outliers in the Y16 reason colours.
   Source: results/rmsd_entry.npy, results/ext/leaf_order.csv, results/sv/SV06/outliers_on_map.csv, the Y07 bar
   inputs (chain_table() of scripts/ext_common.py, results/ext/chain_metadata.csv).
Nothing is written outside figures/sv/SV06 and figures/panels/sv/SV06 (no figure card in results/ext).
Run: .venv/bin/python -B -I scripts/sv/SV06_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.sv_fig import (AQUA, BLUE, GREY, INK, INK2, MUTED, ORANGE, RED, VIOLET,  # noqa: E402
                           mpl_setup, save_figure)
import ext_dendrogram as y07  # noqa: E402  (scripts/ is on sys.path via sv_fig; only bars()/table() are used)

RES = ROOT / "results/sv/SV06"
STEP = "SV06"
REASON_COLOR = {"quality": RED, "conformation": VIOLET}          # the Y16 verdict colours
GROUP_MARKER = {"cl7": "o", "cl4": "s", "minor": "^"}
GROUP_COLOR = {"cl7": BLUE, "cl4": ORANGE, "minor": AQUA}
NO_OUTLIER = "#f4f3f0"


def place_labels(ax, xs, ys, texts, avoid_xy, fontsize=7.0, start=16.0, iters=200, color=INK, avoid_axes=()):
    """Place one label per point without overlaps (simple repulsion in pixel space) and draw leader lines.
    avoid_axes: inset axes whose area the labels must stay out of."""
    if len(xs) == 0:
        return
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [a.get_window_extent(renderer) for a in avoid_axes]
    tr = ax.transData
    P = tr.transform(np.c_[xs, ys])
    A = tr.transform(np.asarray(avoid_xy, float))
    probe = [ax.text(0, 0, t, fontsize=fontsize) for t in texts]
    wh = np.array([[o.get_window_extent(renderer).width, o.get_window_extent(renderer).height] for o in probe])
    for o in probe:
        o.remove()
    wh += 3.0
    box = ax.get_window_extent(renderer)
    r_pt = 7.0

    def overlap(c, w, h, x0, y0, x1, y1):
        ox = min(c[0] + w / 2, x1) - max(c[0] - w / 2, x0)
        oy = min(c[1] + h / 2, y1) - max(c[1] - h / 2, y0)
        return max(ox, 0.0) * max(oy, 0.0)

    def seg_dist(a, b, q):
        ab = b - a
        tt = np.clip(np.dot(q - a, ab) / max(np.dot(ab, ab), 1e-9), 0, 1)
        return np.linalg.norm(a + tt * ab - q)

    def crosses(a, b, c, d):
        def orient(p, q, r):
            return np.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]))
        return orient(a, b, c) != orient(a, b, d) and orient(c, d, a) != orient(c, d, b)

    # greedy start: for each label (most crowded first) the candidate position around its point with the
    # least overlap with placed labels, points and insets, and no leader through another point or leader
    crowd = np.array([(np.linalg.norm(A - q, axis=1) < 40).sum() for q in P])
    L = P.copy()
    placed = []
    for i in np.argsort(-crowd, kind="stable"):
        w, h = wh[i]
        best, best_s = P[i] + [start + w / 2, 0], np.inf
        for dist in (start, start * 2, start * 3.2):
            for ang in np.deg2rad(np.arange(0, 360, 30)):
                u = np.array([np.cos(ang), np.sin(ang)])
                c = P[i] + u * (dist + np.abs(u) @ np.array([w / 2, h / 2]))
                sc = dist * 0.5
                if c[0] - w / 2 < box.x0 or c[0] + w / 2 > box.x1 or c[1] - h / 2 < box.y0 or c[1] + h / 2 > box.y1:
                    sc += 1e6
                for j in placed:
                    sc += 10 * overlap(c, w, h, L[j, 0] - wh[j, 0] / 2, L[j, 1] - wh[j, 1] / 2,
                                       L[j, 0] + wh[j, 0] / 2, L[j, 1] + wh[j, 1] / 2)
                    if crosses(P[i], c, P[j], L[j]):
                        sc += 300
                for q in A:
                    sc += 10 * overlap(c, w, h, q[0] - r_pt, q[1] - r_pt, q[0] + r_pt, q[1] + r_pt)
                    if np.linalg.norm(q - P[i]) > 1e-6 and seg_dist(P[i], c, q) < r_pt:
                        sc += 300
                for bb in boxes:
                    sc += 10 * overlap(c, w, h, bb.x0, bb.y0, bb.x1, bb.y1)
                if sc < best_s:
                    best, best_s = c, sc
        L[i] = best
        placed.append(i)
    for _ in range(iters):
        F = np.zeros_like(L)
        n = len(L)
        for i in range(n):
            for j in range(i + 1, n):
                d = L[j] - L[i]
                ox = (wh[i, 0] + wh[j, 0]) / 2 - abs(d[0])
                oy = (wh[i, 1] + wh[j, 1]) / 2 - abs(d[1])
                if ox > 0 and oy > 0:
                    if ox < oy:
                        s = np.sign(d[0]) or 1.0
                        F[i, 0] -= s * ox / 2; F[j, 0] += s * ox / 2
                    else:
                        s = np.sign(d[1]) or 1.0
                        F[i, 1] -= s * oy / 2; F[j, 1] += s * oy / 2
            for a in A:
                d = L[i] - a
                ox = wh[i, 0] / 2 + r_pt - abs(d[0])
                oy = wh[i, 1] / 2 + r_pt - abs(d[1])
                if ox > 0 and oy > 0:
                    if ox < oy:
                        F[i, 0] += (np.sign(d[0]) or 1.0) * ox * 0.6
                    else:
                        F[i, 1] += (np.sign(d[1]) or 1.0) * oy * 0.6
            for bb in boxes:
                ox1 = L[i, 0] + wh[i, 0] / 2 - bb.x0
                ox2 = bb.x1 - (L[i, 0] - wh[i, 0] / 2)
                oy1 = L[i, 1] + wh[i, 1] / 2 - bb.y0
                oy2 = bb.y1 - (L[i, 1] - wh[i, 1] / 2)
                if min(ox1, ox2, oy1, oy2) > 0:
                    m = int(np.argmin([ox1, ox2, oy1, oy2]))
                    F[i] += [(-ox1, 0), (ox2, 0), (0, -oy1), (0, oy2)][m]
            # weak pull back towards the anchor when far away
            d = P[i] - L[i]
            dist = np.linalg.norm(d)
            if dist > start * 4:
                F[i] += d / dist * 0.5
        if np.abs(F).max() < 0.05:
            break
        L += F
        L[:, 0] = np.clip(L[:, 0], box.x0 + wh[:, 0] / 2, box.x1 - wh[:, 0] / 2)
        L[:, 1] = np.clip(L[:, 1], box.y0 + wh[:, 1] / 2, box.y1 - wh[:, 1] / 2)
    LD = tr.inverted().transform(L)
    for (x, y), (lx, ly), t in zip(np.c_[xs, ys], LD, texts):
        ax.annotate(t, xy=(x, y), xytext=(lx, ly), ha="center", va="center", fontsize=fontsize, color=color,
                    zorder=6, arrowprops=dict(arrowstyle="-", lw=0.5, color=INK2, shrinkA=1, shrinkB=3))


def dense_mask(x, y, r=0.35, k=3):
    """Points with at least k other outliers within r (plot units): drawn and labelled in a zoom inset."""
    P = np.c_[np.asarray(x, float), np.asarray(y, float)]
    D = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    return (D < r).sum(1) - 1 >= k


def zoom_inset(ax, o, x, y, dense, rect, background=None, margin=0.22):
    """Inset that zooms on the dense outliers, labels them and marks the zoom window on the parent axes."""
    ins = ax.inset_axes(rect)
    inside = pd.Series(np.asarray(dense, bool), index=o.index)
    while True:                       # grow the window until every outlier inside it has a margin around it
        sub = o[inside]
        mx = margin * (sub[x].max() - sub[x].min())
        my = margin * (sub[y].max() - sub[y].min())
        x0, x1 = sub[x].min() - mx, sub[x].max() + mx
        y0, y1 = sub[y].min() - my, sub[y].max() + my
        grown = (o[x] >= x0) & (o[x] <= x1) & (o[y] >= y0) & (o[y] <= y1)
        if (grown == inside).all():
            break
        inside = grown | inside
    if background is not None:
        bx_, by_, bc = background
        m = (bx_ > x0) & (bx_ < x1) & (by_ > y0) & (by_ < y1)
        ins.scatter(bx_[m], by_[m], s=12, c=np.asarray(bc)[m], alpha=0.45, lw=0, zorder=2)
    outlier_scatter(ins, o, x, y, size=60)
    ins.set_xlim(x0, x1)
    ins.set_ylim(y0, y1)
    ins.tick_params(labelsize=6.5)
    for sp in ins.spines.values():
        sp.set_visible(True); sp.set_color(MUTED)
    ins.set_facecolor("white")
    ax.indicate_inset_zoom(ins, edgecolor=MUTED, alpha=0.9)
    place_labels(ins, list(sub[x]), list(sub[y]), list(sub.label), avoid_xy=np.c_[sub[x], sub[y]],
                 fontsize=6.8, start=12)
    ins.set_title(f"zoom: {int(inside.sum())} outliers", fontsize=7, fontweight="normal", color=INK2)
    return ins, inside.values


def outlier_scatter(ax, o, x, y, size=46):
    for g, mk in GROUP_MARKER.items():
        for rsn, col in REASON_COLOR.items():
            m = (o.cluster_group == g) & (o.reason == rsn)
            if m.any():
                ax.scatter(o.loc[m, x], o.loc[m, y], s=size, marker=mk, c=col, edgecolor=INK, lw=0.7, zorder=5)


def map_figure(o, nums, plt):
    from matplotlib.lines import Line2D
    emb = pd.read_csv(ROOT / "results/ext/embedding_rmsd.csv")
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", usecols=["chain_key", "s08_cluster", "auth_asym_id"])
    assert (emb.chain_key.values == ct.chain_key.values).all()
    grp = ct.s08_cluster.map(lambda c: {7: "cl7", 4: "cl4"}.get(int(c), "minor"))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(15.5, 7.4), gridspec_kw=dict(width_ratios=[1.12, 1], wspace=0.18))

    # ---------------- A: MDS map
    for g in ("cl7", "cl4", "minor"):
        m = (grp == g).values
        ax.scatter(emb.mds1[m], emb.mds2[m], s=9 if g != "minor" else 16, c=GROUP_COLOR[g],
                   marker="o" if g != "minor" else "^", alpha=0.35 if g != "minor" else 0.6, lw=0, zorder=2)
    outlier_scatter(ax, o, "mds1", "mds2")
    medk = [nums["medoid_cl7"], nums["medoid_cl4"]]
    me = emb.set_index("chain_key").loc[medk, ["mds1", "mds2"]]
    lab = ct.set_index("chain_key").loc[medk]
    med_pdb = [f"{k.split('|')[0]}|{lab.loc[k, 'auth_asym_id']}" for k in medk]
    med_lbl = [f"cl7 medoid {med_pdb[0]}", f"cl4 medoid {med_pdb[1]}"]
    ax.scatter(me.mds1, me.mds2, s=170, marker="*", c=[BLUE, ORANGE], edgecolor=INK, lw=0.8, zorder=5)
    x0, x1 = emb.mds1.min(), emb.mds1.max()
    y0, y1 = emb.mds2.min(), emb.mds2.max()
    ax.set_xlim(x0 - 0.6, x1 + 1.9)
    ax.set_ylim(y0 - 0.6, y1 + 1.0)
    ax.set_xlabel("MDS 1 (metric MDS of the 788 x 788 RMSD matrix, Y03)")
    ax.set_ylabel("MDS 2")
    ax.set_aspect("equal", adjustable="box")
    dA = dense_mask(o.mds1, o.mds2)
    bcol = grp.map(GROUP_COLOR).values
    insA, inA = zoom_inset(ax, o, "mds1", "mds2", dA, [0.70, 0.05, 0.28, 0.34],
                           background=(emb.mds1.values, emb.mds2.values, bcol))
    keep = ~inA
    xs = list(o.mds1[keep]) + list(me.mds1)
    ys = list(o.mds2[keep]) + list(me.mds2)
    place_labels(ax, xs, ys, list(o.label[keep]) + med_lbl, avoid_xy=np.c_[list(o.mds1) + list(me.mds1),
                 list(o.mds2) + list(me.mds2)], fontsize=7.0, start=16, avoid_axes=[insA])
    n = nums["outliers_by_group"]
    ax.set_title(f"A  Shared MDS map: {nums['n_outliers']} consensus outliers - {n['cl7']} in cl7, {n['cl4']} in cl4, "
                 f"{n['minor']} in minor clusters", fontsize=10)
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc=BLUE, mec="none", alpha=0.6,
                label=f"cl7 chain (n={nums['chains_by_group']['cl7']})"),
         Line2D([], [], ls="", marker="o", ms=5, mfc=ORANGE, mec="none", alpha=0.6,
                label=f"cl4 chain (n={nums['chains_by_group']['cl4']})"),
         Line2D([], [], ls="", marker="^", ms=5, mfc=AQUA, mec="none", alpha=0.8,
                label=f"minor-cluster chain (n={nums['chains_by_group']['minor']})"),
         Line2D([], [], ls="", marker="*", ms=11, mfc=GREY, mec=INK, label="complete-H12 medoid (Y10)")]
    rq = (o.reason == "quality").sum()
    h += [Line2D([], [], ls="", marker="o", ms=7, mfc=REASON_COLOR["quality"], mec=INK,
                 label=f"outlier, Y16 reason quality ({rq})"),
          Line2D([], [], ls="", marker="o", ms=7, mfc=REASON_COLOR["conformation"], mec=INK,
                 label=f"outlier, Y16 reason conformation ({len(o) - rq})")]
    h += [Line2D([], [], ls="", marker=GROUP_MARKER[g], ms=6, mfc="white", mec=INK, label=f"outlier marker: {g}")
          for g in ("cl7", "cl4", "minor")]
    ax.legend(handles=h, loc="lower left", fontsize=7, handletextpad=0.3, borderaxespad=0.2)

    # ---------------- B: own vs other medoid RMSD
    lim = float(np.ceil(max(o.rmsd_own_medoid.max(), o.rmsd_other_medoid.max()) + 0.5))
    bx.plot([0, lim], [0, lim], color=MUTED, lw=1, ls="--", zorder=1)
    bx.text(lim * 0.97, lim * 0.97, "y = x (equidistant)", ha="right", va="bottom", rotation=45, fontsize=7.5,
            color=INK2, rotation_mode="anchor", transform_rotates_text=True)
    outlier_scatter(bx, o, "rmsd_own_medoid", "rmsd_other_medoid")
    bx.set_xlim(0, lim)
    bx.set_ylim(0, lim)
    bx.set_aspect("equal", adjustable="box")
    bx.set_xlabel("RMSD to own main-cluster medoid (A)")
    bx.set_ylabel("RMSD to the other main-cluster medoid (A)")
    dB = dense_mask(o.rmsd_own_medoid, o.rmsd_other_medoid)
    insB, inB = zoom_inset(bx, o, "rmsd_own_medoid", "rmsd_other_medoid", dB, [0.47, 0.04, 0.51, 0.38])
    keep = ~inB
    place_labels(bx, list(o.rmsd_own_medoid[keep]), list(o.rmsd_other_medoid[keep]), list(o.label[keep]),
                 avoid_xy=np.c_[o.rmsd_own_medoid, o.rmsd_other_medoid], fontsize=7.0, start=12,
                 avoid_axes=[insB])
    nc, nm = nums["n_main_closer_to_other"], nums["n_main"]
    head = (f"all {nm} main-cluster outliers lie closer to their own medoid" if nc == 0 else
            f"{nc} of {nm} main-cluster outliers lie closer to the other medoid")
    bx.set_title(f"B  Own vs other main-cluster medoid: {head}", fontsize=10)
    bx.text(0.02, 0.98,
            f"medoids (complete H12, Y10): cl7 {med_pdb[0]}, cl4 {med_pdb[1]}, {nums['rmsd_medoid_cl7_cl4']:.2f} A "
            f"apart\nminor-cluster chains: own = the nearer main medoid (Y15 home cluster)\n"
            f"above y = x: closer to own medoid; on y = x: equidistant from both",
            transform=bx.transAxes, ha="left", va="top", fontsize=7, color=INK2)
    fig.subplots_adjust(bottom=0.1)
    fig.text(0.01, 0.015, "Outliers: Y15 consensus (>= 3 of 5 detectors), reasons from Y16 (results/ext/"
             "outlier_verdicts.csv). Labels: PDB|author chain. Source: results/sv/SV06/outliers_on_map.csv.",
             fontsize=7.5, color=INK2)
    return fig, {"A_mds_map": ax, "B_own_vs_other_medoid": bx}


def heatmap_figure(o, nums, plt):
    """Y07 heatmap layout (scripts/ext_dendrogram.py) plus a 'Y16 outlier' side bar."""
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Patch
    d = np.load(ROOT / "results/rmsd_entry.npy")
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    leaf = pd.read_csv(ROOT / "results/ext/leaf_order.csv")
    assert (idx.chain_key.values[leaf.row.values] == leaf.chain_key.values).all()
    leaves = leaf.row.values
    t = y07.table()
    assert (t.chain_key.values == idx.chain_key.values).all()
    defs = y07.bars(t)
    sel = [defs[i] for i in (0, 1, 4, 5, 6)]        # the five bars of the Y07 heatmap
    reason = idx.chain_key.map(o.set_index("chain_key").reason).fillna("none")
    ocol = [to_rgb(REASON_COLOR.get(r, NO_OUTLIER)) for r in reason]
    sel.append(("Y16 outlier", ocol, {"quality": RED, "conformation": VIOLET}))
    is_out = reason.values[leaves] != "none"
    out_pos = np.where(is_out)[0]
    out_col = [REASON_COLOR[r] for r in reason.values[leaves][is_out]]

    fig = plt.figure(figsize=(13, 12.4))
    k = len(sel)
    gs = fig.add_gridspec(1 + k + 1, 1 + k + 1, height_ratios=[0.35] * k + [10, 0.4],
                          width_ratios=[0.35] * k + [10, 0.35], hspace=0.03, wspace=0.03)
    ax = fig.add_subplot(gs[k, k])
    im = ax.imshow(d[np.ix_(leaves, leaves)], cmap="esr_seq", vmin=0, vmax=6, interpolation="nearest", aspect="auto")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("788 chains in dendrogram leaf order (Y07)")
    axes = [ax]
    for i, (name, col, cmap) in enumerate(sel):
        a = fig.add_subplot(gs[i, k])
        b = fig.add_subplot(gs[k, i])
        if name == "Y16 outlier":
            a.set_facecolor(NO_OUTLIER); b.set_facecolor(NO_OUTLIER)
            a.set_xlim(-0.5, len(leaves) - 0.5); a.set_ylim(-0.5, 0.5)
            b.set_ylim(len(leaves) - 0.5, -0.5); b.set_xlim(-0.5, 0.5)
            for p, c in zip(out_pos, out_col):     # >= 1 column wide lines so single chains stay visible
                a.axvline(p, color=c, lw=1.6)
                b.axhline(p, color=c, lw=1.6)
            a.set_xticks([]); a.set_yticks([0], [f"{name} ({len(out_pos)})"], fontsize=7.5, fontweight="bold")
            b.set_yticks([]); b.set_xticks([0], [f"{name} ({len(out_pos)})"], rotation=90, fontsize=7.5,
                                           fontweight="bold")
        else:
            a.imshow(np.array(col)[leaves][None, :, :], aspect="auto", interpolation="nearest")
            b.imshow(np.array(col)[leaves][:, None, :], aspect="auto", interpolation="nearest")
            a.set_xticks([]); a.set_yticks([0], [name], fontsize=7.5)
            b.set_yticks([]); b.set_xticks([0], [name], rotation=90, fontsize=7.5)
        a.yaxis.tick_right()
        for s in list(a.spines.values()) + list(b.spines.values()):
            s.set_visible(False)
        axes += [a, b]
    cax = fig.add_subplot(gs[k + 1, k])
    fig.colorbar(im, cax=cax, orientation="horizontal", label="C-alpha RMSD (A), clipped at 6")
    axes.append(cax)
    hs = []
    for name, col, cmap in sel:
        if cmap:
            hs += [Patch(color=v, label=f"{name}: {kk}") for kk, v in list(cmap.items())[:5]]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.04, 0.965), fontsize=6.5, ncol=1)
    n = nums["outliers_by_group"]
    fig.suptitle(f"RMSD matrix in Y07 leaf order with the {len(out_pos)} Y16 consensus outliers marked: "
                 f"{n['cl4']} in the cl4 block, {n['cl7']} in the cl7 block, {n['minor']} in minor clusters",
                 x=0.01, ha="left", fontsize=11, fontweight="bold", y=0.995)
    return fig, {"heatmap_with_bars": axes}


def main():
    o = pd.read_csv(RES / "outliers_on_map.csv")
    nums = json.loads((RES / "numbers.json").read_text())
    plt = mpl_setup()
    fig, panels = map_figure(o, nums, plt)
    out = save_figure(fig, STEP, "SV06_outliers_map", panels)
    fig, panels = heatmap_figure(o, nums, plt)
    out += save_figure(fig, STEP, "SV06_rmsd_heatmap_outlier_bar", panels)
    for p in out:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
