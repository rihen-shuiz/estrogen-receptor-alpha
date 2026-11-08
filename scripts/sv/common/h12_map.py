"""H12 map (SV03 step 6), reused by SV04 (partial-H12 chains as hollow markers) and SV13 (marker = action class).

x = distance from the H12 centroid to the H11 axis (A), y = angle H12-H3 (deg); one point per chain with an H12
axis; colour = S08 cluster group (cl7 blue, cl4 orange, minor clusters aqua: the project palette, rule 1.8);
marker = a categorical column (default `ligand_chemotype` of the SV00 chain table).
Columns come from results/sv/SV03/geometry_per_chain.csv (dist_H12c_H11axis_A, angle_H12_H3_deg, group,
ligand_chemotype).

Usage:
    from common.h12_map import plot_h12_map, add_legends
    h = plot_h12_map(ax, df)                                   # SV03
    h2 = plot_h12_map(ax, partial, hollow=True, legend=False)  # SV04 overlay, then add_legends(ax, ...)
    plot_h12_map(ax, df, marker_col="action_class", marker_map=..., marker_labels=..., marker_title="action class")
"""
from __future__ import annotations

import pandas as pd

from .sv_fig import AQUA, BLUE, GREY, INK2, ORANGE

X_COL, Y_COL = "dist_H12c_H11axis_A", "angle_H12_H3_deg"
X_LABEL = "distance H12 centroid to H11 axis (Å)"
Y_LABEL = "angle H12–H3 (°)"
GROUP_ORDER = ("cl7", "cl4", "minor")
GROUP_COLOR = {"cl7": BLUE, "cl4": ORANGE, "minor": AQUA}
GROUP_LABEL = {"cl7": "cl7", "cl4": "cl4", "minor": "minor clusters"}
CHEMOTYPE_MARKER = {"other": "o", "amine": "^", "acid": "s", "both": "D", "none": "X"}
CHEMOTYPE_LABEL = {"other": "neither amine nor acid", "amine": "basic amine", "acid": "carboxylic acid",
                   "both": "amine + acid", "none": "no ligand of interest"}


def cluster_group(s08_cluster) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(s08_cluster), "minor")


def plot_h12_map(ax, df: pd.DataFrame, x: str = X_COL, y: str = Y_COL, group_col: str = "group",
                 marker_col: str = "ligand_chemotype", marker_map: dict | None = None,
                 marker_labels: dict | None = None, marker_title: str = "ligand chemotype (SV00)",
                 color_map: dict | None = None, group_labels: dict | None = None, hollow: bool = False,
                 size: float = 22, alpha: float = 0.75, lw: float = 0.5, zorder: float = 2, legend: bool = True,
                 legend_kw: dict | None = None, label_axes: bool = True, default_marker: str = "o",
                 edge_map: dict | None = None) -> dict:
    """Scatter the H12 map on `ax`. Rows with NaN in x or y are skipped.

    hollow=True draws open markers (edge = group colour). edge_map (optional, added by SV13): marker category ->
    edge colour of filled markers (default white), e.g. the SV12 class colours. Returns {'groups': [(label, color, n)],
    'markers': [(label, marker, n)], 'n': points drawn} so a caller can build its own legend (add_legends).
    """
    marker_map = marker_map or CHEMOTYPE_MARKER
    marker_labels = marker_labels or CHEMOTYPE_LABEL
    color_map = color_map or GROUP_COLOR
    group_labels = group_labels or GROUP_LABEL
    d = df[df[x].notna() & df[y].notna()]
    groups = [g for g in GROUP_ORDER if g in set(d[group_col])] + \
             sorted(set(d[group_col]) - set(GROUP_ORDER), key=str)
    cats = [c for c in marker_map if c in set(d[marker_col])] + \
           sorted(set(d[marker_col].dropna()) - set(marker_map), key=str)
    for g in groups:
        col = color_map.get(g, GREY)
        for c in cats:
            s = d[(d[group_col] == g) & (d[marker_col] == c)]
            if not len(s):
                continue
            m = marker_map.get(c, default_marker)
            kw = dict(marker=m, s=size, alpha=alpha, zorder=zorder + (0.2 if g == "minor" else 0))
            if hollow:
                kw.update(facecolors="none", edgecolors=col, linewidths=max(lw, 0.9))
            elif m in ("x", "X", "+", "1", "2", "3", "4"):
                kw.update(color=col, linewidths=lw)
            else:
                kw.update(color=col, edgecolors=(edge_map or {}).get(c, "white"), linewidths=lw)
            ax.scatter(s[x], s[y], **kw)
    if label_axes:
        ax.set_xlabel(X_LABEL if x == X_COL else x)
        ax.set_ylabel(Y_LABEL if y == Y_COL else y)
    out = {"groups": [(group_labels.get(g, str(g)), color_map.get(g, GREY), int((d[group_col] == g).sum()))
                      for g in groups],
           "markers": [(marker_labels.get(c, str(c)), marker_map.get(c, default_marker),
                        int((d[marker_col] == c).sum())) for c in cats],
           "n": int(len(d)), "hollow": hollow}
    if legend:
        add_legends(ax, out, marker_title=marker_title, **(legend_kw or {}))
    return out


def add_legends(ax, handles: dict, marker_title: str = "ligand chemotype (SV00)", extra: list | None = None,
                anchor_groups: tuple = (1.02, 1.0), anchor_markers: tuple = (1.02, 0.6), fontsize: float = 7):
    """Two legends outside the axes on the right: colour = cluster group, marker = category (grey).
    `extra` = additional Line2D handles appended to the cluster legend (e.g. SV04 hollow = partial H12)."""
    from matplotlib.lines import Line2D
    hollow = handles.get("hollow", False)
    gh = [Line2D([], [], ls="", marker="o", markersize=6, color=c, markerfacecolor="none" if hollow else c,
                 label=f"{lab} (n={n})") for lab, c, n in handles["groups"]]
    leg1 = ax.legend(handles=gh + (extra or []), title="S08 cluster", loc="upper left", bbox_to_anchor=anchor_groups,
                     fontsize=fontsize, title_fontsize=fontsize, frameon=False, alignment="left")
    ax.add_artist(leg1)
    mh = [Line2D([], [], ls="", marker=m, markersize=6, color=INK2, markerfacecolor="none" if hollow else INK2,
                 label=f"{lab} (n={n})") for lab, m, n in handles["markers"]]
    ax.legend(handles=mh, title=marker_title, loc="upper left", bbox_to_anchor=anchor_markers, fontsize=fontsize,
              title_fontsize=fontsize, frameon=False, alignment="left")
