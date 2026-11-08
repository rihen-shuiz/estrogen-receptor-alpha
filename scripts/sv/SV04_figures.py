"""SV04 figure (step 6), drawn only from the saved results/sv/SV04/h12_map_coords.csv (rule 1.8).

h12_map_partial: the SV03 H12 map (x = distance from the H12 point to the H11 axis, y = angle H12-H3; colour = S08
cluster group, marker = ligand_chemotype; scripts/sv/common/h12_map.plot_h12_map) with the partial-H12 chains
(SV01a h12_status 'partial') added as hollow markers at their SV04 estimated positions (H12 midpoint = axis point of
residue 542.5; H12, H3 and H10/11 eligible).
  (a) as the step specifies: complete-H12 chains filled at their SV03 positions (H12 centroid = mean of the
      4-residue window centroids of the observed helical C-alpha), partial-H12 chains hollow at SV04 estimates;
  (b) the same chains with one estimator for all: complete-H12 chains filled and partial-H12 chains hollow, both at
      their SV04 midpoint estimates. Added because the H12 point of the two methods differs along the H12 axis
      (SV04 midpoint - SV03 centroid, numbers_midpoints.json), which shifts x by a cluster-dependent offset
      (numbers_cluster.json complete_sv04_minus_sv03_map); y is identical in both.
Each graph is also saved as its own cropped PNG in figures/panels/sv/SV04/.
Run: .venv/bin/python -B -I scripts/sv/SV04_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.h12_map import (CHEMOTYPE_LABEL, CHEMOTYPE_MARKER, GROUP_COLOR, GROUP_LABEL,  # noqa: E402
                            GROUP_ORDER, X_LABEL, Y_LABEL, add_legends, plot_h12_map)
from common.sv_fig import INK2, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV04"
STEP = "SV04"


def panel(ax, filled, hollow, fx, fy, title, note):
    from matplotlib.lines import Line2D
    hf = plot_h12_map(ax, filled, x=fx, y=fy, legend=False, label_axes=False)
    hh = plot_h12_map(ax, hollow, x="map_x_sv04", y="map_y_sv04", hollow=True, legend=False, label_axes=False,
                      size=30, zorder=3)
    extra = [Line2D([], [], ls="", marker="o", markersize=6, color=GROUP_COLOR[g], markerfacecolor="none",
                    markeredgewidth=1.0, label=f"{GROUP_LABEL[g]}, partial H12 (n={int((hollow.group == g).sum())})")
             for g in GROUP_ORDER if (hollow.group == g).any()]
    groups = [(f"{lab}, complete H12", c, n) for lab, c, n in hf["groups"]]
    cnt = pd.concat([filled.ligand_chemotype, hollow.ligand_chemotype]).value_counts()
    markers = [(CHEMOTYPE_LABEL.get(c, c), CHEMOTYPE_MARKER.get(c, "o"), int(cnt[c]))
               for c in CHEMOTYPE_MARKER if c in cnt.index]
    add_legends(ax, {"groups": groups, "markers": markers, "hollow": False}, extra=extra,
                anchor_markers=(1.02, 0.52))
    ax.set_xlabel(X_LABEL.replace("H12 centroid", "H12 point"))
    ax.set_ylabel(Y_LABEL)
    ax.set_title(title, fontsize=9.2)
    ax.text(0.01, 0.99, note, transform=ax.transAxes, ha="left", va="top", fontsize=6.6, color=INK2)
    ax.grid(alpha=0.25, lw=0.5)
    return hf["n"], hh["n"]


def main():
    plt = mpl_setup()
    m = pd.read_csv(RES / "h12_map_coords.csv", keep_default_na=False, na_values=[""])
    filled = m[m.map_layer.str.startswith("filled")]
    hollow = m[m.map_layer.str.startswith("hollow")]
    filled_sv04 = filled[filled.map_x_sv04.notna() & filled.map_y_sv04.notna()]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9.6, 11.4))
    fig.subplots_adjust(right=0.70, hspace=0.32)
    n1, n2 = panel(a1, filled, hollow, "map_x_sv03", "map_y_sv03",
                   f"(a) H12 map: {len(filled)} complete-H12 chains at SV03 positions (filled), "
                   f"{len(hollow)} partial-H12 chains at SV04 estimates (hollow)",
                   "filled: x from the H12 centroid of the observed helical Cα (SV03)\n"
                   "hollow: x from the estimated H12 midpoint (axis point of residue 542.5, SV04)\n"
                   "core v3 frame; H11 = SV02 consensus H10/11 (497-526)")
    panel(a2, filled_sv04, hollow, "map_x_sv04", "map_y_sv04",
          f"(b) same chains, all at SV04 estimates: {len(filled_sv04)} complete-H12 (filled), "
          f"{len(hollow)} partial-H12 (hollow)",
          "x from the estimated H12 midpoint (axis point of residue 542.5) for every chain\n"
          "y identical to (a); core v3 frame; H11 = SV02 consensus H10/11 (497-526)")
    xs = pd.concat([filled.map_x_sv03, filled_sv04.map_x_sv04, hollow.map_x_sv04])
    ys = pd.concat([filled.map_y_sv03, hollow.map_y_sv04])
    for ax in (a1, a2):
        ax.set_xlim(xs.min() - 1.0, xs.max() + 1.0)
        ax.set_ylim(ys.min() - 4.0, ys.max() + 4.0)
    out = save_figure(fig, STEP, "h12_map_partial", panels={"a_sv03_positions": a1, "b_sv04_all": a2})
    for p in out:
        print(p.relative_to(ROOT), f"{p.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main()
