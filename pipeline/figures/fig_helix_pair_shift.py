"""Without superposition: which helix centres move relative to each other between the two lid groups?

Data: results/07_geometry_free/geometry_free_features.csv (median cl4 minus cl7 of every helix-helix and helix-chain COM distance,
with the ordinary-variation threshold per feature).
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import style as s

ft = s.csv("results/07_geometry_free/geometry_free_features.csv")
ft = ft[ft.kind == "distance"]
parts = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12", "COM"]
lab = {**{p: p for p in parts}, "COM": "chain centre"}
n = len(parts)
m = np.full((n, n), np.nan)
real = np.zeros((n, n), bool)
for _, r in ft.iterrows():
    i, j = parts.index(r.part_1), parts.index(r.part_2)
    i, j = max(i, j), min(i, j)
    m[i, j] = r.diff_cl4_minus_cl7
    real[i, j] = str(r.real) == "True"
cmap = LinearSegmentedColormap.from_list("div", ["#184f95", "#3987e5", "#9ec5f4", "#f0efec", "#f2a9a8", "#e34948",
                                                 "#a32a2a"])
lim = np.nanmax(np.abs(m))
fig, ax = s.plt.subplots(figsize=(7.6, 6.6))
im = ax.pcolormesh(np.arange(n + 1) - 0.5, np.arange(n + 1) - 0.5, np.ma.masked_invalid(m), cmap=cmap, vmin=-lim,
                   vmax=lim, linewidth=0)
ax.set_ylim(n - 0.5, -0.5)
ax.set_aspect("equal")
for i in range(n):
    for j in range(i):
        v = m[i, j]
        if np.isfinite(v):
            txt = f"{v:+.1f}" if abs(v) >= 0.95 else f"{v:+.2f}"
            ax.text(j, i, txt + ("\n●" if real[i, j] else ""), ha="center", va="center", fontsize=7.5,
                    color="white" if abs(v) > 0.55 * lim else s.INK)
ax.set_xticks(range(n), [lab[p] for p in parts], rotation=45, ha="right")
ax.set_yticks(range(n), [lab[p] for p in parts])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.tick_params(length=0)
cb = fig.colorbar(im, ax=ax, shrink=0.75)
cb.set_label("Distance in cl4 minus cl7 (Å, medians)")
cb.outline.set_visible(False)
cb.solids.set_rasterized(False)   # vector colour bar: no embedded image data in the SVG
nr = int(real.sum())
other = ft[~ft.feature.str.contains("H12")].diff_cl4_minus_cl7.abs().max()
s.title(fig, f"Without superposition, only distances to H12 change between the lid groups beyond ordinary variation "
             f"(● = {nr} of {int(np.isfinite(m).sum())}); every other distance changes by {other:.2f} Å or less.",
        width=75)
s.save(fig, "helix_pair_shift")
