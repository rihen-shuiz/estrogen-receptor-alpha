"""Do the 788 chains fall into groups? The RMSD matrix ordered by group (bins of up to 10 chains within a group).

Data: results/02_lid_groups/rmsd.npy, groups.csv, numbers.json.
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import style as s

d = np.load(s._path("results/02_lid_groups/rmsd.npy"))
g = s.csv("results/02_lid_groups/groups.csv")
n = s.numbers("02_lid_groups")
grp = g["cut_1.5"].map({7: "cl7", 4: "cl4"}).fillna("minor")
order_key = grp.map({"cl7": 0, "minor": 1, "cl4": 2})
o = np.lexsort((g["cut_1.0"].to_numpy(), order_key.to_numpy()))
# bins of up to 10 chains inside each group (drawn as vector cells): mean RMSD per pair of bins
bins = []
for k in ("cl7", "minor", "cl4"):
    members = [i for i in o if grp.iloc[i] == k]
    bins += [members[j:j + 10] for j in range(0, len(members), 10)]
edges = np.cumsum([0] + [len(b) for b in bins])
m = np.array([[d[np.ix_(a, b)].mean() for b in bins] for a in bins])
cmap = LinearSegmentedColormap.from_list("seq", ["#0d366b", "#256abf", "#6da7ec", "#cde2fb", "#f7f6f3"])
fig, ax = s.plt.subplots(figsize=(7.4, 6.6))
im = ax.pcolormesh(edges, edges, m, cmap=cmap, vmin=0, vmax=5, linewidth=0, rasterized=False)
ax.set_xlim(0, edges[-1])
ax.set_ylim(edges[-1], 0)
ax.set_aspect("equal")
bounds = np.cumsum([0] + [int((grp == k).sum()) for k in ("cl7", "minor", "cl4")])
for b0, b1, k in zip(bounds[:-1], bounds[1:], ("cl7", "minor", "cl4")):
    ax.plot([-12, -12], [b0, b1 - 1], color=s.GROUP_COLOR[k], lw=6, solid_capstyle="butt", clip_on=False)
    ax.text(-22, (b0 + b1) / 2, f"{k}\n{b1 - b0}", ha="right", va="center", fontsize=8.5, color=s.INK)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
cb = fig.colorbar(im, ax=ax, shrink=0.75)
cb.set_label("C-alpha RMSD (Å), mean over bins of up to 10 chains")
cb.outline.set_visible(False)
cb.solids.set_rasterized(False)   # vector colour bar: no embedded image data in the SVG
s.title(fig, f"The 788 chains form two tight blocks: median RMSD {n['median_rmsd_within_cl7']:.2f} Å within cl7, "
             f"{n['median_rmsd_within_cl4']:.2f} Å within cl4 and {n['median_rmsd_cl7_cl4']:.2f} Å between them.",
        width=80)
s.save(fig, "rmsd_blocks")
