"""How many clusters does each cut of the tree give? Average-linkage dendrogram of the 788 chains.

Data: results/02_lid_groups/rmsd.npy (C-alpha RMSD after the core fit), results/02_lid_groups/groups.csv (cluster
numbers at 1.0, 1.5, 2.0 and 2.5 A).
"""
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

import style as s

d = np.load(s._path("results/02_lid_groups/rmsd.npy"))
g = s.csv("results/02_lid_groups/groups.csv")
z = linkage(squareform(d, checks=False), "average")
fig, ax = s.plt.subplots(figsize=(10, 4.8))
dendrogram(z, ax=ax, no_labels=True, color_threshold=0, above_threshold_color=s.INK2, truncate_mode="lastp", p=60,
           show_contracted=True)
for cut in (1.0, 1.5, 2.0, 2.5):
    k = g[f"cut_{cut}"].nunique()
    ax.axhline(cut, color=s.CL4 if cut == 1.5 else s.MUTED, lw=1.2, ls="-" if cut == 1.5 else "--")
    ax.text(ax.get_xlim()[1], cut, f"  {cut} Å: {k} clusters", va="center", fontsize=10, color=s.INK)
ax.set_ylabel("Average-linkage distance (C-alpha RMSD, Å)")
ax.set_xlabel("788 chains (tree cut to its last 60 branches)")
ax.spines["bottom"].set_visible(False)
s.grid(ax, "y")
top3 = g["cut_1.5"].value_counts().head(3).sum() / len(g)
s.title(fig, f"cl7 and cl4 join far above any cut we use; the 1.5 Å cut used throughout gives "
             f"{g['cut_1.5'].nunique()} clusters, and the three largest hold {top3:.0%} of the chains.", width=90)
s.save(fig, "dendrogram")
