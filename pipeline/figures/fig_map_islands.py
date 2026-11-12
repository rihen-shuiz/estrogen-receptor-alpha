"""Are the islands that the maps break into separated by more than copies of one crystal?

Data: results/05_maps/maps_islands.csv (HDBSCAN islands per map, within/between median RMSD, copy-to-copy 95th pct).
"""
import numpy as np

import style as s

isl = s.csv("results/05_maps/maps_islands.csv")
names = {"sammon": "Sammon map", "som": "self-organising map", "umap_on_pca": "UMAP on PCA scores"}
copy = float(isl.copy_rmsd_p95.iloc[0])
fig, ax = s.plt.subplots(figsize=(9, 4.6))
rng = np.random.default_rng(42)
for y, mp in enumerate(["umap_on_pca", "som", "sammon"]):
    d = isl[isl["map"] == mp]
    for grp in ("cl7", "cl4"):
        e = d[d.main_group == grp]
        ax.scatter(e.gap_A, y + rng.uniform(-0.15, 0.15, len(e)), s=12 + e.n_chains * 0.6, alpha=0.75,
                   color=s.GROUP_COLOR[grp], lw=0, label=f"island mostly {grp}" if y == 0 else None)
    ax.text(-0.42, y + 0.32, f"{len(d)} islands", fontsize=8, color=s.INK2)
ax.axvline(copy, color=s.INK, lw=1.2)
ax.text(copy, 2.55, f"copies of one crystal differ by\nup to {copy:.2f} Å (95th pct)", fontsize=8, color=s.INK,
        ha="right", va="top")
ax.axvline(0, color=s.MUTED, lw=0.8)
ax.set_yticks([0, 1, 2], [names[k] for k in ["umap_on_pca", "som", "sammon"]])
ax.set_xlim(-0.45, 1.45)
ax.set_ylim(-0.5, 2.6)
ax.set_xlabel("Gap: median RMSD to the other islands of the same lid group minus median RMSD inside the island (Å)")
s.grid(ax, "x")
ax.legend(loc="lower right")
s.title(fig, f"No island on any map is set apart by more than copies of one crystal differ (largest gap "
             f"{isl.gap_A.max():.2f} Å vs {copy:.2f} Å); the extra islands follow crystal form and construct.")
s.save(fig, "map_islands")
