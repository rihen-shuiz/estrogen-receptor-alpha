"""Do a Sammon map, a self-organising map and UMAP on PCA scores show anything beyond the two lid groups?

Data: results/05_maps/maps_coords.csv, results/03_lid_outliers/lid_outliers.csv, results/05_maps/numbers.json.
"""
import numpy as np

import style as s

m = s.csv("results/05_maps/maps_coords.csv")
fl = s.csv("results/03_lid_outliers/lid_outliers.csv")
n = s.numbers("05_maps")
rng = np.random.default_rng(42)
m["som_x"] = m.som_col + rng.uniform(-0.32, 0.32, len(m))
m["som_y"] = m.som_row + rng.uniform(-0.32, 0.32, len(m))
panels = [("sammon1", "sammon2", f"Sammon map of the RMSD matrix (stress {n['sammon_stress']:.3f})"),
          ("som_x", "som_y", f"Self-organising map, {n['som_grid'][0]} × {n['som_grid'][1]} units (points jittered)"),
          ("umap1", "umap2", f"UMAP of the first {n['pca_components_95pct']} PCA scores")]
fig, axs = s.plt.subplots(1, 3, figsize=(15, 5.4))
flagged = m.chain_key.isin(fl.chain_key)
for ax, (x, y, ttl) in zip(axs, panels):
    for grp in ("cl7", "cl4", "minor"):
        d = m[(m.group == grp) & ~flagged]
        ax.scatter(d[x], d[y], s=12 if grp != "minor" else 30, color=s.GROUP_COLOR[grp], alpha=0.55, lw=0,
                   marker=s.GROUP_MARKER[grp], label=s.GROUP_LABEL[grp])
    d = m[flagged]
    ax.scatter(d[x], d[y], s=40, facecolor=[s.GROUP_COLOR[g] for g in d.group], edgecolor=s.INK, linewidth=1.0,
               marker="D", zorder=3, label="lid outliers and grey chains (examined)")
    ax.set_title(ttl, fontsize=9.5, fontweight="normal", color=s.INK2)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(s.GRID)
axs[0].legend(loc="upper left", bbox_to_anchor=(0, -0.02), ncol=4, fontsize=8.5, markerscale=1.3)
q = n["silhouette_cl7_cl4"]
s.title(fig, "All three new maps separate the same two lid groups (silhouette cl7 vs cl4: Sammon "
             f"{q['sammon']:.2f}, UMAP {q['umap_on_pca']:.2f}, self-organising map {q['som']:.2f}); the examined "
             "chains sit at the edges or between the groups.", width=150)
s.save(fig, "maps_three")
