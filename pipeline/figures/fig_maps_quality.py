"""Which 2-D map keeps the RMSD distances and the nearest neighbours best?

Data: results/05_maps/maps_quality.csv (three new maps and the five October maps).
"""
import style as s

q = s.csv("results/05_maps/maps_quality.csv")
names = {"sammon": "Sammon", "som": "self-organising map", "umap_on_pca": "UMAP on PCA scores", "pca": "PCA",
         "cmds": "classical MDS", "mds": "metric MDS", "tsne": "t-SNE", "umap": "UMAP on RMSD"}
new = {"sammon", "som", "umap_on_pca"}
fig, ax = s.plt.subplots(figsize=(7.5, 5.6))
for _, r in q.iterrows():
    is_new = r["map"] in new
    ax.scatter(r.spearman_dist, r.trustworthiness, s=70 if is_new else 46, color=s.INK if is_new else "white",
               edgecolor=s.INK if is_new else s.MUTED, linewidth=1.3, zorder=3)
    ax.annotate(names[r["map"]], (r.spearman_dist, r.trustworthiness), (7, 4), textcoords="offset points",
                fontsize=8.5, color=s.INK if is_new else s.INK2)
ax.scatter([], [], s=70, color=s.INK, label="new maps (11 Nov)")
ax.scatter([], [], s=46, color="white", edgecolor=s.MUTED, linewidth=1.3, label="maps made in October")
ax.legend(loc="lower left")
ax.set_xlabel("Distances kept: Spearman correlation of map distance with RMSD")
ax.set_ylabel("Neighbours kept: trustworthiness (10 nearest)")
s.grid(ax, "both")
ax.set_xlim(0.55, 0.96)
ax.set_ylim(0.76, 0.95)
best = q.sort_values("spearman_dist", ascending=False).iloc[0]
s.title(fig, f"The Sammon map keeps distances best (Spearman {q.set_index('map').loc['sammon', 'spearman_dist']:.2f}"
             f", like metric MDS); UMAP and t-SNE keep neighbours but bend distances.", width=80)
s.save(fig, "maps_quality")
