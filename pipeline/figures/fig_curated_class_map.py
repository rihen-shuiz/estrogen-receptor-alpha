"""Where do the action classes sit on the structure map (UMAP of the 788 chains)?

Data: results/05_maps/maps_coords.csv (UMAP of the first 6 PCA scores of the C-alpha coordinates, one point per
chain), results/12_curated_action/tables/chains_curated.csv (curated action class per chain; in_core3_strict = chain
whose single ligand has an agonist / partial agonist / antagonist label).
"""
import os

import style as s

SLIDE = bool(os.environ.get("ESR1_NO_TITLE"))  # wide layout, legend at the right, larger type

m = s.csv("results/05_maps/maps_coords.csv")
c = s.csv(f"{s.CURATED}/chains_curated.csv")[["chain_key", "action_class", "in_core3_strict"]]
m = m.merge(c, on="chain_key", how="left", validate="one_to_one")
lab = m.in_core3_strict.astype(str).eq("True")
fig, ax = s.plt.subplots(figsize=(11.5, 5.6) if SLIDE else (8.6, 6.2))
ax.scatter(m.umap1[~lab], m.umap2[~lab], s=10, color="#d9d8d3", lw=0, label=f"no class ({(~lab).sum()} chains)")
for cls in s.CLASSES:
    d = m[lab & (m.action_class == cls)]
    ax.scatter(d.umap1, d.umap2, s=48 if SLIDE else 34, marker=s.CLASS_MARKER[cls], facecolor="none",
               edgecolor=s.CLASS_COLOR[cls], lw=1.1, label=f"{s.CLASS_LABEL[cls]} ({len(d)} chains)")
for g, txt in (("cl7", "cl7: lid over the pocket"), ("cl4", "cl4: lid in the coactivator groove")):
    d = m[m.group == g]
    ax.text(d.umap1.median(), d.umap2.max() + 0.6, txt, ha="center", va="bottom", fontsize=14 if SLIDE else 11,
            color=s.INK, fontweight="bold")
ax.set_xticks([])
ax.set_yticks([])
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(s.GRID)
ax.set_xlabel("UMAP 1")
ax.set_ylabel("UMAP 2")
ax.margins(0.08)
if SLIDE:
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=13, markerscale=1.8)
else:
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=9.5, markerscale=1.3)
n7 = {cls: (lab & (m.action_class == cls) & (m.group == "cl7")).sum() for cls in s.CLASSES}
n4 = {cls: (lab & (m.action_class == cls) & (m.group == "cl4")).sum() for cls in s.CLASSES}
s.title(fig, f"On the structure map the classes follow the lid: agonist chains sit with cl7 ({n7['agonist']} vs "
             f"{n4['agonist']} in cl4), antagonist chains mostly with cl4 ({n4['antagonist']} vs {n7['antagonist']}); "
             "partial agonists fall in both.", width=80)
s.save(fig, "curated_class_map")
