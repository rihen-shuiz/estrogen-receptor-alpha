"""Do steroidal ligands form their own chemical cluster, and which action classes do they carry?

Data: results/12_curated_action/tables/chem_ligands_order.csv (240 ligands with a curated class; UMAP of ECFP4 bits
with Jaccard distance, seed 42; steroid = carbon four-ring 6-6-6-5 nucleus; majority class over the ligand's pairs).
"""
import style as s

lig = s.csv(f"{s.CURATED}/chem_ligands_order.csv")
st = lig.steroid.astype(str).eq("True")
fig, ax = s.plt.subplots(figsize=(8.6, 6.2))
for cls in s.CLASSES:
    d = lig[(lig.action_class_majority == cls) & ~st]
    ax.scatter(d.umap_x, d.umap_y, s=40, marker=s.CLASS_MARKER[cls], facecolor="none", edgecolor=s.CLASS_COLOR[cls],
               lw=1.1, label=f"{s.CLASS_LABEL[cls]}")
    d = lig[(lig.action_class_majority == cls) & st]
    ax.scatter(d.umap_x, d.umap_y, s=70, marker=s.CLASS_MARKER[cls], color=s.CLASS_COLOR[cls], edgecolor=s.INK,
               lw=1.0, zorder=3)
sx, sy = lig[st].umap_x, lig[st].umap_y
ax.annotate(f"steroids ({st.sum()} ligands, filled)", (sx.median(), sy.max()), (sx.median(), sy.max() + 1.2),
            ha="center", fontsize=11, fontweight="bold", color=s.INK, arrowprops=dict(arrowstyle="-", color=s.MUTED))
ax.set_xticks([])
ax.set_yticks([])
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(s.GRID)
ax.set_xlabel("UMAP 1 (ECFP4 fingerprints)")
ax.set_ylabel("UMAP 2")
ax.margins(0.1)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.16), ncol=3, fontsize=10, markerscale=1.3)
k = lig[st].action_class_majority.value_counts()
s.title(fig, f"The {st.sum()} steroidal ligands occupy one region of chemical space (one K-means cluster); "
             f"{k.get('agonist', 0)} are agonists, "
             f"{k.get('partial_agonist', 0)} partial agonists and only {k.get('antagonist', 0)} antagonists.",
        width=80)
s.save(fig, "curated_steroid_map")
