"""Does ESR2 share the two ESR1 conformations? Joint landscape of the ESR1 and ESR2 chains.

Data: results/11_esr2/joint_map.csv (classical MDS of the joint 874 x 874 C-alpha RMSD matrix, ESR1 numbering, same
reference frame), results/11_esr2/numbers.json (group sizes, peptide shares, median RMSD between groups).
"""
import style as s

j = s.csv("results/11_esr2/joint_map.csv")
n = s.numbers("11_esr2")
fig, ax = s.plt.subplots(figsize=(8.4, 5.8))
e1 = j[j.receptor == "ESR1"]
for cl, c, lab in (("ESR1 cl7", s.CL7, "ESR1 cl7 (lid over the pocket)"), ("ESR1 cl4", s.CL4, "ESR1 cl4 (lid in the groove)")):
    d = e1[e1.cluster == cl]
    ax.scatter(d.mds1, d.mds2, s=14, color=c, alpha=0.45, lw=0, label=f"{lab}, {len(d)}")
d = e1[~e1.cluster.isin(["ESR1 cl7", "ESR1 cl4"])]
ax.scatter(d.mds1, d.mds2, s=30, color=s.MINOR, marker="^", lw=0, label=f"ESR1 minor groups, {len(d)}")
e2 = j[j.receptor == "ESR2"]
sizes = e2.cluster.value_counts()
for cl, mk, face in ((sizes.index[0], "D", s.INK), (sizes.index[1], "D", "white")):
    d = e2[e2.cluster == cl]
    lab = "ESR2 main group" if cl == sizes.index[0] else "ESR2 second group"
    ax.scatter(d.mds1, d.mds2, s=46, marker=mk, facecolor=face, edgecolor=s.INK, lw=1.1, zorder=3, label=f"{lab}, {len(d)}")
d = e2[~e2.cluster.isin(sizes.index[:2])]
ax.scatter(d.mds1, d.mds2, s=46, marker="D", facecolor=s.MINOR, edgecolor=s.INK, lw=0.8, zorder=3,
           label=f"ESR2 small groups, {len(d)}")
ax.set_xlabel("MDS 1 (Å)")
ax.set_ylabel("MDS 2 (Å)")
ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=9.5)
s.grid(ax, "both")
bm = n["block_medians"]
s.title(fig, f"ESR2 has the same two states: its main group (peptide in {n['peptide_share_main']:.0%}) lies "
             f"{bm['ESR1 cl7 - ESR2 main']:.2f} Å from ESR1 cl7 and its second group (no peptide) "
             f"{bm['ESR1 cl4 - ESR2 second']:.2f} Å from cl4, against {bm['ESR1 cl7 - ESR1 cl4']:.2f} Å between cl7 and cl4.",
        width=85)
s.save(fig, "esr2_map")
