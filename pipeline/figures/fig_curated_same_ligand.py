"""Can one ligand hold the lid in both positions? Ligands crystallised in cl7 and in cl4.

Data: results/12_curated_action/tables/pairs_by_group.csv (one row per PDB entry x ligand x group: median lid-to-H11
distance, peptide bound, construct; ligands matched by full InChIKey).
"""
import numpy as np

import style as s

pg = s.csv(f"{s.CURATED}/pairs_by_group.csv")
pg = pg[pg.action_class.isin(s.CLASSES)].dropna(subset=["dist_H12c_H11axis_A"])
both = pg.groupby("inchikey").group.agg(lambda g: {"cl7", "cl4"} <= set(g))
d = pg[pg.inchikey.isin(both[both].index)]
order = d.groupby("lig_norm").size().sort_values().index.tolist()
rng = np.random.default_rng(42)
fig, ax = s.plt.subplots(figsize=(9.0, 4.2))
for y, lig in enumerate(order):
    e = d[d.lig_norm == lig]
    for _, r in e.iterrows():
        pep = str(r.peptide_bound) == "True"
        ax.scatter(r.dist_H12c_H11axis_A, y + rng.uniform(-0.15, 0.15), s=60, marker="o",
                   color=s.GROUP_COLOR.get(r.group, s.MINOR) if pep else "white",
                   edgecolor=s.GROUP_COLOR.get(r.group, s.MINOR), lw=1.6, zorder=3)
ax.set_yticks(range(len(order)), [f"{lig} ({s.CLASS_LABEL[d[d.lig_norm == lig].action_class.mode()[0]]})"
                                  for lig in order])
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlabel("Lid (H12) centre to H11 axis (Å), one point per PDB entry")
ax.scatter([], [], s=60, color=s.MUTED, label="coactivator peptide bound")
ax.scatter([], [], s=60, facecolor="white", edgecolor=s.MUTED, lw=1.6, label="no peptide")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=10)
ax.text(12.4, len(order) - 0.45, "cl7", ha="center", fontsize=11, color=s.CL7, fontweight="bold")
ax.text(22.3, len(order) - 0.45, "cl4", ha="center", fontsize=11, color=s.CL4, fontweight="bold")
ax.set_ylim(-0.6, len(order) - 0.1)
s.grid(ax, "x")
c7 = d[d.group == "cl7"]
k = int(c7.peptide_bound.astype(str).eq("True").sum())
s.title(fig, f"{len(order)} ligands are crystallised in both lid positions; {k} of their {len(c7)} closed-lid (cl7) "
             "entries carry a coactivator peptide: the crystal context, not only the ligand, goes with the lid position.",
        width=85)
fig.subplots_adjust(top=0.72)
s.save(fig, "curated_same_ligand")
