"""Are the ligands of the two groups different in size and in basic amines?

Data: data/ligands/ligand_descriptors.csv (one row per ligand code; group = where >= 80 % of its chains fall).
"""
import numpy as np

import style as s

ld = s.csv("data/ligands/ligand_descriptors.csv")
ld = ld[ld.structural_majority.isin(["cl7", "cl4"])]
fig, ax = s.plt.subplots(figsize=(8.5, 4.6))
rng = np.random.default_rng(42)
stats = {}
for y, grp in enumerate(("cl4", "cl7")):
    d = ld[ld.structural_majority == grp]
    amine = d.basic_amine.astype(str).eq("True").to_numpy()
    yy = y + rng.uniform(-0.2, 0.2, len(d))
    ax.scatter(d.heavy_atoms[~amine], yy[~amine], s=18, color=s.GROUP_COLOR[grp], alpha=0.55, lw=0)
    ax.scatter(d.heavy_atoms[amine], yy[amine], s=26, facecolor="white", edgecolor=s.GROUP_COLOR[grp], lw=1.2)
    med = d.heavy_atoms.median()
    ax.plot([med, med], [y - 0.3, y + 0.3], color=s.INK, lw=2)
    stats[grp] = (len(d), med, amine.mean())
    ax.text(70, y, f"{len(d)} ligands\nmedian {med:.0f} heavy atoms\nbasic amine {amine.mean():.0%}", va="center",
            fontsize=8.5, color=s.INK)
ax.set_yticks([0, 1], ["cl4\n(lid in groove)", "cl7\n(lid over pocket)"])
ax.set_xlim(0, 85)
ax.set_xlabel("Ligand heavy atoms (open circle: ligand with a basic amine)")
s.grid(ax, "x")
s.title(fig, f"Ligands of cl4 chains are larger (median {stats['cl4'][1]:.0f} vs {stats['cl7'][1]:.0f} heavy atoms) "
             f"and usually carry a basic amine ({stats['cl4'][2]:.0%} vs {stats['cl7'][2]:.0%}).", width=90)
s.save(fig, "ligand_size")
