"""Does expressing affinity relative to E2 measured in the same assay make different assays agree better?

Data: data/ligands/assay_spread.csv (ligands with at least two assays that also measured E2: SD of the absolute
pChEMBL across assays vs SD of dpX = pX(ligand) - pX(E2) in the same assay).
"""
import style as s

a = s.csv("data/ligands/assay_spread.csv").sort_values("sd_px_abs")
fig, ax = s.plt.subplots(figsize=(7.5, 4.6))
for y, (_, r) in enumerate(a.iterrows()):
    c = s.GROUP_COLOR.get(r.ligand_group, s.MUTED)
    ax.plot([r.sd_px_abs, r.sd_dpx], [y, y], color=s.GRID, lw=3, solid_capstyle="round", zorder=1)
    ax.scatter(r.sd_px_abs, y, s=60, facecolor="white", edgecolor=c, lw=1.6, zorder=3)
    ax.scatter(r.sd_dpx, y, s=60, color=c, zorder=3)
    ax.text(max(r.sd_px_abs, r.sd_dpx) + 0.06, y, f"{int(r.n_paired_assays)} assays", va="center", fontsize=8,
            color=s.INK2)
ax.set_yticks(range(len(a)), [f"{r.ccd} ({r.action_class}, {r.ligand_group})" for _, r in a.iterrows()])
ax.scatter([], [], s=60, facecolor="white", edgecolor=s.MUTED, lw=1.6, label="absolute pChEMBL")
ax.scatter([], [], s=60, color=s.MUTED, label="relative to E2 in the same assay")
ax.legend(loc="lower right")
ax.set_xlabel("Spread across assays (SD, log units)")
s.grid(ax, "x")
n_better = int((a.sd_dpx < a.sd_px_abs).sum())
s.title(fig, f"Referencing each affinity to E2 in the same assay narrows the spread between assays for {n_better} of "
             f"{len(a)} ligands (median SD {a.sd_px_abs.median():.2f} → {a.sd_dpx.median():.2f}).", width=80)
s.save(fig, "assay_spread")
