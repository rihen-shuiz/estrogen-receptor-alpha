"""Are steroidal ligands antagonists? Action class of steroidal against non-steroidal ligands.

Data: results/12_curated_action/tables/steroid_summary.csv (counts per class at the level of unique ligands and of
PDB entry x ligand pairs; Fisher exact test of antagonist vs other classes, steroidal vs non-steroidal).
"""
import style as s

t = s.csv(f"{s.CURATED}/steroid_summary.csv")
rows = [("pairs_core3", "True", "steroidal"), ("pairs_core3", "False", "non-steroidal"),
        ("unique_ligands", "True", "steroidal"), ("unique_ligands", "False", "non-steroidal")]
fig, ax = s.plt.subplots(figsize=(9.0, 3.9))
ylab = []
for y, (lvl, st, lab) in enumerate(rows):
    r = t[(t.level == lvl) & (t.steroid.astype(str) == st)].iloc[0]
    left = 0
    for cls in s.CLASSES:
        k = int(r[cls])
        ax.barh(y, 100 * k / r.n, left=left, height=0.6, color=s.CLASS_COLOR[cls], edgecolor="white", linewidth=1.5)
        if k / r.n > 0.07:
            ax.text(left + 50 * k / r.n, y, str(k), ha="center", va="center", fontsize=10.5, color="white",
                    fontweight="bold")
        left += 100 * k / r.n
    ylab.append(f"{lab} {'pairs' if lvl == 'pairs_core3' else 'ligands'} ({int(r.n)})")
ax.set_yticks(range(len(rows)), ylab)
ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlim(0, 100)
ax.set_xlabel("Share of the set (%)")
for cls in s.CLASSES:
    ax.scatter([], [], marker="s", s=80, color=s.CLASS_COLOR[cls], label=s.CLASS_LABEL[cls])
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=10)
p = t[(t.level == "pairs_core3") & (t.steroid.astype(str) == "True")].iloc[0]
o = t[(t.level == "pairs_core3") & (t.steroid.astype(str) == "False")].iloc[0]
s.title(fig, f"Steroidal ligands are mostly agonists, not antagonists: {int(p.antagonist)} of {int(p.n)} steroid pairs "
             f"are antagonists against {int(o.antagonist)} of {int(o.n)} pairs of other ligands "
             f"(Fisher p = {p.antagonist_fisher_p:.1g}).", width=85)
fig.subplots_adjust(top=0.70)
s.save(fig, "curated_steroid_class")
