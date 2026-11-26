"""Where does each action class put the lid (H12)?

Data: results/12_curated_action/tables/pairs_core3_with_chem.csv (291 PDB entry x ligand pairs with a curated
agonist / partial agonist / antagonist label; h12_group_excl = cl7, cl4, minor, or 'multiple' when the chains of one
pair fall in different groups).
"""
import style as s

p = s.csv(f"{s.CURATED}/pairs_core3_with_chem.csv")
groups = [("cl7", "lid over the pocket (cl7)", s.CL7), ("cl4", "lid in the groove (cl4)", s.CL4),
          ("minor", "minor groups", s.MINOR), ("multiple", "chains in two groups", s.INK2)]
fig, ax = s.plt.subplots(figsize=(9.2, 3.6))
for y, cls in enumerate(s.CLASSES):
    d = p[p.action_class == cls]
    left = 0
    for g, _, c in groups:
        k = int((d.h12_group_excl == g).sum())
        if k:
            ax.barh(y, 100 * k / len(d), left=left, height=0.6, color=c, edgecolor="white", linewidth=1.5)
            if k / len(d) > 0.06:
                ax.text(left + 50 * k / len(d), y, str(k), ha="center", va="center", fontsize=11, color="white",
                        fontweight="bold")
            left += 100 * k / len(d)
ax.set_yticks(range(3), [f"{s.CLASS_LABEL[c]}\n{(p.action_class == c).sum()} pairs" for c in s.CLASSES])
ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlim(0, 100)
ax.set_xlabel("Share of the class's PDB entry × ligand pairs (%)")
for g, lab, c in groups:
    ax.scatter([], [], marker="s", s=80, color=c, label=lab)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=9.5, handletextpad=0.3, columnspacing=1.2)
n = {c: (p.action_class == c).sum() for c in s.CLASSES}
k = {c: ((p.action_class == c) & (p.h12_group_excl == "cl7")).sum() for c in s.CLASSES}
k4 = ((p.action_class == "partial_agonist") & (p.h12_group_excl == "cl4")).sum()
s.title(fig, f"Agonists close the lid ({k['agonist']} of {n['agonist']} in cl7); partial agonists mostly do "
             f"({k['partial_agonist']} cl7, {k4} cl4); {k['antagonist']} of {n['antagonist']} antagonists still sit "
             "in the closed position.", width=82)
fig.subplots_adjust(top=0.72)
s.save(fig, "curated_class_h12")
