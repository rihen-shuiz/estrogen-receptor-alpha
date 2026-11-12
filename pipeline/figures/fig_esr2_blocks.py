"""Does ESR2 repeat the two ESR1 lid positions?

Data: results/11_esr2/joint_block_medians.csv (median RMSD between groups in one joint frame), results/11_esr2/numbers.json.
"""
import style as s

b = s.csv("results/11_esr2/joint_block_medians.csv")
n = s.numbers("11_esr2")
b["label"] = b.group_a + "  vs  " + b.group_b
b = b.sort_values("median_rmsd_A")
fig, ax = s.plt.subplots(figsize=(8, 4))
y = range(len(b))
cols = [s.INK if "ESR1 cl4" in r.label and "ESR1 cl7" in r.label else s.MUTED for _, r in b.iterrows()]
ax.hlines(list(y), 0, b.median_rmsd_A, color=s.GRID, lw=3)
ax.scatter(b.median_rmsd_A, list(y), s=60, color=cols, zorder=3)
for yy, v in zip(y, b.median_rmsd_A):
    ax.text(v + 0.06, yy, f"{v:.2f} Å", va="center", fontsize=8.5, color=s.INK)
ax.set_yticks(list(y), b.label)
ax.set_xlim(0, 4.2)
ax.set_xlabel("Median C-alpha RMSD between the groups (Å, joint superposition)")
s.grid(ax, "x")
m = b.set_index("label").median_rmsd_A
s.title(fig, f"ESR2 shows the same two lid positions: its peptide-bound main group ({n['peptide_share_main']:.0%} "
             f"with peptide) sits {m['ESR1 cl7  vs  ESR2 main']:.2f} Å from ESR1 cl7, closer than the two ESR1 groups "
             f"are to each other ({m['ESR1 cl7  vs  ESR1 cl4']:.2f} Å).", width=85)
s.save(fig, "esr2_blocks")
