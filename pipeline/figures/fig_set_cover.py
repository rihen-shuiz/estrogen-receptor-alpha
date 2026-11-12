"""How many structures are needed to stand for all 788 chains?

Data: results/10_redundancy/set_cover.csv (representatives with a complete H12 needed for 90 % coverage per radius),
results/10_redundancy/cover_curve_1A.csv.
"""
import style as s

sc = s.csv("results/10_redundancy/set_cover.csv")
sc = sc[sc.coverage_target == 0.9]
fig, ax = s.plt.subplots(figsize=(7.5, 4.6))
ax.plot(sc.radius_A, sc.representatives, color=s.INK, lw=2, marker="o", ms=6)
for _, r in sc.iterrows():
    ax.annotate(f"{int(r.representatives)}", (r.radius_A, r.representatives), (0, 8), textcoords="offset points",
                ha="center", fontsize=9, color=s.INK)
ax.set_yscale("log")
ax.set_yticks([2, 3, 5, 10, 20, 50], ["2", "3", "5", "10", "20", "50"])
ax.minorticks_off()
ax.set_xlabel("Coverage radius (Å): a chain counts as covered within this RMSD of a representative")
ax.set_ylabel("Representatives for 90 % of chains (log)")
s.grid(ax, "y")
r1 = int(sc[sc.radius_A == 1.0].representatives.iloc[0])
r05 = int(sc[sc.radius_A == 0.5].representatives.iloc[0])
s.title(fig, f"{r1} structures with a complete H12 stand for 90 % of all 788 chains within 1.0 Å; at 0.5 Å it takes "
             f"{r05}: most of the PDB set repeats two conformations.", width=80)
s.save(fig, "set_cover")
