"""Where along the chain do the two groups differ?

Data: results/02_lid_groups/displacement_profile.csv (distance between the cl7 and cl4 mean C-alpha positions),
results/02_lid_groups/numbers.json.
"""
import style as s

p = s.csv("results/02_lid_groups/displacement_profile.csv")
p = p[(p.n_cl7 >= 5) & (p.n_cl4 >= 5)]
n = s.numbers("02_lid_groups")
fig, ax = s.plt.subplots(figsize=(11, 4.4))
ax.axvspan(527.5, 550.5, color=s.BAND, lw=0)
top = float(p.distance_cl7_cl4_A.max()) * 1.08
ax.text(526.5, top * 0.97, "528-550: H11-H12 loop, H12, tail", ha="right", va="top", fontsize=8, color=s.INK2)
ax.plot(p.position, p.distance_cl7_cl4_A, color=s.INK, lw=1.6)
ax.set_xlim(305, 550)
ax.set_ylim(0, top)
ax.set_xlabel("ESR1 residue")
ax.set_ylabel("Distance between group-mean C-alpha (Å)")
s.grid(ax, "y")
s.title(fig, f"The two groups differ at the lid: residues 528-550 carry {n['share_sq_difference_528_550']:.1%} of the "
             f"squared difference; elsewhere the largest shift is {n['max_distance_outside_528_550_A']:.1f} Å, in loops.",
        width=110)
s.save(fig, "where_groups_differ")
