"""With H12 hidden, where do the two groups that remain differ, and is that more than copies of one crystal differ?

Data: results/06_lid_hidden/no_lid_residue_profile.csv (split spectral_k2: distance between the group-mean C-alpha positions,
within-group RMS spread and copy-to-copy 95th percentile per residue), results/06_lid_hidden/numbers.json.
"""
import style as s

p = s.csv("results/06_lid_hidden/no_lid_residue_profile.csv")
p = p[(p.split == "spectral_k2") & (p.n_obs_group1 >= 5) & (p.n_obs_group2 >= 5)].sort_values("position")
n = s.numbers("06_lid_hidden")["splits"]["spectral_k2"]
lid = p.in_lid_region.astype(str) == "True"
q = p[~lid]
fig, ax = s.plt.subplots(figsize=(11, 4.8))
ax.axvspan(527.5, 550.5, color=s.BAND, lw=0)
ax.text(539, ax.get_ylim()[1] if False else 9.4, "hidden\n(528-550)", ha="center", va="top", fontsize=8,
        color=s.INK2)
ax.fill_between(q.position, 0, q.copy_p95_A, color=s.MUTED, alpha=0.18, lw=0, step="mid",
                label="copies of one crystal differ by up to (95th pct)")
ax.plot(q.position, q.within_rms_A, color=s.MUTED, lw=1, label="spread inside the groups (RMS)")
ax.plot(q.position, q.dist_A, color=s.INK, lw=1.6, label="distance between the two group means")
pk = q.loc[q.dist_A.idxmax()]
ax.annotate(f"loop 331-340: {pk.dist_A:.1f} Å at {int(pk.position)}", (pk.position, pk.dist_A),
            (pk.position + 25, pk.dist_A + 1.2), fontsize=8.5, color=s.INK,
            arrowprops=dict(arrowstyle="-", color=s.MUTED, lw=0.8))
ax.set_xlim(305, 550)
ax.set_ylim(0, 10)
ax.set_xlabel("ESR1 residue")
ax.set_ylabel("C-alpha distance (Å)")
s.grid(ax, "y")
ax.legend(loc="upper right", bbox_to_anchor=(0.93, 1.0))
share = n["loop_331_340_share"]
s.title(fig, f"With H12 hidden, the two groups that remain ({n['sizes'].replace(';', ' and ')} chains) differ almost "
             f"only at loop 331-340 ({share:.0%} of the squared difference), where copies of one crystal differ "
             f"as much.", width=110)
s.save(fig, "no_lid_profile")
