"""Where do ESR1 and ESR2 differ residue by residue, and do the pocket selectivity positions move?

Data: data/esr2/esr1_esr2_displacement_profile.csv (median and 90th percentile C-alpha distance per ESR1 position over
135 matched ESR1/ESR2 chain pairs with the same ligand, after a core fit; residue identity of each pair of positions).
Selectivity positions: Pike et al. 1999, EMBO J 18:4608 (ESR1 L384 / ESR2 M336, ESR1 M421 / ESR2 I373).
"""
import numpy as np

import style as s

p = s.csv("data/esr2/esr1_esr2_displacement_profile.csv")
full = p[p.n_pairs == p.n_pairs.max()]
fig, ax = s.plt.subplots(figsize=(11, 4.4))
ax.fill_between(p.esr1_position, 0, p.p90_displacement_A.clip(upper=6), color=s.BAND, lw=0, step="mid")
ax.plot(p.esr1_position, p.median_displacement_A.clip(upper=6), color=s.INK, lw=1.4, drawstyle="steps-mid")
diff = p[~p.same_residue.astype(bool)]
ax.scatter(diff.esr1_position, np.full(len(diff), -0.25), marker="|", s=40, color=s.INK2, lw=1)
ax.text(p.esr1_position.min(), -0.45, "residue differs between ESR1 and ESR2", fontsize=9, color=s.INK2, va="top")
for pos in (384, 421):
    r = p[p.esr1_position == pos].iloc[0]
    ax.annotate(f"{r.esr1_aa}{pos} / {r.esr2_aa}{int(r.esr2_position)}\n{r.median_displacement_A:.2f} Å", (pos, r.median_displacement_A),
                (pos, 3.2), ha="center", fontsize=10, color=s.INK, arrowprops=dict(arrowstyle="-", color=s.MUTED))
ax.set_ylim(-0.7, 6.2)
ax.set_xlim(p.esr1_position.min() - 2, p.esr1_position.max() + 2)
ax.set_xlabel("ESR1 residue")
ax.set_ylabel("C-alpha distance ESR1-ESR2 (Å)\nmedian; shaded: 90th percentile")
s.grid(ax, "y")
s.title(fig, f"After a core fit ESR1 and ESR2 backbones agree to a median {full.median_displacement_A.median():.2f} Å "
             f"over {int(p.n_pairs.max())} matched pairs; the pocket selectivity positions move only 0.25-0.28 Å, so "
             "selectivity comes from the side chains (L/M, M/I), not from backbone shifts.", width=100)
s.save(fig, "esr2_displacement")
