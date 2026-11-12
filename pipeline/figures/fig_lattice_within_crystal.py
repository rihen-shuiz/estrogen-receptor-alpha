"""Inside one crystal (same drop, same additives, same solvent), does the copy with more lattice contacts at the lid
have the more ordered or the less mobile lid?

Data: results/08_crystal_environment/lattice_within_crystal.csv (sign counts per entry; contacts = symmetry-mate atoms on 528-550).
"""
import numpy as np

import style as s

w = s.csv("results/08_crystal_environment/lattice_within_crystal.csv")
w = w[w.contact == "sym_atoms"]
rows = [("lid_order", "cl4", "lid order, cl4"), ("lid_order", "cl7", "lid order, cl7"),
        ("mobility_h12_z", "cl4", "lid mobility, cl4"), ("mobility_h12_z", "cl7", "lid mobility, cl7")]
fig, ax = s.plt.subplots(figsize=(9.5, 4.2))
for y, (o, g, lab) in enumerate(rows[::-1]):
    r = w[(w.outcome == o) & (w.group == g)].iloc[0]
    tot = r.more_contacts_better + r.more_contacts_worse + r.ties
    segs = [(r.more_contacts_better, s.GROUP_COLOR[g], 1.0), (r.ties, s.BAND, 1.0),
            (r.more_contacts_worse, s.GROUP_COLOR[g], 0.35)]
    left = 0
    for v, c, a in segs:
        ax.barh(y, v / tot, left=left, height=0.55, color=c, alpha=a, edgecolor="white", linewidth=2)
        left += v / tot
    ax.text(1.01, y, f"{int(r.more_contacts_better)} better / {int(r.more_contacts_worse)} worse / {int(r.ties)} equal"
                     f"   q {r.q:.2g}", va="center", fontsize=8.5, color=s.INK)
ax.set_yticks(range(len(rows)), [r[2] for r in rows[::-1]])
ax.set_xlim(0, 1)
ax.set_xticks(np.linspace(0, 1, 6), [f"{v:.0%}" for v in np.linspace(0, 1, 6)])
ax.set_xlabel("Share of crystals: copy with more lattice contacts at the lid has the better (dark), equal (light) "
              "or worse (pale) lid")
r = w[(w.outcome == "lid_order") & (w.group == "cl4")].iloc[0]
s.title(fig, f"Within one crystal, the cl4 copy with more lattice contacts at the lid has more of the lid modelled "
             f"in {int(r.more_contacts_better)} crystals and less in {int(r.more_contacts_worse)}; mobility shows "
             f"no clear rule.", width=100)
s.save(fig, "lattice_within_crystal")
