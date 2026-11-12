"""Do additives, the crystallisation recipe, solvent content or lattice contacts change lid order or mobility by more
than ordinary variation?

Data: results/08_crystal_environment/additives_lattice_tests.csv (entry-level effects within cl7 and cl4, thresholds, q values).
"""
import numpy as np

import style as s

t = s.csv("results/08_crystal_environment/additives_lattice_tests.csv").dropna(subset=["effect"])
t["ratio"] = t.effect / t.threshold
labels = {"additives_per_chain": "bound additives per chain", "other_het_per_chain": "other bound small molecules",
          "recipe_classes": "recipe component classes", "solvent_pct": "solvent content of the crystal",
          "sym_atoms": "lattice contacts at the lid", "n_residues_sym_contact": "lattice contacts, whole chain"}
order = list(labels)
fig, axs = s.plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, (o, ttl) in zip(axs, (("lid_order", "lid order (C-alpha of 538-547 modelled)"),
                              ("mobility_h12_z", "lid mobility (B z-score)"))):
    ax.axvspan(-1, 1, color=s.BAND, lw=0)
    ax.axvline(0, color=s.MUTED, lw=0.8)
    for grp, off in (("cl7", -0.12), ("cl4", 0.12)):
        d = t[(t.outcome == o) & (t.group == grp)].set_index("predictor").reindex(order).dropna(subset=["ratio"])
        y = np.array([order.index(p) for p in d.index]) + off
        sig = (d.q < 0.05).to_numpy()
        ax.scatter(d.ratio[sig], y[sig], s=46, color=s.GROUP_COLOR[grp], zorder=3, label=f"{grp}, q < 0.05")
        ax.scatter(d.ratio[~sig], y[~sig], s=46, facecolor="white", edgecolor=s.GROUP_COLOR[grp], lw=1.4, zorder=3,
                   label=f"{grp}, q ≥ 0.05")
    ax.set_yticks(range(len(order)), [labels[p] for p in order])
    ax.set_xlim(-1.3, 1.3)
    ax.set_title(ttl, fontsize=9.5, fontweight="normal", color=s.INK2)
    ax.set_xlabel("Effect ÷ ordinary variation (shaded ±1 = copies of one crystal)")
    s.grid(ax, "x")
axs[1].legend(loc="lower left", fontsize=8)
mx = t.loc[t.ratio.abs().idxmax()]
both = sorted({labels[p] for p in t[(t.group == "cl4") & (t.q < 0.05)].groupby("predictor").outcome.nunique()
               .loc[lambda x: x == 2].index})
s.title(fig, "Additives, recipe, solvent content and lattice contacts all stay inside ordinary variation of lid order "
             f"and mobility (largest: {labels[mx.predictor]}, {mx.group}, {abs(mx.ratio):.2f} of it); in cl4 only "
             f"{'the two lattice-contact measures' if all('lattice' in b for b in both) else ' and '.join(both)} "
             f"go with both outcomes (q < 0.05).", width=110)
fig.tight_layout()
s.save(fig, "crystal_components")
