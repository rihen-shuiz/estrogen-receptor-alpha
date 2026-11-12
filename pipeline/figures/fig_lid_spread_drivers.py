"""Does any covariate move the lid within a group by more than ordinary variation?

Data: results/03_lid_outliers/lid_spread_drivers.csv (entry-level effects, q values, thresholds).
"""
import numpy as np

import style as s

dr = s.csv("results/03_lid_outliers/lid_spread_drivers.csv").dropna(subset=["effect"])
names = {"peptide_bound": "coactivator peptide bound", "L536S": "L536S construct", "Y537S": "Y537S construct",
         "basic_amine": "ligand has a basic amine", "agonist_vs_SERM": "agonist vs SERM (ChEMBL class)",
         "h12_complete": "H12 fully modelled", "heavy_atoms": "ligand size (heavy atoms)",
         "resolution_A": "resolution", "sym_atoms": "lattice contacts at the lid"}
dr["label"] = dr.covariate.map(lambda c: names.get(c, "largest 1.0 Å sub-group vs the rest"))
fig, axs = s.plt.subplots(1, 2, figsize=(10.5, 5), sharex=True)
for ax, grp in zip(axs, ("cl7", "cl4")):
    d = dr[dr.group == grp].sort_values("effect")
    thr = float(d.threshold.iloc[0])
    ax.axvspan(-thr, thr, color=s.GROUP_COLOR[grp], alpha=0.08, lw=0)
    ax.axvline(0, color=s.MUTED, lw=0.8)
    y = np.arange(len(d))
    sig = (d.q < 0.05).to_numpy()
    ax.scatter(d.effect[sig], y[sig], s=46, color=s.GROUP_COLOR[grp], zorder=3, label="q < 0.05")
    ax.scatter(d.effect[~sig], y[~sig], s=46, facecolor="white", edgecolor=s.GROUP_COLOR[grp], linewidth=1.4,
               zorder=3, label="q ≥ 0.05 or untested")
    ax.set_yticks(y, d.label)
    ax.set_title(f"{grp}: shaded = ± {thr:.2f} Å (copies of one crystal)", fontsize=9.5, fontweight="normal",
                 color=s.INK2)
    s.grid(ax, "x")
    ax.set_xlabel("Shift of the lid-to-H11 distance (Å, entry medians)")
axs[1].legend(loc="upper left")
top = dr.loc[dr.effect.abs().idxmax()]
s.title(fig, "Inside each lid group, no ligand, construct, peptide or crystal property shifts the lid by more than "
             f"copies of one crystal differ; the largest, {top.label} in {top.group}, is {abs(top.effect):.2f} Å.")
fig.tight_layout()
s.save(fig, "lid_spread_drivers")
