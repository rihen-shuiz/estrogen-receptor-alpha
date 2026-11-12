"""Is a longer H12 (phi/psi) a less mobile H12?

Data: results/07_geometry_free/geometry_free_per_chain.csv (len_H12 from phi/psi, H12 B z-score), entry means within each group.
"""
import numpy as np
from scipy.stats import spearmanr

import style as s

per = s.csv("results/07_geometry_free/geometry_free_per_chain.csv")
fig, axs = s.plt.subplots(1, 2, figsize=(10.5, 4.8), sharey=True)
rng = np.random.default_rng(42)
rhos = {}
for ax, grp in zip(axs, ("cl7", "cl4")):
    e = per[per.group == grp].groupby("entry_id")[["len_H12", "mobility_h12_z"]].mean().dropna()
    r, p = spearmanr(e.len_H12, e.mobility_h12_z)
    rhos[grp] = (r, p, len(e))
    ax.scatter(e.len_H12 + rng.uniform(-0.15, 0.15, len(e)), e.mobility_h12_z, s=22, color=s.GROUP_COLOR[grp],
               alpha=0.6, lw=0)
    gb = e.groupby(e.len_H12.round()).mobility_h12_z
    med = gb.median()[gb.size() >= 5]          # median only where a length holds at least 5 crystals
    ax.plot(med.index, med.values, color=s.INK, lw=1.6, label="median per length (≥ 5 crystals)")
    ax.set_title(f"{grp}: Spearman ρ {r:+.2f} (p {p:.1g}), {len(e)} crystals", fontsize=9.5, fontweight="normal",
                 color=s.INK2)
    ax.set_xlabel("H12 length from phi/psi (residues, mean over copies)")
    s.grid(ax, "y")
axs[0].set_ylabel("H12 mobility (C-alpha B z-score against the core)")
axs[1].legend(loc="upper right")
s.title(fig, f"In cl4 crystals a longer H12 is a less mobile H12 (ρ {rhos['cl4'][0]:+.2f}); in cl7, where H12 is "
             f"almost always complete, there is no link (ρ {rhos['cl7'][0]:+.2f}).", width=100)
fig.tight_layout()
s.save(fig, "h12_length_mobility")
