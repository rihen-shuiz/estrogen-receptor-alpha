"""Which chains sit unusually far from or near to H11, and are they beyond ordinary variation?

Data: results/03_lid_outliers/lid_distance_per_chain.csv (lid-to-H11 distance per chain), results/03_lid_outliers/lid_distance_by_group.csv
(medians and ordinary-variation thresholds), results/03_lid_outliers/lid_outliers.csv (examined chains).
"""
import numpy as np

import style as s

g = s.csv("results/03_lid_outliers/lid_distance_per_chain.csv")
summ = s.csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group")
fl = s.csv("results/03_lid_outliers/lid_outliers.csv")
col = "lid_to_h11_A"
rows = {"cl7": 2, "minor": 1, "cl4": 0}
rng = np.random.default_rng(42)

fig, ax = s.plt.subplots(figsize=(10, 5.2))
for grp in ("cl7", "cl4"):
    m, thr = summ.loc[grp, "median"], summ.loc[grp, "ordinary_threshold"]
    ax.axvspan(m - thr, m + thr, color=s.GROUP_COLOR[grp], alpha=0.08, lw=0)
    ax.text(m, 2.62, f"{grp} median {m:.1f} Å ± {thr:.2f} Å\n(copies of one crystal, 95th pct)", ha="center",
            va="bottom", fontsize=8, color=s.INK2)
flagged = set(fl.chain_key)
for grp, y in rows.items():
    d = g[(g.group == grp) & g[col].notna()]
    jit = rng.uniform(-0.18, 0.18, len(d))
    normal = ~d.chain_key.isin(flagged)
    ax.scatter(d[col][normal], y + jit[normal.to_numpy()], s=14, color=s.GROUP_COLOR[grp], alpha=0.45, lw=0,
               marker=s.GROUP_MARKER[grp])
    ax.scatter(d[col][~normal], y + jit[(~normal).to_numpy()], s=46, color=s.GROUP_COLOR[grp],
               marker=s.GROUP_MARKER[grp], edgecolor=s.INK, linewidth=1.0, zorder=3)
    lab = d[~normal].assign(yy=y + jit[(~normal).to_numpy()])
    s.place_labels(ax, lab[col].tolist(), lab.yy.tolist(),
                   [k.rsplit("|", 1)[0].replace("|", " ") for k in lab.chain_key], y, min_gap=0.75)
ax.set_yticks([0, 1, 2], ["cl4\n(lid in groove)", "minor groups\n(grey)", "cl7\n(lid over pocket)"])
ax.set_ylim(-0.6, 3.1)
ax.set_xlabel("Distance from the H12 centre to the H11 axis (Å)")
s.grid(ax, "x")
n_beyond = int(fl.beyond_ordinary.astype(str).eq("True").sum())
n_between = int((fl.lid_state == "between cl7 and cl4").sum())
s.title(fig, f"{n_beyond} chains sit beyond the ordinary spread of their lid group; {n_between} of them, "
             f"mostly grey, hold the lid between the two positions. Ringed: chains examined one by one.")
s.save(fig, "lid_outliers")
