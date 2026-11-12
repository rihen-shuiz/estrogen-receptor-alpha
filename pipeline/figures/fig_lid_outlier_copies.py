"""Do the other copies in the same crystal share an outlier's unusual lid position?

Data: results/03_lid_outliers/lid_distance_per_chain.csv, results/03_lid_outliers/lid_outliers.csv, results/03_lid_outliers/lid_distance_by_group.csv.
"""
import numpy as np

import style as s

g = s.csv("results/03_lid_outliers/lid_distance_per_chain.csv")
fl = s.csv("results/03_lid_outliers/lid_outliers.csv")
summ = s.csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group")
col = "lid_to_h11_A"
order = {"cl7": 0, "minor": 1, "cl4": 2}
ents = fl.assign(o=fl.group.map(order)).sort_values(["o", "lid_to_h11_A"]).entry_id.drop_duplicates().tolist()
ents = [e for e in ents if (g.entry_id == e).sum() > 1]  # crystals with more than one copy

fig, ax = s.plt.subplots(figsize=(9.5, 6.4))
for grp in ("cl7", "cl4"):
    m, thr = summ.loc[grp, "median"], summ.loc[grp, "ordinary_threshold"]
    ax.axvspan(m - thr, m + thr, color=s.GROUP_COLOR[grp], alpha=0.08, lw=0)
flagged = set(fl.chain_key)
n_alone = 0
for y, e in enumerate(ents[::-1]):
    d = g[g.entry_id == e]
    v = d[col].dropna()
    if len(v) > 1:
        ax.plot([v.min(), v.max()], [y, y], color=s.GRID, lw=2, zorder=1, solid_capstyle="round")
    for _, r in d.iterrows():
        if np.isnan(r[col]):
            continue
        f = r.chain_key in flagged
        ax.scatter(r[col], y, s=58 if f else 34, color=s.GROUP_COLOR[r.group] if f else "white",
                   edgecolor=s.INK if f else s.GROUP_COLOR[r.group], linewidth=1.2 if f else 1.4,
                   marker=s.GROUP_MARKER[r.group], zorder=3)
        if f or len(d) <= 4:   # copy letters only where they stay readable
            ax.annotate(r.chain_key.split("|")[1], (r[col], y), (0, 7), textcoords="offset points", fontsize=7,
                        color=s.INK2, ha="center")
    n_na = int(d[col].isna().sum())
    # other copies inside a main group's ordinary range -> the outlier is copy-specific
    others = d[~d.chain_key.isin(flagged)][col].dropna()
    inside = others.apply(lambda x: any(abs(x - summ.loc[k, "median"]) <= summ.loc[k, "ordinary_threshold"]
                                        for k in ("cl7", "cl4")))
    if len(others) and inside.all():
        n_alone += 1
    lab = f"{e}" + (f"  ({n_na} {'copy' if n_na == 1 else 'copies'} without H12 axis)" if n_na else "")
    ax.text(10.3, y, lab, ha="left", va="center", fontsize=8, color=s.INK)
ax.set_yticks([])
ax.spines["left"].set_visible(False)
ax.set_xlim(10.2, 24.8)
ax.set_ylim(-0.7, len(ents) - 0.3)
ax.set_xlabel("Distance from the H12 centre to the H11 axis (Å); shaded: cl7 (blue) and cl4 (orange) ordinary range")
s.grid(ax, "x")
ax.scatter([], [], s=58, color=s.MUTED, edgecolor=s.INK, label="examined chain (filled)")
ax.scatter([], [], s=34, color="white", edgecolor=s.MUTED, linewidth=1.4, label="other copy in the same crystal")
ax.legend(loc="upper center", bbox_to_anchor=(0.62, 1.0))
sp = g[g.entry_id == "9B25"][col]
s.title(fig, f"In {n_alone} of {len(ents)} crystals the other copies keep an ordinary lid, so the unusual position "
             f"belongs to one copy; the {sp.notna().sum()} copies of 9B25 spread from {sp.min():.0f} to {sp.max():.0f} Å.")
s.save(fig, "lid_outlier_copies")
