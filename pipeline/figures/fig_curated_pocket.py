"""Does the pocket differ between action classes inside one lid position?

Data: results/12_curated_action/tables/pairs_by_group.csv (pocket volume as % of the E2 pocket, per PDB entry x
ligand x group; pairs whose ligand is fully modelled), results/12_curated_action/tables/pocket_by_class.csv (medians,
Kruskal-Wallis across classes with at least 3 pairs).
"""
import numpy as np

import style as s

pg = s.csv(f"{s.CURATED}/pairs_by_group.csv")
kw = s.csv(f"{s.CURATED}/pocket_by_class.csv").set_index(["group", "metric"])
pg = pg[pg.action_class.isin(s.CLASSES) & pg.pocket_complete_primary.astype(str).eq("True")]
col = "pocket_pct_of_E2_pocket"
rng = np.random.default_rng(42)
fig, axs = s.plt.subplots(1, 2, figsize=(10.5, 4.6), sharey=True)
for ax, g in zip(axs, ("cl7", "cl4")):
    d = pg[pg.group == g]
    ax.axhline(100, color=s.MUTED, lw=1, ls="--")
    for x, cls in enumerate(s.CLASSES):
        v = d[d.action_class == cls][col].dropna()
        if v.empty:
            continue
        ax.scatter(x + rng.uniform(-0.2, 0.2, len(v)), v, s=24, marker=s.CLASS_MARKER[cls], facecolor="none",
                   edgecolor=s.CLASS_COLOR[cls], lw=1.0)
        ax.plot([x - 0.3, x + 0.3], [v.median()] * 2, color=s.INK, lw=2.4)
        ax.text(x, 182, f"n = {len(v)}", ha="center", fontsize=9.5, color=s.INK2)
    p = kw.loc[(g, col), "kruskal_p"]
    ax.set_title(f"{g}: Kruskal-Wallis p = {p:.1g}", fontsize=10.5, fontweight="normal", color=s.INK2)
    ax.set_xticks(range(3), [s.CLASS_LABEL[c] for c in s.CLASSES])
    ax.set_xlim(-0.6, 2.6)
    s.grid(ax, "y")
axs[0].set_ylabel("Pocket volume (% of the E2 pocket)")
axs[0].set_ylim(20, 195)
axs[0].text(2.55, 101, "E2 pocket", ha="right", va="bottom", fontsize=9, color=s.INK2)
r = kw.loc[("cl7", col)]
s.title(fig, f"In the closed position the pocket grows from agonist to antagonist (medians {r.median_agonist:.0f}%, "
             f"{r.median_partial_agonist:.0f}% and {r.median_antagonist:.0f}% of the E2 pocket), while the share the "
             "ligand fills stays about one third.", width=95)
fig.subplots_adjust(top=0.86, wspace=0.08)
s.save(fig, "curated_pocket")
