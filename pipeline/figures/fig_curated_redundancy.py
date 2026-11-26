"""Is each action class as redundant as the whole set? Conformations found as PDB entries are added at random.

Data: results/12_curated_action/tables/redundancy_accumulation_curves.csv (mean and 10th-90th percentile of the number
of 1.5 A clusters seen after adding n PDB entries in random order, 200 orders; an entry brings the clusters of all its
chains), results/12_curated_action/tables/redundancy_hill_summary.csv (Hill q1 with each PDB entry weighted equally).
"""
import style as s

cur = s.csv(f"{s.CURATED}/redundancy_accumulation_curves.csv")
hill = s.csv(f"{s.CURATED}/redundancy_hill_summary.csv").set_index(["subset", "cut_A"])
cur = cur[cur.cut_A == 1.5]
fig, ax = s.plt.subplots(figsize=(8.2, 4.6))
ends = []
for cls in s.CLASSES:
    d = cur[cur.subset == cls]
    ax.fill_between(d.n_pdbs, d.p10_clusters, d.p90_clusters, color=s.CLASS_COLOR[cls], alpha=0.12, lw=0)
    ax.plot(d.n_pdbs, d.mean_clusters, color=s.CLASS_COLOR[cls], lw=2.4)
    h = hill.loc[(cls, 1.5)]
    ends.append((d.mean_clusters.iloc[-1], f"{s.CLASS_LABEL[cls]}: {int(h.n_pdbs)} entries, q1 = {h.hill_q1:.2f}",
                 d.n_pdbs.iloc[-1]))
ys = s.spread_labels([e[0] for e in ends], 0.32)
for (y0, lab, x), y in zip(ends, ys):
    ax.text(x + 2, y, lab, va="center", fontsize=10, color=s.INK)
ax.set_xlim(0, 175)
ax.set_xlabel("PDB entries added in random order")
ax.set_ylabel("Clusters (1.5 Å) seen, mean of 200 orders")
s.grid(ax, "y")
al = hill.loc[("all_351_pdbs", 1.5)]
s.title(fig, f"Each class holds one to two effective conformations (Hill q1, entries weighted equally), against "
             f"{al.hill_q1:.2f} for all {int(al.n_pdbs)} entries: most new structures repeat a known state.", width=85)
s.save(fig, "curated_redundancy")
