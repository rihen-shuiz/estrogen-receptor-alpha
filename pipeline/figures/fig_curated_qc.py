"""Do the clusters survive when low-quality structures are removed?

Data: results/12_curated_action/tables/qc_clustering_sensitivity.csv (average linkage on the saved C-alpha RMSD matrix
restricted to the kept chains, at each cut; ARI against the full 788-chain clustering on the shared chains).
"""
import style as s

q = s.csv(f"{s.CURATED}/qc_clustering_sensitivity.csv")
sets = [("rfree_le_0p30", "R-free ≤ 0.30", s.INK), ("qc_all3_tierA", "resolution, R-free and B-factor all pass",
                                                     s.VIOLET)]
fig, ax = s.plt.subplots(figsize=(7.6, 4.4))
ys = s.spread_labels([q[q.subset == k].ari_vs_full.iloc[-1] for k, _, _ in sets], 0.02)
for (key, lab, c), y in zip(sets, ys):
    d = q[q.subset == key]
    ax.plot(d.cut_A, d.ari_vs_full, color=c, marker="o", ms=7, lw=2.4)
    ax.text(2.55, y, f"{lab}\n({int(d.n_chains.iloc[0])} of 788 chains)", va="center", fontsize=10, color=s.INK)
ax.set_xticks([1.0, 1.5, 2.0, 2.5], ["1.0", "1.5", "2.0", "2.5"])
ax.set_xlim(0.9, 2.5)
ax.set_ylim(0.88, 1.005)
ax.set_xlabel("Clustering cut (Å, average linkage on C-alpha RMSD)")
ax.set_ylabel("Agreement with the full clustering (ARI)")
s.grid(ax, "y")
a = q[(q.subset == "rfree_le_0p30") & (q.cut_A == 1.5)].ari_vs_full.iloc[0]
b = q[(q.subset == "qc_all3_tierA") & (q.cut_A == 1.5)].ari_vs_full.iloc[0]
s.title(fig, f"Dropping low-quality structures leaves the clusters in place (ARI {a:.3f} and {b:.3f} at 1.5 Å): the "
             "groups are not built by poor models.", width=75)
s.save(fig, "curated_qc")
