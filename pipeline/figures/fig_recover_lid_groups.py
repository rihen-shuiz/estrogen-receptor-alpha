"""Which clusterings recover the lid groups: with or without H12, with or without superposition?

Data: results/06_lid_hidden/no_lid_k_scan.csv (average linkage on RMSD without 528-550) and
results/07_geometry_free/geometry_free_kscan.csv (Ward on superposition-free features, all and without H12).
"""
import style as s

a = s.csv("results/06_lid_hidden/no_lid_k_scan.csv")
b = s.csv("results/07_geometry_free/geometry_free_kscan.csv")
series = [
    (b[b.features == "all_features"], "ari_vs_rmsd_groups", "helix COM distances and lengths, with H12 (no superposition)",
     s.INK, "-"),
    (b[b.features == "without_H12"], "ari_vs_rmsd_groups", "helix COM distances and lengths, without H12", s.MUTED,
     "-"),
    (a, "ari_vs_lid_groups", "RMSD without residues 528-550 (superposed)", s.VIOLET, "-"),
]
fig, ax = s.plt.subplots(figsize=(8.5, 5))
for d, col, lab, c, ls in series:
    ax.plot(d.k, d[col], color=c, lw=2, marker="o", ms=5, label=lab, linestyle=ls)
ax.set_xlabel("Number of clusters cut (K)")
ax.set_ylabel("Agreement with the lid groups (adjusted Rand index)")
ax.set_ylim(-0.03, 1.0)
ax.set_xticks(range(2, 11))
s.grid(ax, "y")
ax.legend(loc="upper right")
top = b[(b.features == "all_features") & (b.k == 2)].ari_vs_rmsd_groups.iloc[0]
nl = b[b.features == "without_H12"].ari_vs_rmsd_groups.max()
s.title(fig, f"Superposition-free helix geometry recovers the lid groups (ARI {top:.2f} at K = 2) only while H12 is "
             f"included; without H12 no cut does better than ARI {max(nl, a.ari_vs_lid_groups.max()):.2f}.", width=85)
s.save(fig, "recover_lid_groups")
