"""Do helix lengths read from phi/psi differ between the lid groups by more than copies of one crystal differ?

Data: results/07_geometry_free/geometry_free_per_chain.csv (phi/psi helix length per chain) and
results/07_geometry_free/geometry_free_features.csv (copy-to-copy 95th percentile per length).
"""
import numpy as np

import style as s

per = s.csv("results/07_geometry_free/geometry_free_per_chain.csv")
ft = s.csv("results/07_geometry_free/geometry_free_features.csv").set_index("feature")
hel = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]
fig, ax = s.plt.subplots(figsize=(9.5, 5))
x = np.arange(len(hel))
for k, (grp, off) in enumerate((("cl7", -0.14), ("cl4", 0.14))):
    d = per[per.group == grp]
    med = [d[f"len_{h}"].median() for h in hel]
    q1 = [d[f"len_{h}"].quantile(0.25) for h in hel]
    q3 = [d[f"len_{h}"].quantile(0.75) for h in hel]
    ax.vlines(x + off, q1, q3, color=s.GROUP_COLOR[grp], lw=2.4, alpha=0.6)
    ax.scatter(x + off, med, s=56, color=s.GROUP_COLOR[grp], zorder=3, label=f"{s.GROUP_LABEL[grp]} (median, IQR)")
for i, h in enumerate(hel):
    r = ft.loc[f"len_{h}"]
    ax.text(i, 1.5, f"Δ {r.diff_cl4_minus_cl7:+.0f}\ncopies ≤ {r.copy_p95:.0f}", ha="center", fontsize=7.5,
            color=s.INK if str(r.real) == "True" else s.INK2)
ax.set_xticks(x, hel)
ax.set_ylim(0, 38)
ax.set_ylabel("Helix length from phi/psi (residues)")
s.grid(ax, "y")
ax.legend(loc="upper left")
big = ft[(ft.kind == "length")].assign(a=lambda d: d.diff_cl4_minus_cl7.abs()).sort_values("a", ascending=False)
names = ", ".join(f"{i[4:]} {r.diff_cl4_minus_cl7:+.0f}" for i, r in big.head(3).iterrows())
s.title(fig, f"Helix ends shift by 3-5 residues between the lid groups ({names}), but copies of one crystal differ "
             f"by as much, so no length difference passes the ordinary-variation rule.", width=95)
s.save(fig, "helix_lengths")
