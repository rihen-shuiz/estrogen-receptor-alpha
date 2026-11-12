"""How much does "incomplete H12" depend on where the H12 window ends?

Data: data/lid/h12_window_status.csv (chain status on 538-547 and on the data-derived 538-544).
"""
import style as s

st = s.csv("data/lid/h12_window_status.csv")
st = st[st.level == "chain"]
fig, ax = s.plt.subplots(figsize=(7.5, 5))
vals = {}
for grp in ("cl4", "cl7", "minor"):
    d = st[st.cluster_group == grp]
    vals[grp] = (int(d[d.h12_status_538_547 != "complete"].n.sum()), int(d[d.h12_status_538_544 != "complete"].n.sum()),
                 int(d.n.sum()))
    ax.plot([0, 1], vals[grp][:2], color=s.GROUP_COLOR[grp], marker=s.GROUP_MARKER[grp], ms=8, lw=2)
grps = list(vals)
ya = s.spread_labels([vals[g][0] for g in grps], 9)
yb = s.spread_labels([vals[g][1] for g in grps], 9)
for g, a_, b_ in zip(grps, ya, yb):
    ax.text(-0.04, a_, f"{g}: {vals[g][0]} of {vals[g][2]}", ha="right", va="center", fontsize=9, color=s.INK)
    ax.text(1.04, b_, f"{vals[g][1]}", ha="left", va="center", fontsize=9, color=s.INK)
ax.set_xticks([0, 1], ["H12 window 538-547\n(primary)", "H12 window 538-544\n(from consensus secondary structure)"])
ax.set_xlim(-0.55, 1.25)
ax.set_ylabel("Chains with an incomplete H12 (partial or not modelled)")
s.grid(ax, "y")
d = st[st.cluster_group == "cl4"]
a = int(d[d.h12_status_538_547 != "complete"].n.sum())
b = int(d[d.h12_status_538_544 != "complete"].n.sum())
s.title(fig, f"Ending the H12 window at 544 instead of 547 makes an incomplete H12 rare (cl4 {a} → {b} chains): most "
             f"'partial' lids lack only residues 545-547.", width=80)
s.save(fig, "h12_window")
