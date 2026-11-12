"""How far does a helix midpoint move when residues are missing at its ends?

Data: data/lid/midpoint_error_by_k.csv (median midpoint shift of complete helices truncated by k residues,
N and C ends pooled; summary of the 28 Oct truncation test).
"""
import style as s

v = s.csv("data/lid/midpoint_error_by_k.csv")
order = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]
cols = {"H12": s.CL4, "H4/5": s.VIOLET}
fig, ax = s.plt.subplots(figsize=(8.5, 5))
ends = []
for h in order:
    d = v[v.helix == h].set_index("k_removed").median_midpoint_error_A
    c = cols.get(h, s.MUTED)
    ax.plot(d.index, d.values, color=c, lw=2 if h in cols else 1.2, marker="o", ms=4)
    ends.append((h, d.index[-1], d.values[-1]))
by_x = {}
for h, x, y in ends:
    by_x.setdefault(x, []).append((h, y))
for x, items in by_x.items():
    ys = s.spread_labels([y for _, y in items], 0.05)
    for (h, _), y in zip(items, ys):
        ax.text(x + 0.1, y, h, fontsize=8.5, va="center", color=s.INK if h in cols else s.INK2)
ax.axhline(1.0, color=s.INK, lw=1)
ax.text(0.6, 1.03, "1 Å acceptance limit", fontsize=8.5, color=s.INK, va="bottom")
ax.set_xlabel("Residues removed from one end of a complete helix")
ax.set_ylabel("Median shift of the helix midpoint (Å)")
ax.set_xlim(0.5, 7)
s.grid(ax, "y")
h45 = v[v.helix == "H4/5"].set_index("k_removed").median_midpoint_error_A
first_fail = int(h45[h45 > 1.0].index.min())
s.title(fig, f"A helix midpoint stays within 1 Å when up to 4 end residues are missing; only H4/5 crosses 1 Å, at "
             f"{first_fail} removed, so midpoints let incomplete helices be compared.", width=85)
s.save(fig, "midpoint_error")
