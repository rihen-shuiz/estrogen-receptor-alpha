"""Can action class be compared independently of lid position? How the uncurated ChEMBL classes fall into the groups.

Data: data/ligands/chain_class.csv (provisional ChEMBL-mechanism class per chain; an entry counts in a cell if at
least one of its chains has that class and group).
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

import style as s

cc = s.csv("data/ligands/chain_class.csv")
cc["action_class"] = cc.action_class.fillna("unlabelled")
classes = ["agonist", "SERM", "SERD", "unlabelled"]
groups = ["cl7", "cl4", "minor"]
m = np.array([[cc[(cc.action_class == k) & (cc.cluster == g)].entry_id.nunique() for g in groups] for k in classes])
tot = np.array([cc[cc.action_class == k].entry_id.nunique() for k in classes])
cmap = LinearSegmentedColormap.from_list("seq", ["#f7f6f3", "#cde2fb", "#6da7ec", "#256abf", "#0d366b"])
share = m / tot[:, None]
fig, ax = s.plt.subplots(figsize=(7.2, 4.6))
ax.pcolormesh(np.arange(4) - 0.5, np.arange(5) - 0.5, share, cmap=cmap, vmin=0, vmax=1, linewidth=0)
ax.set_ylim(3.5, -0.5)
for i in range(len(classes)):
    for j in range(len(groups)):
        ax.text(j, i, f"{m[i, j]}\n({share[i, j]:.0%})", ha="center", va="center", fontsize=9,
                color="white" if share[i, j] > 0.6 else s.INK)
ax.set_xticks(range(3), ["cl7\nlid over pocket", "cl4\nlid in groove", "minor groups"])
ax.set_yticks(range(4), [f"{k} ({t} entries)" for k, t in zip(classes, tot)])
ax.tick_params(length=0)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xlabel("Entries per class and lid group (share of the class's entries)")
s.title(fig, f"Labelled action classes almost coincide with the lid groups (agonist {m[0, 0]} of {tot[0]} entries in "
             f"cl7, SERM {m[1, 1]} of {tot[1]} and SERD {m[2, 1]} of {tot[2]} in cl4), so class cannot yet be "
             f"tested apart from lid position.", width=80)
s.save(fig, "class_vs_group")
