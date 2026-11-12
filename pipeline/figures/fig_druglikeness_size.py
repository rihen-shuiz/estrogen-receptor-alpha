"""Which drug-likeness differences between cl7 and cl4 ligands remain once ligand size is accounted for?

Data: data/ligands/druglikeness_by_group.csv (q of the cl7-cl4 difference before and after adjusting for heavy
atoms, ligand and entry level; label = adds beyond size / restates size).
"""
import numpy as np

import style as s

d = s.csv("data/ligands/druglikeness_by_group.csv")
d = d[d.identical_to.isna()].copy()            # drop QED terms that repeat a plain parameter
d["q_raw"] = d[["lig_q", "ent_q"]].max(1)
d["q_adj"] = d[["lig_adj_q", "ent_adj_q"]].max(1)
d = d.sort_values("q_adj", ascending=False)
names = {"mw": "molecular weight", "clogp": "cLogP", "hbd": "H-bond donors", "hba": "H-bond acceptors",
         "lipinski_violations": "Lipinski violations", "rotb": "rotatable bonds", "tpsa": "polar surface area",
         "veber_violations": "Veber violations", "qed": "QED score", "qed_hba": "QED acceptor term",
         "qed_arom": "aromatic rings", "qed_alerts": "QED structural alerts", "brenk_n": "Brenk alerts",
         "pains_n": "PAINS alerts", "le": "ligand efficiency", "lipe": "lipophilic efficiency"}
fig, ax = s.plt.subplots(figsize=(8.5, 5.6))
for y, (_, r) in enumerate(d.iterrows()):
    adds = r.label == "adds beyond size"
    c = s.INK if adds else s.MUTED
    ax.plot([r.q_raw, r.q_adj], [y, y], color=s.GRID, lw=3, zorder=1, solid_capstyle="round")
    ax.scatter(r.q_raw, y, s=40, facecolor="white", edgecolor=c, lw=1.4, zorder=3)
    ax.scatter(r.q_adj, y, s=46, color=c, zorder=3)
ax.set_xscale("log")
ax.axvline(0.05, color=s.INK, lw=1)
ax.text(0.05, len(d) - 0.3, " q = 0.05", fontsize=8.5, color=s.INK, va="bottom")
ax.set_yticks(range(len(d)), [names.get(p, p) for p in d.parameter])
ax.scatter([], [], s=40, facecolor="white", edgecolor=s.MUTED, lw=1.4, label="cl7 vs cl4, unadjusted")
ax.scatter([], [], s=46, color=s.INK, label="adjusted for ligand size: adds beyond size")
ax.scatter([], [], s=46, color=s.MUTED, label="adjusted for ligand size: restates size")
ax.legend(loc="lower left", fontsize=8)
ax.set_xlabel("q value (larger of ligand and entry level; log scale)")
s.grid(ax, "x")
ax.set_xlim(1e-28, 3)
n_add = int((d.label == "adds beyond size").sum())
s.title(fig, f"Ligand size explains most drug-likeness differences between cl7 and cl4 ligands; {n_add} distinct "
             f"quantities (H-bond counts, rotatable bonds, polar surface, alert counts) differ beyond size.", width=85)
s.save(fig, "druglikeness_size")
