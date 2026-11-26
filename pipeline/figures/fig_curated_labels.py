"""How many ligands carry an action class now, against the ChEMBL labels used in October?

Data: data/ligands/action_class.csv (ChEMBL mechanism labels per entry x ligand), results/12_curated_action/tables/
join_audit_summary.json (hand labels and their link to the structures).
"""
import json

import style as s

old = s.csv("data/ligands/action_class.csv")
j = json.loads(s._path(f"{s.CURATED}/join_audit_summary.json").read_text())
n_old = int(old.action_class.isin(["agonist", "SERM", "SERD"]).sum())
rows = [("ChEMBL labels (October)", n_old, len(old), s.MUTED),
        ("Hand labels: agonist, partial agonist\nor antagonist", j["curated_core3_rows"], j["curated_real_rows"], s.INK),
        ("... linked to one ligand\nin the crystal structure", j["core3_strict_pairs"], j["curated_real_rows"], s.INK2)]
fig, ax = s.plt.subplots(figsize=(8.2, 3.4))
for y, (lab, k, n, c) in enumerate(rows):
    ax.barh(y, n, height=0.55, color=s.BAND)
    ax.barh(y, k, height=0.55, color=c)
    ax.text(k + 4, y, f"{k} of {n}", va="center", fontsize=10.5, color=s.INK)
ax.set_yticks(range(len(rows)), [r[0] for r in rows])
ax.invert_yaxis()
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlabel("PDB entry × ligand pairs (grey: all pairs considered)")
s.grid(ax, "x")
by = j["core3_strict_by_class"]
s.title(fig, f"Hand labelling raised the labelled pairs from {n_old} to {j['curated_core3_rows']}; {j['core3_strict_pairs']} "
             f"are tied to a structure ({by['agonist']} agonist, {by['partial_agonist']} partial agonist, "
             f"{by['antagonist']} antagonist).", width=85)
s.save(fig, "curated_labels")
