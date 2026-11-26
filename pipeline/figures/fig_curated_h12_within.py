"""Inside one lid position, does the action class move the lid?

Data: results/12_curated_action/tables/pairs_by_group.csv (median lid-to-H11 distance of the chains of each PDB entry
x ligand pair in a group), results/12_curated_action/tables/effects_summary.csv (median differences, DOI-block
bootstrap 95% CI), results/03_lid_outliers/lid_distance_by_group.csv (ordinary variation of the distance).
Strata: cl7 with the Y537S construct and a bound coactivator peptide; cl4 with an L536S construct and no peptide.
"""
import numpy as np

import style as s

pg = s.csv(f"{s.CURATED}/pairs_by_group.csv")
eff = s.csv(f"{s.CURATED}/effects_summary.csv").set_index("contrast_id")
ov = s.csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group").ordinary_threshold
col = "dist_H12c_H11axis_A"
pep = pg.peptide_bound.astype(str)
strata = [("cl7", "cl7: Y537S construct, peptide bound",
           pg[(pg.group == "cl7") & (pg.construct_mutations == "Y537S") & pep.eq("True")]),
          ("cl4", "cl4: L536S construct, no peptide",
           pg[(pg.group == "cl4") & (pg.mutation_status == "L536S") & pep.eq("False")])]
rng = np.random.default_rng(42)
fig, axs = s.plt.subplots(1, 2, figsize=(10.5, 4.8))
for ax, (g, ttl, d) in zip(axs, strata):
    d = d[d.action_class.isin(s.CLASSES)].dropna(subset=[col])
    ref = d[d.action_class == ("agonist" if g == "cl7" else "partial_agonist")][col].median()
    for yy in (ref - ov[g], ref + ov[g]):
        ax.axhline(yy, color=s.MUTED, lw=1.2, ls="--")
    ax.set_ylim(ref - ov[g] - 0.45, ref + ov[g] + 0.55)
    for x, cls in enumerate(s.CLASSES):
        v = d[d.action_class == cls][col]
        if v.empty:
            continue
        ax.scatter(x + rng.uniform(-0.2, 0.2, len(v)), v, s=26, marker=s.CLASS_MARKER[cls], facecolor="none",
                   edgecolor=s.CLASS_COLOR[cls], lw=1.0)
        ax.plot([x - 0.3, x + 0.3], [v.median()] * 2, color=s.INK, lw=2.4)
        ax.text(x, ref + ov[g] + 0.12, f"n = {len(v)}", ha="center", va="bottom", fontsize=9.5, color=s.INK2)
    ax.set_xticks(range(3), [s.CLASS_LABEL[c] for c in s.CLASSES])
    ax.set_xlim(-0.6, 2.6)
    ax.set_title(ttl, fontsize=10.5, fontweight="normal", color=s.INK2)
    s.grid(ax, "y")
axs[0].set_ylabel("Lid (H12) centre to H11 axis (Å), median per pair")
for ax, (g, _, _) in zip(axs, strata):
    ax.text(2.55, ax.get_ylim()[0] + 0.08, f"dashed: reference median ± {ov[g]:.2f} Å (ordinary variation)",
            ha="right", va="bottom", fontsize=9, color=s.INK2)
a = eff.loc["cl7_Y537S_pep: ant vs ago"]
b = eff.loc["cl7_Y537S_pep: part vs ago"]
c4 = eff.loc["cl4_L536S_nopep: ant vs part"]
s.title(fig, f"Inside one lid position the class moves H12 by at most {abs(a.diff_a_minus_b):.2f} Å (cl7: antagonist "
             f"{s.signed(a.diff_a_minus_b)}, partial agonist {s.signed(b.diff_a_minus_b)} Å vs agonist; cl4: antagonist "
             f"{s.signed(c4.diff_a_minus_b)} Å vs partial agonist), far below copy-to-copy variation.", width=95)
fig.subplots_adjust(top=0.88, wspace=0.16)
s.save(fig, "curated_h12_within")
