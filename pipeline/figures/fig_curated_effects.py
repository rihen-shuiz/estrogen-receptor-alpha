"""How large is each class contrast inside one lid position, against ordinary variation?

Data: results/12_curated_action/tables/effects_summary.csv (difference of medians of the lid-to-H11 distance between
two classes, PDB entry x ligand pairs, 95% CI from a bootstrap over DOI blocks; empty CI when either class has fewer
than 4 DOI blocks), results/03_lid_outliers/lid_distance_by_group.csv (ordinary variation per group).
"""
import style as s

eff = s.csv(f"{s.CURATED}/effects_summary.csv")
eff = eff[eff.metric == "dist_H12c_H11axis_A"].reset_index(drop=True)
ov = s.csv("results/03_lid_outliers/lid_distance_by_group.csv").set_index("group").ordinary_threshold
short = {"agonist": "agonist", "partial_agonist": "partial agonist", "antagonist": "antagonist"}
stratum = {"Y537S + peptide (exact)": "Y537S + peptide", "L536S group, no peptide": "L536S, no peptide",
           "L372S;L536S (exact), no pep": "L372S/L536S, no peptide", "all cl7 (crude)": "all chains",
           "all cl4 (crude)": "all chains"}
fig, ax = s.plt.subplots(figsize=(9.6, 5.0))
for y, r in eff.iterrows():
    thr = ov[r.group]
    ax.fill_betweenx([y - 0.4, y + 0.4], -thr, thr, color=s.BAND, lw=0)
    if r.ci_valid in (True, "True"):
        ax.plot([r.ci95_lo, r.ci95_hi], [y, y], color=s.INK2, lw=2)
    ax.scatter(r.diff_a_minus_b, y, s=60, color=s.INK if r.exploratory in (False, "False") else "white",
               edgecolor=s.INK, lw=1.4, zorder=3)
    ax.text(1.32, y, f"{int(r.n_a)} vs {int(r.n_b)} pairs", va="center", fontsize=9.5, color=s.INK2)
ax.set_yticks(range(len(eff)), [f"{r.group}, {stratum[r.stratum]}: {short[r.class_a]} − {short[r.class_b]}"
                                for r in eff.itertuples()])
ax.invert_yaxis()
ax.axvline(0, color=s.MUTED, lw=0.8)
ax.set_xlim(-1.3, 1.3)
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_xlabel("Difference of median lid-to-H11 distance (Å); grey band: ordinary variation of the group")
s.grid(ax, "x")
s.title(fig, "No class contrast reaches ordinary variation (±0.91 Å in cl7, ±1.08 Å in cl4); filled: the planned "
             "comparable strata, open: all chains of the group; bars: 95% CI over DOI blocks.", width=95)
s.save(fig, "curated_effects")
