"""SV01b: write results/sv/SV01b/SUMMARY.md (rule 1.7) from the saved SV01b (and SV01a) files only.
Run: .venv/bin/python -B scripts/sv/SV01b_summary.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV01b"
NOT_DONE = (
    "Not done / substituted: (1) step 8 (medoid H12 clash test) not run: the step runs it only if the verdict is "
    "'mixed' and the re-issued verdict is '{step}' (scripts/sv/SV01b_clash.py is written and gated, "
    "clash_test_not_run.md; it was run once outside results/ as a smoke test on the 235 partial/not-modelled chains "
    "and that output was discarded, not used); (2) choices the run file leaves open, fixed before fitting (docstrings of "
    "SV01b_models.py and common/lattice_contacts.py): residue B = C-alpha B (altloc: highest occupancy, then blank/A); "
    "mobility = mean z over every observed residue of 528-547; core SD with ddof 0 (existing regional-B code); "
    "core v3 used unchanged, so 528-531 are in both the core and the window; lattice terms = number of heavy atoms "
    "of 528-550 within 4.0 A of an ASU copy / a symmetry copy of another chain (polymer chains; waters and ligands "
    "removed as in SV00), with the chain's own bound peptide copy left out because it is 'peptide' in block B; "
    "(3) re-issued step verdict = 'artifact'/'biology' only if the [DEFAULT] condition holds for both outcomes "
    "(complete, mobility) at both levels, else 'mixed' (the run file names one verdict for two outcomes); "
    "(4) mobility models: chain level OLS with entry-clustered sandwich errors, entry level OLS on the first matrix "
    "chain of the entry among chains with the outcome, block tests cluster-robust Wald chi2 / nested F, Wald 95 % CI "
    "b +- 1.96 SE; (5) inherited unchanged from SV01a: median imputation + 'not reported' indicators, entry mean B "
    "as 'mean B', the drop rule extended to collinearity and joint separation, McFadden gains, one BH family.")


def rd(name, d=OUT):
    return pd.read_csv(d / name, low_memory=False)


def size(p):
    s = p.stat().st_size
    return f"{s / 1e3:.1f} kB" if s < 1e6 else f"{s / 1e6:.2f} MB"


def main():
    nums = json.loads((OUT / "numbers.json").read_text())
    lchk = json.loads((OUT / "lattice_checks.json").read_text())
    rep = json.loads((OUT / "reproduce_SV01a.json").read_text())
    mb, orr, vt = rd("model_blocks.csv"), rd("odds_ratios.csv"), rd("verdict_table.csv")
    lg, mob = rd("lattice_contacts_by_group.csv"), rd("mobility_per_chain.csv")
    pr = rd("lattice_contacts_per_residue.csv")
    a_mb = rd("model_blocks.csv", ROOT / "results/sv/SV01a")
    step, ov, lv = nums["step_verdict"], nums["outcome_verdicts"], nums["level_verdicts"]
    L = [NOT_DONE.format(step=step), "", "# SV01b - H12 completeness and mobility with coordinates (lattice contacts)", ""]
    L += [f"**Re-issued verdict ([DEFAULT] rule, lattice terms in block X): {step}** (SV01a: "
          f"{nums['sv01a_step_verdict']['complete']}). Complete H12 (538-547): {ov['complete']} (chain "
          f"{lv['complete|chain']}, entry {lv['complete|entry']}); mobility (B z-score 528-547): {ov['mobility']} "
          f"(chain {lv['mobility|chain']}, entry {lv['mobility|entry']}). Sensitivity: fraction observed "
          f"{ov['fraction_observed']}, mobility without TLS-refined entries {ov['mobility_no_tls']}. Step 8 not run "
          "(verdict not mixed). The change from SV01a comes from the lattice terms, which for the completeness "
          "outcome are counted on modelled atoms (see Reading).", ""]
    L += ["## Inputs", "",
          "- SV00 (DONE): results/sv/SV00/chain_table.csv (788 chains, tls_refined, auth_asym_id), observed_mask.npy, "
          "data/sv/backbone/backbone.npz (C-alpha B), restored mmCIFs data/raw/, regenerated Q04 tables "
          "data/sv/lbd_extract/ (residue mapping, chain identities).",
          "- SV01a (DONE): scripts/sv/SV01a_completeness.py imported (load, predictors, imputation, drop rule, fit); "
          f"refitting the SV01a predictor set with the SV01b code reproduces results/sv/SV01a: {rep['model_block_rows']}"
          f"/{rep['model_block_rows_sv01a']} block rows, max |d gain| {rep['max_abs_diff_gain']:.1e}, max |d p| "
          f"{rep['max_abs_diff_block_p']:.1e}; {rep['odds_rows']} odds rows, status equal {rep['status_equal']}, "
          f"max |d log OR| {rep['max_abs_diff_log_or']:.1e} (reproduce_SV01a.json).",
          "- Core v3: esr_analysis/resources/config/alignment_core_v3.json core_ranges (the S03 core) "
          + ", ".join(f"{a}-{b}" for a, b in nums["core_ranges"]) + ".",
          f"- C-alpha B in backbone.npz equals residue_mapping B_factor_A2 for {nums['b_crosscheck']['n_compared']} "
          f"observed C-alpha (max |d| {nums['b_crosscheck']['max_abs_diff']}).", ""]

    # step 6
    odd = pr.set_index(["chain_key", "esr1_position"])
    odd_keys = [tuple([k.split(":")[0], int(k.split(":")[1])]) for k in lchk["found_not_ca_observed"]]
    odd_contacts = int(odd.loc[odd_keys, ["asu_atoms", "sym_atoms"]].to_numpy().sum()) if odd_keys else 0
    L += ["## 6. Lattice contacts (scripts/sv/common/lattice_contacts.py, lattice_contacts_per_residue.csv / _per_chain.csv)", "",
          f"- {lchk['n_chains']} chains in {lchk['n_entries']} entries; crystal frame {lchk['crystal_frame']}; "
          f"own bound peptide copy agrees with SV00 chain_peptide.csv for {lchk['own_peptide_agrees_with_SV00']}/"
          f"{lchk['n_chains']} chains ({lchk['own_peptide_source']}).",
          f"- Residues 528-550 with atoms: {lchk['residues_528_550_found']}; with C-alpha observed (SV00 mask): "
          f"{lchk['residues_528_550_ca_observed_mask']}; {len(lchk['found_not_ca_observed'])} residues have atoms but no "
          f"observed C-alpha ({', '.join(lchk['found_not_ca_observed'])}; lattice contact atoms on them: {odd_contacts}).",
          f"- Chains with any ASU-partner contact on 528-550: {lchk['chains_with_asu_contact']}; with any symmetry-mate "
          f"contact: {lchk['chains_with_sym_contact']}.", "",
          "By cluster and H12 status (lattice_contacts_by_group.csv; mean heavy atoms of 528-550 in contact, and where):", "",
          "| group | H12 | chains | entries | any ASU | any sym | mean ASU atoms | mean sym atoms | sym 528-537 | sym 538-545 | sym 546-550 | modelled atoms 546-550 |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in lg[lg.cluster_group.isin(["cl7", "cl4"])].itertuples():
        L.append(f"| {r.cluster_group} | {r.h12_status} | {r.n_chains} | {r.n_entries} | {r.chains_asu_contact} | "
                 f"{r.chains_sym_contact} | {r.mean_asu_atoms:.2f} | {r.mean_sym_atoms:.2f} | {r.mean_sym_atoms_528_537:.2f} | "
                 f"{r.mean_sym_atoms_538_545:.2f} | {r.mean_sym_atoms_546_550:.2f} | {r.mean_modelled_heavy_atoms_546_550:.1f} |")
    L.append("")

    # step 5
    m = mob[mob.s08_cluster.isin([7, 4])]
    L += ["## 5. Local mobility (mobility_per_chain.csv, mobility_per_residue.csv)", "",
          f"- mobility_z defined for cl7 {nums['mobility_defined']['cl7']}/{int((m.cluster_group == 'cl7').sum())} and "
          f"cl4 {nums['mobility_defined']['cl4']}/{int((m.cluster_group == 'cl4').sum())} chains (undefined = no residue "
          f"of 528-547 observed: {', '.join(nums['mobility_undefined_chains'])}); core with < 2 valid B: "
          f"{len(nums['core_not_ok'])}.",
          f"- TLS-refined (SV00 tls_refined): cl7 {int(m[m.cluster_group == 'cl7'].tls_refined.sum())}, cl4 "
          f"{int(m[m.cluster_group == 'cl4'].tls_refined.sum())} chains; the TLS sensitivity keeps "
          f"{int((m.mobility_defined & ~m.tls_refined & (m.cluster_group == 'cl7')).sum())} cl7 and "
          f"{int((m.mobility_defined & ~m.tls_refined & (m.cluster_group == 'cl4')).sum())} cl4 chains.", "",
          "| group | H12 | chains with mobility | median mobility_z | median loop 528-537 z | median H12 538-547 z |",
          "|---|---|---|---|---|---|"]
    for (g, st), x in m[m.mobility_defined].groupby(["cluster_group", "h12_status"]):
        f2 = lambda v: "-" if pd.isna(v) else f"{v:.2f}"  # noqa: E731
        L.append(f"| {g} | {st} | {len(x)} | {f2(x.mobility_z.median())} | {f2(x.mobility_loop_z.median())} | "
                 f"{f2(x.mobility_h12_z.median())} |")
    L.append("")

    # step 7
    L += ["## 7. Refit with lattice terms (model_blocks.csv, odds_ratios.csv, verdict_table.csv)", "",
          "Gain entered first / entered last (McFadden R2 for complete and fraction observed, OLS R2 for mobility); "
          "test and BH q of the block entered last:", "",
          "| outcome | cluster | level | n | Q | X | B | Q+X | test Q+X (last) | q Q+X | test B (last) | q B |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (o, c, lvl), d in mb.groupby(["outcome", "cluster", "level"], sort=False):
        g = {(r.block, r.order): r.gain for r in d.itertuples()}
        last = d[d.order == "last"].set_index("block")

        def tst(b):
            r = last.loc[b]
            return (f"{r.test.split(' (')[0]} {r.statistic:.1f} (df {int(r.df)}), p {r.p_value:.2g}"
                    if pd.notna(r.p_value) else "-"), (f"{r.q_value:.2g}" if pd.notna(r.q_value) else "-")
        tqx, qqx = tst("Q+X")
        tb, qb = tst("B")
        L.append(f"| {o} | cl{c} | {lvl} | {int(d.n.iloc[0])} | " + " | ".join(
            f"{g[(b, 'first')]:.3f} / {g[(b, 'last')]:.3f}" for b in ("Q", "X", "B", "Q+X")) +
            f" | {tqx} | {qqx} | {tb} | {qb} |")
    L += ["", f"BH: one family over every p-value of SV01b ({nums['n_tests_bh']} tests).", ""]
    L.append("Q+X gain entered first, complete outcome, SV01a (no lattice terms) -> SV01b: " + "; ".join(
        f"cl{c} {lvl} {a_mb[(a_mb.outcome == 'complete') & (a_mb.cluster == c) & (a_mb.level == lvl) & (a_mb.order == 'first') & (a_mb.block == 'Q+X')].mcfadden_gain.iloc[0]:.3f} -> "
        f"{mb[(mb.outcome == 'complete') & (mb.cluster == c) & (mb.level == lvl) & (mb.order == 'first') & (mb.block == 'Q+X')].gain.iloc[0]:.3f}"
        for lvl in ("chain", "entry") for c in (7, 4)) + "; fraction observed: " + "; ".join(
        f"cl{c} {lvl} {a_mb[(a_mb.outcome == 'fraction_observed') & (a_mb.cluster == c) & (a_mb.level == lvl) & (a_mb.order == 'first') & (a_mb.block == 'Q+X')].mcfadden_gain.iloc[0]:.3f} -> "
        f"{mb[(mb.outcome == 'fraction_observed') & (mb.cluster == c) & (mb.level == lvl) & (mb.order == 'first') & (mb.block == 'Q+X')].gain.iloc[0]:.3f}"
        for lvl in ("chain", "entry") for c in (7, 4)) + ".")
    L.append("")
    lat = orr[orr.predictor.str.startswith("lattice")]
    L += ["Lattice terms (per 1 SD; odds ratio for complete / fraction observed, coefficient in z units for mobility):", "",
          "| outcome | cluster | level | term | estimate [95% CI] | q |", "|---|---|---|---|---|---|"]
    for r in lat.itertuples():
        est = (f"{r.estimate:.3g} [{r.ci_low:.3g}, {r.ci_high:.3g}]" if r.status == "fitted"
               else r.status.split(" (")[0])
        L.append(f"| {r.outcome} | cl{r.cluster} | {r.level} | {r.predictor.split('_')[1]} | {est} | "
                 f"{'' if pd.isna(r.q_value) else f'{r.q_value:.2g}'} |")
    L.append("")
    sig = orr[(orr.status == "fitted") & (orr.q_value < 0.05)].sort_values(["outcome", "cluster", "level", "q_value"])
    L.append("All effects with BH q < 0.05 (full model): " + "; ".join(
        f"{r.outcome} cl{r.cluster} {r.level} {r.label} ({r.block}) {r.estimate:.3g} [{r.ci_low:.3g}, {r.ci_high:.3g}], "
        f"q {r.q_value:.2g}" for r in sig.itertuples()) + ".")
    L.append("")
    dr = orr[orr.status != "fitted"]
    L.append("Dropped predictors (cluster/level): " + " | ".join(
        f"{o} cl{c} {lvl}: " + ", ".join(f"{r.predictor} [{r.status.replace('dropped: ', '').split(' (')[0]}]"
                                         for r in g.itertuples())
        for (o, c, lvl), g in dr.groupby(["outcome", "cluster", "level"], sort=False)) + ".")
    L.append("")
    ev = mb[(mb.order == "first") & (mb.block == "Q")]
    L.append("Rows / predictors / full R2: " + "; ".join(
        f"{r.outcome} cl{r.cluster} {r.level} n {r.n}"
        + (f" ({int(r.n_incomplete)} incomplete)" if pd.notna(r.n_incomplete) else "")
        + f", {r.n_predictors_full} predictors, R2 {r.r2_full:.3f}" for r in ev.itertuples()) + ".")
    L.append("")

    # verdict
    L += ["## Verdict (verdict.md)", "", (OUT / "verdict.md").read_text().strip(), "",
          "## 8. Medoid H12 clash test", "", (OUT / "clash_test_not_run.md").read_text().strip(), ""]

    a_or = rd("odds_ratios.csv", ROOT / "results/sv/SV01a")
    diffs = []
    for (o, c, lvl), d in orr[orr.outcome.isin(["complete", "fraction_observed"])].groupby(
            ["outcome", "cluster", "level"], sort=False):
        fb = set(d[d.status == "fitted"].predictor) - {"lattice_asu", "lattice_sym"}
        fa = set(a_or[(a_or.outcome == o) & (a_or.cluster == c) & (a_or.level == lvl) & (a_or.status == "fitted")].predictor)
        if fa != fb:
            diffs.append(f"{o} cl{c} {lvl} loses {', '.join(sorted(fa - fb))} to joint separation")
    mob_b = orr[(orr.outcome == "mobility") & (orr.cluster == 4) & (orr.level == "chain") & (orr.block == "B")
                & (orr.status == "fitted") & (orr.q_value < 0.05)]
    ent = mb[(mb.outcome == "mobility_no_tls") & (mb.level == "entry") & (mb.order == "first") & (mb.block == "Q")]
    cl7c = mb[(mb.outcome == "complete") & (mb.cluster == 7) & (mb.order == "first") & (mb.block == "Q")].set_index("level")
    g = lg.set_index(["cluster_group", "h12_status"])
    c4c, c4p = g.loc[("cl4", "complete")], g.loc[("cl4", "partial")]
    L += ["## Reading", "",
          "- Adding the lattice terms moves the completeness verdict from mixed (SV01a) to artifact. Apart from the "
          "lattice terms the fitted predictors are those of SV01a" + (f" (except: {'; '.join(diffs)})" if diffs else "")
          + ", so the rise of the Q+X gain is the lattice terms' share.",
          f"- For completeness this share is partly mechanical: contacts are counted on modelled atoms. In cl4, complete "
          f"chains have {c4c.mean_sym_atoms:.1f} mean symmetry-contact atoms on 528-550, partial chains "
          f"{c4p.mean_sym_atoms:.1f}; {c4c.mean_sym_atoms_546_550:.1f} of the former sit on 546-550, where partial chains "
          f"model {c4p.mean_modelled_heavy_atoms_546_550:.1f} heavy atoms on average (complete "
          f"{c4c.mean_modelled_heavy_atoms_546_550:.1f}). On 538-545, which both model, the means are "
          f"{c4c.mean_sym_atoms_538_545:.2f} vs {c4p.mean_sym_atoms_538_545:.2f}. The data cannot tell whether the "
          "contact holds the C-terminal end in place or exists only because the end is modelled; step 8 would "
          "separate the two but its condition (verdict mixed) is not met.",
          "- Mobility is z-scored on the residues that are modelled, so it is not tied to missing atoms in the same "
          "way: more symmetry-mate contact goes with lower local B relative to the core in both clusters at chain "
          "level (lattice table above), the direction expected if contacts restrain the segment. The B effects with "
          "q < 0.05 in cl4 mobility (chain level) are " + ", ".join(
              f"{r.label} ({int(r.n_with)} chains) {r.estimate:+.2f}" if pd.notna(r.n_with) else f"{r.label} {r.estimate:+.2f}"
              for r in mob_b.itertuples()) + ".",
          f"- Small samples: cl7 has {int(cl7c.loc['chain'].n_incomplete)} incomplete chains "
          f"({int(cl7c.loc['entry'].n_incomplete)} entries) for {int(cl7c.loc['chain'].n_predictors_full)} fitted "
          "predictors; the TLS-free mobility models at entry level fit " + " and ".join(
              f"{int(r.n_predictors_full)} predictors on {int(r.n)} cl{r.cluster} entries (R2 {r.r2_full:.2f})"
              for r in ent.itertuples()) + ", so their gains are optimistic.", ""]

    files = sorted(OUT.glob("*")) + sorted((ROOT / "figures/sv/SV01b").glob("*")) + \
        sorted((ROOT / "figures/panels/sv/SV01b").glob("*")) + \
        [ROOT / f"scripts/sv/{n}" for n in ("SV01b_lattice.py", "SV01b_models.py", "SV01b_clash.py",
                                            "SV01b_figures.py", "SV01b_summary.py")] + \
        [ROOT / "scripts/sv/common/lattice_contacts.py"]
    L += ["## Files", ""] + [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files if p.name != "SUMMARY.md"] + \
        ["- docs/EXPLORATION.md: appended 'Suggested, not done (SV01b, 2026-11-08)'",
         "- results/sv/SV01b/SUMMARY.md (this file)"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
