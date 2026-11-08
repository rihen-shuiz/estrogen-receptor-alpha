"""SV01c: write results/sv/SV01c/SUMMARY.md (rule 1.7) from the saved SV01c (and SV01b) files only.
Run: .venv/bin/python -B scripts/sv/SV01c_summary.py (after SV01c_window.py and, if run, SV01c_clash.py)
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV01c"
SV01B = ROOT / "results/sv/SV01b"


def rd(name, d=OUT):
    return pd.read_csv(d / name, low_memory=False)


def size(p):
    s = p.stat().st_size
    return f"{s / 1e3:.1f} kB" if s < 1e6 else f"{s / 1e6:.2f} MB"


def f3(x):
    return f"{x:.3f}"


def g3(x):
    return f"{x:.3g}"


def main():
    n = json.loads((OUT / "numbers.json").read_text())
    rep = json.loads((OUT / "reproduce_SV01b.json").read_text())
    mb, orr, sbs, st = rd("model_blocks.csv"), rd("odds_ratios.csv"), rd("side_by_side.csv"), rd("status_transitions.csv")
    lg = rd("lattice_contacts_by_group.csv")
    w0, w1 = n["window_proposed"]
    p0, p1 = n["window_primary"]
    wt, pt = f"{w0}_{w1}", f"{p0}_{p1}"
    step8 = n["run_step8"]
    if step8:
        cs, ck, ct = rd("clash_summary.csv"), json.loads((OUT / "clash_checks.json").read_text()), rd("clash_test.csv")
        fix = ck["distance_fix"]
    ov, lv, bov, blv = n["outcome_verdicts"], n["level_verdicts"], n["sv01b_outcome_verdicts"], n["sv01b_level_verdicts"]

    line1 = (
        f"Window propagation, fixed before fitting (docstring of scripts/sv/SV01c_window.py): H12 complete = all "
        f"{w1 - w0 + 1} C-alpha of {w0}-{w1} observed (fraction observed over {w0}-{w1}; status complete / partial / not "
        f"modelled re-derived on {w0}-{w1}); mobility = mean C-alpha B z-score over the observed residues of "
        f"{n['mobility_range'][0]}-{n['mobility_range'][1]} (SV01b: 528-{p1}; loop 528-537, core v3 reference, B and "
        f"altloc rules unchanged); lattice terms keep {n['lattice_range'][0]}-{n['lattice_range'][1]}, the SV01b counts "
        "reused unchanged (step 6 defines them on 528-550, not on the H12 window). Not done / substituted: no figure "
        "(step 9 asks for the verdict side by side, given as tables); "
        + (f"step 8 run because the SV01c verdict is mixed, with scripts/sv/SV01b_clash.py unchanged except that the "
           f"distance it reports is taken to the nearest lattice image (the shared lattice_contacts.lattice_hits reports "
           f"the distance to the neighbour's stored image, which is wrong for {fix['n_hits_stored_image_dist_gt_cutoff']}"
           f"/{fix['n_hits']} hits found across a lattice translation, up to {fix['max_stored_image_dist_A']:.1f} A; "
           f"hit sets identical, so clash counts and 'blocked' are those of the shared code; the module is not edited); "
           if step8 else "step 8 not run (verdict not mixed); ")
        + f"everything else inherited unchanged from SV01b (blocks Q/X/B with lattice terms in X, median imputation + "
        f"'not reported' indicators, drop rule, McFadden / OLS R2 gains, block tests, one BH family over the "
        f"{n['n_tests_bh']} SV01c p-values, [DEFAULT] verdict rule with the SV01b step rule).")
    L = [line1, "", f"# SV01c - H12 completeness and mobility on the data-derived H12 window {w0}-{w1}", ""]
    L += [f"**SV01c verdict ({w0}-{w1}): {n['step_verdict']}; primary SV01b verdict ({p0}-{p1}): "
          f"{n['sv01b_step_verdict']}.** Complete H12: {bov['complete']} -> {ov['complete']} (chain "
          f"{blv['complete|chain']} -> {lv['complete|chain']}, entry {blv['complete|entry']} -> {lv['complete|entry']}); "
          f"mobility: {bov['mobility']} -> {ov['mobility']} (chain {blv['mobility|chain']} -> {lv['mobility|chain']}, "
          f"entry {blv['mobility|entry']} -> {lv['mobility|entry']}). Sensitivity: fraction observed "
          f"{bov['fraction_observed']} -> {ov['fraction_observed']}, mobility without TLS-refined entries "
          f"{bov['mobility_no_tls']} -> {ov['mobility_no_tls']}. "
          + ("Step 8 run (verdict mixed)." if step8 else "Step 8 not run (verdict not mixed)."), ""]

    L += ["## Inputs", "",
          f"- SV02 (DONE): data-derived H12 window {w0}-{w1} = results/sv/SV02/helix_boundaries.csv row H12 "
          "(global consensus; proposal, not adopted; results/sv/SV02/h12_window.md).",
          "- SV01b (DONE): scripts/sv/SV01b_models.py imported with SV01a_completeness.H12 and SV01b_models.MOB set to "
          "the window (no file edited); lattice terms from results/sv/SV01b/lattice_contacts_per_chain.csv "
          "(scripts/sv/common/lattice_contacts.py output, not recomputed).",
          f"- Check (reproduce_SV01b.json): the same code on {p0}-{p1} reproduces results/sv/SV01b: "
          f"{rep['model_block_rows']}/{rep['model_block_rows_sv01b']} block rows, max |d gain| "
          f"{rep['max_abs_diff_gain']:.1e}, max |d q| {rep['max_abs_diff_block_q']:.1e}; {rep['odds_rows']}/"
          f"{rep['odds_rows_sv01b']} effect rows, status equal {rep['status_equal']}, max |d coef| "
          f"{rep['max_abs_diff_coef']:.1e}; BH tests {rep['n_tests_bh']} = {rep['n_tests_bh_sv01b']}; mobility z "
          f"for {rep['mobility_rows']} chains, max |d| {rep['max_abs_diff_mobility_z']:.1e}. Only the window differs.",
          f"- C-alpha B cross-check (backbone.npz vs residue mapping): {n['b_crosscheck']['n_compared']} C-alpha, max "
          f"|d| {n['b_crosscheck']['max_abs_diff']}.", ""]

    L += [f"## Outcomes on the two windows (chain_outcomes.csv, status_transitions.csv)", "",
          f"| group | level | {p0}-{p1} status | {w0}-{w1} status | n |", "|---|---|---|---|---|"]
    for r in st[st.cluster_group.isin(["cl7", "cl4", "minor", "all"])].itertuples():
        L.append(f"| {r.cluster_group} | {r.level} | {getattr(r, f'h12_status_{pt}')} | "
                 f"{getattr(r, f'h12_status_{wt}')} | {r.n} |")
    ic, ie = n["n_incomplete_chains"], n["n_incomplete_entries"]
    L += ["", f"Incomplete H12 (model events) {p0}-{p1} -> {w0}-{w1}: cl7 {ic[pt]['cl7']} -> {ic[wt]['cl7']} chains "
          f"({ie[pt]['cl7']} -> {ie[wt]['cl7']} entries at entry level); cl4 {ic[pt]['cl4']} -> {ic[wt]['cl4']} chains "
          f"({ie[pt]['cl4']} -> {ie[wt]['cl4']} entries). Mobility defined for cl7 {n['mobility_defined']['cl7']}, cl4 "
          f"{n['mobility_defined']['cl4']} chains (no TLS: {n['mobility_defined_no_tls']['cl7']}, "
          f"{n['mobility_defined_no_tls']['cl4']}); undefined: {', '.join(n['mobility_undefined_chains'])}.", ""]

    L += ["## Verdict side by side (side_by_side.csv; gain = McFadden R2 for complete / fraction observed, OLS R2 for "
          "mobility)", "",
          f"| outcome | level | cluster | order | Q+X {p0}-{p1} | B {p0}-{p1} | Q+X {w0}-{w1} | B {w0}-{w1} | level "
          f"verdict {p0}-{p1} -> {w0}-{w1} |", "|---|---|---|---|---|---|---|---|---|"]
    for r in sbs.itertuples():
        L.append(f"| {r.outcome} | {r.level} | cl{r.cluster} | {r.order} | {f3(getattr(r, f'gain_QX_{pt}'))} | "
                 f"{f3(getattr(r, f'gain_B_{pt}'))} | {f3(getattr(r, f'gain_QX_{wt}'))} | "
                 f"{f3(getattr(r, f'gain_B_{wt}'))} | {getattr(r, f'level_verdict_{pt}')} -> "
                 f"{getattr(r, f'level_verdict_{wt}')} |")
    L.append("")

    L += [f"## Models on {w0}-{w1} (model_blocks.csv, odds_ratios.csv)", "",
          "| outcome | cluster | level | n | events | fitted predictors | full R2 | test Q+X (last) | q Q+X | "
          "test B (last) | q B |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    last = mb[mb.order == "last"]
    for (o, c, lvl), m in last.groupby(["outcome", "cluster", "level"], sort=False):
        qx, bb = m[m.block == "Q+X"].iloc[0], m[m.block == "B"].iloc[0]

        def tst(x):
            if pd.isna(x.statistic):
                return "-"
            nm = "Wald" if "Wald" in x.test else "LR" if x.test.startswith("LR") else "F"
            return f"{nm} {x.statistic:.1f} (df {int(x.df)}), p {g3(x.p_value)}"
        ev = "-" if pd.isna(qx.n_incomplete) else int(qx.n_incomplete)
        L.append(f"| {o} | cl{c} | {lvl} | {qx.n} | {ev} | {qx.n_predictors_full} | {f3(qx.r2_full)} | {tst(qx)} | "
                 f"{g3(qx.q_value) if pd.notna(qx.q_value) else '-'} | {tst(bb)} | "
                 f"{g3(bb.q_value) if pd.notna(bb.q_value) else '-'} |")
    L += ["", "Dropped predictors, complete outcome (cluster / level: predictor [reason]): " + " | ".join(
        f"cl{c} {lvl}: " + ", ".join(f"{r.predictor} [{r.status.replace('dropped: ', '').split(' (')[0]}]"
                                     for r in g.itertuples())
        for (c, lvl), g in orr[(orr.outcome == "complete") & (orr.status != "fitted")].groupby(
            ["cluster", "level"], sort=False)) + "."]
    sig = orr[orr.q_value < 0.05].sort_values(["outcome", "cluster", "level", "q_value"])
    L += ["", "Effects with BH q < 0.05 (full model; odds ratio per 1 SD or present vs absent for complete / fraction "
          "observed, coefficient in z units for mobility): " + "; ".join(
              f"{r.outcome} cl{r.cluster} {r.level} {r.label} ({r.block}) {g3(r.estimate)} [{g3(r.ci_low)}, "
              f"{g3(r.ci_high)}], q {g3(r.q_value)}" for r in sig.itertuples()) + ".", ""]
    L += ["Lattice contacts of 528-550 by H12 status on the new window (lattice_contacts_by_group.csv; mean heavy "
          "atoms in contact):", "",
          "| group | H12 (538-544) | chains | entries | any ASU | any sym | mean ASU atoms | mean sym atoms | modelled "
          "heavy atoms 528-550 |", "|---|---|---|---|---|---|---|---|---|"]
    for r in lg.itertuples():
        L.append(f"| {r.cluster_group} | {r.h12_status} | {r.n_chains} | {r.n_entries} | {r.chains_asu_contact} | "
                 f"{r.chains_sym_contact} | {r.mean_asu_atoms:.2f} | {r.mean_sym_atoms:.2f} | "
                 f"{r.mean_modelled_heavy_atoms_528_550:.1f} |")
    L.append("")

    if step8:
        L += [f"## 8. Medoid H12 clash test on {w0}-{w1} (clash_test.csv, clash_summary.csv, clash_checks.json)", "",
              f"- Targets: {ck['n_targets']} chains partial or not modelled on {w0}-{w1} (all clusters: " + ", ".join(
                  f"{d['cluster_group']} {d['h12_status']} {d['n']}" for d in ck["targets_by_group_status"]) +
              f"); moved atoms = medoid residues {w0}-{w1} ({', '.join(f'{k} {v}' for k, v in ck['medoid_h12_atoms'].items())} "
              f"heavy atoms); clash = {ck['clash_cutoff_A']} A; blocked = >= 1 backbone atom in clash.",
              f"- Core fit: >= {ck['core_ca_fit_min']} core v3 C-alpha per fit, RMSD median {ck['core_rmsd_median']:.2f} A, "
              f"max {ck['core_rmsd_max']:.2f} A; direct fit vs SV00 gate transforms on the medoid H12 C-alpha: median "
              f"{ck['frame_check_direct_vs_gate_h12_ca_max_dev_A']['median']:.2f} A, max "
              f"{ck['frame_check_direct_vs_gate_h12_ca_max_dev_A']['max']:.2f} A (n "
              f"{ck['frame_check_direct_vs_gate_h12_ca_max_dev_A']['n']}).", "",
              "| group | H12 | chains | cl7 position blocked | cl4 position blocked | only cl7 | only cl4 | both | "
              "neither |", "|---|---|---|---|---|---|---|---|---|"]
        for (gname, sname), d in cs.groupby(["cluster_group", "h12_status"], sort=False):
            k = d.set_index("measure").k
            L.append(f"| {gname} | {sname} | {int(d.n_chains.iloc[0])} | {int(k['cl7 medoid H12 blocked'])} | "
                     f"{int(k['cl4 medoid H12 blocked'])} | {int(k['only cl7 position blocked'])} | "
                     f"{int(k['only cl4 position blocked'])} | {int(k['both positions blocked'])} | "
                     f"{int(k['neither position blocked'])} |")
        bl = ct[ct.blocked]
        L += ["", "Blocked chain x medoid rows: " + "; ".join(
            f"{r.chain_key} ({r.cluster_group}, {r.h12_status}) {r.medoid_cluster} medoid: {r.clash_backbone_atoms} "
            f"backbone / {r.clash_atoms} heavy atoms in clash, min {r.min_distance_A:.2f} A, partner "
            f"{r.clash_partners} ({r.clash_partner_kinds}, {'ASU' if r.clash_atoms_asu else 'symmetry'})"
            for r in bl.itertuples()) + ".", ""]

    # reading (numbers from the files above)
    oc7 = orr[(orr.outcome == "complete") & (orr.cluster == 7) & (orr.level == "chain")].set_index("predictor")
    lha = oc7.loc["ligand_heavy_atoms"]
    trans = st[(st.level == "chain")].set_index(["cluster_group", f"h12_status_{pt}", f"h12_status_{wt}"]).n
    lgi = lg.set_index(["cluster_group", "h12_status"])
    npf = mb[(mb.outcome == "complete") & (mb.order == "first") & (mb.block == "Q+X")].n_predictors_full
    sm_ = sbs[sbs.outcome.isin(["mobility", "mobility_no_tls"])]
    dmob = max((sm_[f"gain_QX_{wt}"] - sm_[f"gain_QX_{pt}"]).abs().max(),
               (sm_[f"gain_B_{wt}"] - sm_[f"gain_B_{pt}"]).abs().max())
    L += ["## Reading", "",
          f"- On {w0}-{w1}, {int(trans[('all', 'partial', 'complete')])} chains that are partial on {p0}-{p1} become "
          f"complete (cl4 {int(trans[('cl4', 'partial', 'complete')])}, cl7 {int(trans[('cl7', 'partial', 'complete')])}, "
          f"minor {int(trans[('minor', 'partial', 'complete')])}); what stays incomplete is mostly H12 not modelled at "
          f"all (cl7 {int(trans[('cl7', 'not_modelled', 'not_modelled')])} of {ic[wt]['cl7']}, cl4 "
          f"{int(trans[('cl4', 'not_modelled', 'not_modelled')])} of {ic[wt]['cl4']} incomplete chains). The "
          "completeness outcome on the new window is therefore a different, much rarer event than on the primary "
          "window.",
          f"- The completeness verdict changes from artifact to mixed at chain level only: in cl7 the B gain exceeds the "
          f"Q+X gain in both orders, in cl4 the reverse. In cl7 (chain level) both lattice terms and space group P21 are "
          f"dropped by joint separation (not-modelled cl7 chains have no lattice contact on 528-550: "
          f"{int(lgi.loc[('cl7', 'not_modelled')].chains_sym_contact)} of "
          f"{int(lgi.loc[('cl7', 'not_modelled')].n_chains)} with a symmetry contact), so block X loses the terms that "
          f"carried its SV01b gain; the cl7 B block is then led by ligand heavy atoms, OR {g3(lha.estimate)} per 1 SD "
          f"[{g3(lha.ci_low)}, {g3(lha.ci_high)}], q {g3(lha.q_value)}.",
          f"- The events are few: cl7 {ic[wt]['cl7']} incomplete chains ({ie[wt]['cl7']} entries at entry level), "
          f"cl4 {ic[wt]['cl4']} chains ({ie[wt]['cl4']} entries), against {npf.min()}-{npf.max()} fitted predictors "
          f"(table above), so the gains on the new window are unstable and the entry-level 'artifact' rests on "
          f"{ie[wt]['cl7']} (cl7) and {ie[wt]['cl4']} (cl4) incomplete entries. At entry level every incomplete row "
          f"is H12 not modelled, so fraction observed gives the same fit as the binary outcome there.",
          f"- The mobility verdict does not depend on the window: it stays artifact at both levels and without "
          f"TLS-refined entries; dropping {w1 + 1}-{p1} from the mobility mean changes the Q+X and B gains by at most "
          f"{dmob:.3f} (side-by-side table)."]
    if step8:
        k_all = cs[(cs.cluster_group == "all")].set_index("measure").k
        bl_sym_rec = bool((bl.clash_atoms_asu == 0).all() and (bl.clash_partner_kinds == "receptor").all())
        L.append(f"- Step 8: of {ck['n_targets']} partial / not-modelled chains on {w0}-{w1}, "
                 f"{int(k_all['both positions blocked'])} have both medoid H12 positions blocked "
                 f"({'by symmetry copies of receptor chains' if bl_sym_rec else 'see the blocked rows above'}), "
                 f"{int(k_all['only cl7 position blocked'])} only the cl7 position, "
                 f"{int(k_all['only cl4 position blocked'])} only the cl4 position, and "
                 f"{int(k_all['neither position blocked'])} neither. The crystal lattice leaves room for either H12 "
                 "position in most chains whose H12 is missing on the new window, and where it blocks, it blocks both, so "
                 "it does not favour one position over the other.")
    L.append("- Wording: crystal structures show H12 order and position (modelled or not, local B), not receptor "
             "activity.")
    L.append("")

    files = ["verdict.md", "side_by_side.csv", "verdict_table.csv", "model_blocks.csv", "odds_ratios.csv",
             "chain_outcomes.csv", "status_transitions.csv", "mobility_per_chain.csv", "lattice_contacts_by_group.csv",
             "predictor_scaling.csv", "model_data.csv", "reproduce_SV01b.json", "numbers.json", "env.txt"]
    if step8:
        files += ["clash_test.csv", "clash_summary.csv", "clash_checks.json"]
    L += ["## Verdict (verdict.md)", ""] + (OUT / "verdict.md").read_text().splitlines() + [""]
    L += ["## Files", ""] + [f"- results/sv/SV01c/{f} ({size(OUT / f)})" for f in files]
    for s in ("SV01c_window.py", "SV01c_clash.py", "SV01c_summary.py"):
        p = ROOT / "scripts/sv" / s
        if p.exists():
            L.append(f"- scripts/sv/{s} ({size(p)})")
    L += ["- docs/EXPLORATION.md: appended 'Suggested, not done (SV01c, 2026-11-08)'",
          "- results/sv/SV01c/SUMMARY.md (this file)"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
