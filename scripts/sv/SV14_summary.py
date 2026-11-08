"""SV14 SUMMARY.md, verdict.md and env.txt, written from the saved outputs only (rule 1.7).

Inputs: results/sv/SV14/{pocket_per_chain.csv, references.json, correlations.csv, tests.csv, filling_by_group.csv,
numbers.json, e2_reference_pocket_lining.csv, apo_blocking_residues.csv}, data/sv/SV14/{chains.csv,
ligand_instances.csv, apo_check.csv, e2_centroid.json}.
Run: .venv/bin/python -B scripts/sv/SV14_summary.py
"""
from __future__ import annotations

import json
import platform
import sys
import time
from importlib.metadata import version
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/sv/SV14"
DATA = ROOT / "data/sv/SV14"


def f(x, nd=0):
    return f"{x:.{nd}f}"


def qtxt(x):
    return f"{x:.2g}" if x >= 1e-3 else f"{x:.1e}"


def main():
    ref = json.loads((RES / "references.json").read_text())
    n = json.loads((RES / "numbers.json").read_text())
    cor = pd.read_csv(RES / "correlations.csv")
    tests = pd.read_csv(RES / "tests.csv")
    bg = pd.read_csv(RES / "filling_by_group.csv").set_index(["grouping", "group"])
    lin = pd.read_csv(RES / "e2_reference_pocket_lining.csv")
    blk = pd.read_csv(RES / "apo_blocking_residues.csv")
    inst = pd.read_csv(DATA / "ligand_instances.csv")
    e2c = json.loads((DATA / "e2_centroid.json").read_text())
    L, P = ref["ligand_100pct"], ref["pocket_100pct"]
    apo, e2p = ref["apo_baseline"], ref["e2_bound_pocket"]
    passing = [e for e, v in apo["entries"].items() if v["passes"]]
    failing = [e for e, v in apo["entries"].items() if not v["passes"]]
    big = lin[lin.volume_far_from_ligand_A3 >= 100]
    sml = lin[lin.volume_far_from_ligand_A3 < 100]
    top_res = pd.Series([int(t.split(":")[0]) for s in big.lining_residues_far_part_top12 for t in s.split(";")[:8]])
    common = sorted(top_res.value_counts()[lambda s: s == len(big)].index) if len(big) else []

    def c(subset, level):
        r = cor[(cor.subset == subset) & (cor.level == level) & (cor.kind == "ligand vs pocket volume")].iloc[0]
        return r

    def sens(t, level):
        return cor[(cor.kind == f"buriedness sensitivity >= {t} vs >= 10 rays") & (cor.level == level)].iloc[0].rho

    T = tests.set_index(["comparison", "level"])
    tc_e, tc_c = T.loc[("cluster", "entry")], T.loc[("cluster", "chain")]
    ta_e, ta_c = T.loc[("action class", "entry")], T.loc[("action class", "chain")]
    g7, g4 = bg.loc[("cluster", "cl7")], bg.loc[("cluster", "cl4")]
    ag, se = bg.loc[("action_class", "agonist")], bg.loc[("action_class", "SERM")]
    efill = n["e2_ref_filling_values"]
    # ------------------------------------------------------------------ verdict
    verdict = [
        f"- Ligand and pocket volume go together over all chains (Spearman rho {f(c('all','ligand').rho,2)} ligand "
        f"level, {f(c('all','entry').rho,2)} entry level, q {qtxt(c('all','ligand').q_bh)} / "
        f"{qtxt(c('all','entry').q_bh)}) and within cl7 ({f(c('cl7','ligand').rho,2)} / {f(c('cl7','entry').rho,2)}), "
        f"not within cl4 ({f(c('cl4','ligand').rho,2)} / {f(c('cl4','entry').rho,2)}, q >= "
        f"{f(min(c('cl4','ligand').q_bh, c('cl4','entry').q_bh),2)}) nor within agonist or SERM (q >= "
        f"{f(min(c('class agonist','ligand').q_bh, c('class agonist','entry').q_bh, c('class SERM','ligand').q_bh, c('class SERM','entry').q_bh),2)}).",
        f"- Filling (median {f(n['filling_primary_median'],1)} %, IQR {f(n['filling_primary_q1'],1)}-"
        f"{f(n['filling_primary_q3'],1)}) does not differ cl7 vs cl4 (entry medians {f(tc_e.entry_median_cl7,1)} vs "
        f"{f(tc_e.entry_median_cl4,1)} %, q {f(tc_e.q_bh,2)}; chain bootstrap q {f(tc_c.q_bh,2)}) or agonist vs SERM "
        f"(q {f(ta_e.q_bh,2)} / {f(ta_c.q_bh,2)}).",
        f"- References: ligand 100 % = E2 {f(L['median_A3'])} A^3; pocket 100 % = wild-type E2 pocket "
        f"{f(P['median_A3'])} A^3 (fallback: {apo['entries_passing']} apo entries pass, < 5); the E2 pocket is "
        f"bimodal ({f(sml.pocket_volume_A3.min())}-{f(sml.pocket_volume_A3.max())} vs {f(big.pocket_volume_A3.min())}-"
        f"{f(big.pocket_volume_A3.max())} A^3): in {len(big)} of {len(lin)} chains it joins a buried region past "
        f"E353/R394.",
        f"- cl4 ligands are larger (median {f(g4.ligand_volume_A3_chain_median)} vs "
        f"{f(g7.ligand_volume_A3_chain_median)} A^3) and sit less completely inside the pocket (inside share "
        f"{f(g4.inside_share_chain_median,2)} vs {f(g7.inside_share_chain_median,2)}); filling stays similar because "
        f"their pockets are larger too ({f(g4.pocket_volume_A3_chain_median)} vs {f(g7.pocket_volume_A3_chain_median)} A^3).",
        "- These are crystal-structure pocket geometries only; they say nothing about cellular activity or receptor "
        "degradation.",
    ]
    (RES / "verdict.md").write_text("\n".join(verdict) + "\n")
    # ------------------------------------------------------------------ line 1
    line1 = (
        "Not done / substituted: (1) pocket 100 % = median pocket of wild-type EST chains (fallback, [DEFAULT] "
        f"step 5): only {apo['entries_passing']} apo entries pass ({', '.join(passing)}; {', '.join(failing)} fail "
        f"because a PEG (7 heavy atoms) sits in the pocket, centroid "
        f"{f(blk.centroid_to_e2_centroid_A.min(),1)}-{f(blk.centroid_to_e2_centroid_A.max(),1)} A from the E2 centroid), "
        f"fewer than 5; the apo baseline ({f(apo['median_of_entry_medians_A3'])} A^3) is reported but not used as "
        "100 %. (2) No pocket program was used: the grid method is coded in numpy/scipy (scripts/sv/SV14_grid.py) as "
        "the step describes; receptor atoms read from the restored mmCIFs with esr_analysis.read_structure.read_cif "
        "(highest-occupancy altloc), not from data/clean/*.pdb, which keep every altloc atom (same atom selection "
        "otherwise: own polymer only, no ligand, water, ion, peptide or other chain). (3) Kept as computed, flagged: "
        f"pocket_open is True for {n['pocket_open_chains']}/{n['chains']} chains (every component borders some "
        "non-buried cavity point), so the border count, not the flag, carries the information; the >= 10-ray pocket "
        f"reaches the [DEFAULT] box edge in {n['touches_box_b10_chains']} chains ({n['touches_box_b8_chains']} at >= 8, "
        f"{n['touches_box_b12_chains']} at >= 12), so those volumes are bounded by the box. (4) Declared before "
        "the outcomes where the run file is open: ligand instances with different alternate-location ids -> the "
        f"highest-occupancy alt id ({n['alt_instance_rule_chains']} chains); apo check over the asymmetric unit, any "
        "residue other than the chain's own polymer or water; E2 references = wild-type chains with EST alone, "
        "complete; rays marched on the 0.4 A grid (diagonals 14 steps = 9.70 A); ligand volume counted on the "
        "extended grid; several overlapping components -> largest overlap; class = SV12 action_class.csv row of "
        "(entry, analysed CCD), class analyses on primary_set rows; entry cluster group = >= 80 % rule; chain level "
        "= bootstrap over entries (10000, seed 42); BH over all 14 tested p-values (10 correlations + 4 filling "
        "tests). (5) Supporting files not named in the run file: filling_by_group.csv, filling_entry_level.csv, "
        "numbers.json, e2_reference_pocket_lining.csv, apo_blocking_residues.csv (results/sv/SV14/), data/sv/SV14/*; "
        "scripts/sv/common/class_style.py created (did not exist). Nothing installed.")
    # ------------------------------------------------------------------ body
    def cor_row(subset):
        a, b = c(subset, "ligand"), c(subset, "entry")
        return (f"| {subset} | {int(a.n_units)} | {f(a.rho,2)} | {qtxt(a.q_bh)} | {int(b.n_units)} | {f(b.rho,2)} | "
                f"{qtxt(b.q_bh)} |")

    body = [line1, "", "# SV14 - Ligand size vs pocket size, and pocket filling", "", "**Verdict**", ""] + verdict + [
        "", "## Inputs (rule 1.10)", "",
        "- SV00 card: results/sv/SV00/chain_table.csv (788 chains; cluster, construct_mutations = wild type 'none', "
        "ligand_chemotype), restored mmCIFs data/raw/, data/sv/lbd_extract/receptor_chains.csv and residue_mapping.csv.",
        "- Q18 card: results/chain_ligand_of_interest.csv; SV11 card: results/sv/SV11/site_residues.csv (ligand "
        "pocket, in_site: 18 residues); SV03 card: data/sv/SV03/core_v3_transforms.csv (core_frame.py); SV12 card: "
        "results/sv/SV12/action_class.csv, class_coverage.csv, CCD records data/sv/rcsb_chemcomp/ (atom_count_heavy).",
        "", "## Steps 1-4 - grid, pocket, ligand, filling", "",
        "- Grid 0.4 A (voxel 0.064 A^3) over the box of the 18 ligand-pocket residues + 6 A in the core v3 frame; "
        "Bondi radii, no hydrogens; probe 1.4 A; buried = >= 10 of 14 rays hit a receptor vdW sphere within 10 A; "
        "pocket = 6-connected buried cavity component overlapping the ligand (apo: nearest to the E2 centroid).",
        f"- Atoms: receptor heavy atoms per chain median {int(pd.read_csv(DATA / 'chains.csv').n_rec_heavy.median())}; "
        f"ligand instances {len(inst)} in {n['chains_with_ligand']} chains, {int(inst.kept.sum())} kept "
        f"({n['multi_instance_chains']} chains with several Q18 instances).",
        f"- ligand_incomplete (fewer modelled heavy atoms than the CCD): {n['ligand_incomplete_chains']} chains in "
        f"{n['ligand_incomplete_entries']} entries, out of the primary analysis; no ligand: 7 chains. Primary set: "
        f"{n['primary_chains']} chains, {n['primary_entries']} entries, {n['primary_ligands']} ligands (cl7 "
        f"{n['primary_by_cluster'].get('cl7')}, cl4 {n['primary_by_cluster'].get('cl4')}, minor "
        f"{n['primary_by_cluster'].get('minor')}).",
        f"- Primary set: pocket {f(n['pocket_volume_primary_median'])} A^3 (IQR {f(n['pocket_volume_primary_q1'])}-"
        f"{f(n['pocket_volume_primary_q3'])}, range {f(n['pocket_volume_primary_min'])}-"
        f"{f(n['pocket_volume_primary_max'])}); ligand {f(n['ligand_volume_primary_median'])} A^3; filling "
        f"{f(n['filling_primary_median'],1)} % (IQR {f(n['filling_primary_q1'],1)}-{f(n['filling_primary_q3'],1)}); "
        f"inside share median {f(n['inside_share_primary_median'],2)}, < 0.90 in {n['inside_share_lt_0_9_chains']} "
        f"chains; ligand overlaps > 1 buried component in {n['overlap_gt1_component_chains']} chains; ligand vdW "
        f"outside the box in {n['ligand_outside_box_chains']} chains (counted, as outside the pocket).",
        f"- pocket_open: {n['pocket_open_chains']}/{n['chains']}; border points median {f(n['border_points_median'])} "
        f"(IQR {f(n['border_points_q1'])}-{f(n['border_points_q3'])}); E2 reference chains "
        f"{min(n['e2_ref_border_points'])}-{max(n['e2_ref_border_points'])}.",
        f"- Box edge (>= 10 rays) in the primary set: {n['touches_box_b10_primary']} chains (cl7 "
        f"{n['touches_box_b10_by_cluster_primary'].get('cl7', 0)}, cl4 {n['touches_box_b10_by_cluster_primary'].get('cl4', 0)}, "
        f"minor {n['touches_box_b10_by_cluster_primary'].get('minor', 0)}).",
        f"- Buriedness sensitivity (primary set): Spearman rho with the >= 10-ray volume: >= 8 rays {f(sens(8,'chain'),2)} "
        f"chain / {f(sens(8,'entry'),2)} entry; >= 12 rays {f(sens(12,'chain'),2)} / {f(sens(12,'entry'),2)}.",
        "", "## Step 5 - references (references.json)", "",
        f"- Ligand 100 % = {f(L['median_A3'],1)} A^3 (IQR {f(L['q1_A3'],1)}-{f(L['q3_A3'],1)}), median EST volume over "
        f"{L['n_chains']} wild-type EST chains ({e2p['n_entries']} entries).",
        f"- Apo entries: {apo['candidate_entries']} entries ({apo['candidate_chains']} chains) with ligand_chemotype "
        f"none; {apo['entries_passing']} pass ({apo['chains_passing']} chains: {', '.join(passing)}); "
        f"{', '.join(failing)}: PEG in the pocket (apo_blocking_residues.csv). Apo baseline = "
        f"{f(apo['median_of_entry_medians_A3'],1)} A^3 (IQR of entry medians {f(apo['q1_A3'],1)}-{f(apo['q3_A3'],1)}; "
        f"entry medians " + ", ".join(f"{e} {f(v['median_pocket_A3'])}" for e, v in apo["entries"].items() if v["passes"])
        + ").",
        f"- Pocket 100 % = {P['source']}: {f(P['median_A3'],1)} A^3 (IQR {f(P['q1_A3'],1)}-{f(P['q3_A3'],1)}); E2-bound "
        f"pocket (second reference) = the same {e2p['n_chains']} chains. The apo baseline is "
        f"{f(100 * apo['median_of_entry_medians_A3'] / P['median_A3'])} % of it.",
        f"- The E2 reference pockets split: {len(sml)} chains {f(sml.pocket_volume_A3.min())}-"
        f"{f(sml.pocket_volume_A3.max())} A^3 (<= {f(sml.volume_far_from_ligand_A3.max())} A^3 farther than 4 A from "
        f"E2) and {len(big)} chains {f(big.pocket_volume_A3.min())}-{f(big.pocket_volume_A3.max())} A^3, where "
        f"{f(big.volume_far_from_ligand_A3.min())}-{f(big.volume_far_from_ligand_A3.max())} A^3 lie farther than 4 A "
        f"from E2, lined by residues {', '.join(map(str, common))} (among the 8 most frequent in every one of those "
        f"chains; e2_reference_pocket_lining.csv). Same entry can fall in both groups (1ERE A, B, F vs C, D, E). "
        f"E2 reference filling {f(min(efill),1)}-{f(max(efill),1)} %.",
        f"- E2 reference centroid: median of {e2c['n_chains']} EST centroids, all within "
        f"{f(e2c['max_distance_of_chain_centroid_to_median_A'],2)} A of it.",
        "", "## Step 6 - correlations (correlations.csv; Spearman, ligand = median per analysed CCD)", "",
        "| subset | ligands | rho | q | entries | rho | q |", "|---|---|---|---|---|---|---|",
        cor_row("all"), cor_row("cl7"), cor_row("cl4"), cor_row("class agonist"), cor_row("class SERM"), "",
        f"- Class gate (SV12 class_coverage.csv, >= 5 entries): {', '.join(ref['class_gate']['analysed'])} "
        f"(entries in the SV14 primary set: " + ", ".join(f"{k} {v}" for k, v in
                                                       ref['class_gate']['entries_in_sv14_primary_set'].items())
        + "); SERD (4 entries), antagonist, partial agonist, other not testable.",
        "", "## Step 6 - filling tests (tests.csv, BH with the correlations)", "",
        "| comparison | level | medians | effect | p | q |", "|---|---|---|---|---|---|",
        f"| cl7 vs cl4 | entry ({int(tc_e.entries_cl7)} vs {int(tc_e.entries_cl4)}) | {f(tc_e.entry_median_cl7,1)} vs "
        f"{f(tc_e.entry_median_cl4,1)} % | rank-biserial {f(tc_e.effect,2)} | {qtxt(tc_e.p)} | {qtxt(tc_e.q_bh)} |",
        f"| cl7 vs cl4 | chain ({int(tc_c.chains_cl7)} vs {int(tc_c.chains_cl4)}) | {f(tc_c.chain_median_cl7,1)} vs "
        f"{f(tc_c.chain_median_cl4,1)} % | {f(tc_c.effect,1)} pp [{f(tc_c.ci_low,1)}, {f(tc_c.ci_high,1)}] | "
        f"{qtxt(tc_c.p)} | {qtxt(tc_c.q_bh)} |",
        f"| agonist vs SERM | entry ({int(ta_e.entries_agonist)} vs {int(ta_e.entries_SERM)}) | "
        f"{f(ta_e.entry_median_agonist,1)} vs {f(ta_e.entry_median_SERM,1)} % | rank-biserial {f(ta_e.effect,2)} | "
        f"{qtxt(ta_e.p)} | {qtxt(ta_e.q_bh)} |",
        f"| agonist vs SERM | chain ({int(ta_c.chains_agonist)} vs {int(ta_c.chains_SERM)}) | "
        f"{f(ta_c.chain_median_agonist,1)} vs {f(ta_c.chain_median_SERM,1)} % | {f(ta_c.effect,1)} pp "
        f"[{f(ta_c.ci_low,1)}, {f(ta_c.ci_high,1)}] | {qtxt(ta_c.p)} | {qtxt(ta_c.q_bh)} |", "",
        f"- By group (filling_by_group.csv, chain medians): cl7 pocket {f(g7.pocket_volume_A3_chain_median)} A^3, "
        f"ligand {f(g7.ligand_volume_A3_chain_median)} A^3, inside share {f(g7.inside_share_chain_median,2)}; cl4 "
        f"{f(g4.pocket_volume_A3_chain_median)}, {f(g4.ligand_volume_A3_chain_median)}, "
        f"{f(g4.inside_share_chain_median,2)}; agonist {f(ag.pocket_volume_A3_chain_median)}, "
        f"{f(ag.ligand_volume_A3_chain_median)}, {f(ag.inside_share_chain_median,2)}; SERM "
        f"{f(se.pocket_volume_A3_chain_median)}, {f(se.ligand_volume_A3_chain_median)}, "
        f"{f(se.inside_share_chain_median,2)}.",
        "", "## Files", "",
        "- results/sv/SV14/: pocket_per_chain.csv, references.json, correlations.csv, tests.csv, verdict.md, "
        "SUMMARY.md, env.txt, filling_by_group.csv, filling_entry_level.csv, numbers.json, "
        "e2_reference_pocket_lining.csv, apo_blocking_residues.csv",
        "- figures/sv/SV14/: SV14_pocket_vs_ligand.png/.svg (main), SV14_filling.png/.svg; panels "
        "figures/panels/sv/SV14/ (3 PNG)",
        "- data/sv/SV14/: atoms_core.npz, chains.csv, ligand_instances.csv, apo_check.csv, e2_centroid.json, "
        "grid_per_chain.csv, extract_stdout.txt, grid_stdout.txt",
        "- scripts/sv/: SV14_extract.py, SV14_grid.py, SV14_analysis.py, SV14_e2_lining.py, SV14_apo_blockers.py, "
        "SV14_figures.py, SV14_summary.py; scripts/sv/common/class_style.py",
    ]
    (RES / "SUMMARY.md").write_text("\n".join(body) + "\n")
    env = [f"# SV14 environment, {time.strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for p in ("numpy", "scipy", "pandas", "statsmodels", "matplotlib", "biopython"):
        env.append(f"{p} {version(p)}")
    env += ["grid method: own numpy/scipy code (scipy.ndimage.distance_transform_edt, label, binary_dilation); "
            "no external pocket program", "network: none (CCD heavy-atom counts from the SV12 cache "
            "data/sv/rcsb_chemcomp/)", "nothing installed, upgraded or downgraded by SV14"]
    (RES / "env.txt").write_text("\n".join(env) + "\n")
    print("\n".join(body))


if __name__ == "__main__":
    main()
