"""SV15 SUMMARY.md, verdict.md and env.txt, written from the saved outputs only (rule 1.7).

Inputs: results/sv/SV15/{affinity_per_ligand.csv, assay_spread.csv, tests.csv, paired_assays.csv, entry_level.csv,
tests_first_run_type_units.csv, numbers.json}; figure paths under figures/sv/SV15/ and figures/panels/sv/SV15/.
Run: .venv/bin/python -B scripts/sv/SV15_summary.py
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
RES = ROOT / "results/sv/SV15"
REL_KEY = "standard_relation '='"


def qtxt(x):
    return f"{x:.2g}" if x >= 1e-3 else f"{x:.1e}"


def fx(x):
    """2 decimals without a negative zero."""
    v = f"{x:.2f}"
    return "0.00" if v == "-0.00" else v


def size(p: Path) -> str:
    return f"{p.stat().st_size / 1024:.1f} KB"


def row(tests, comparison, level, set_name="primary"):
    return tests[(tests.set == set_name) & (tests.comparison == comparison) & (tests.level == level)].iloc[0]


def main():
    n = json.loads((RES / "numbers.json").read_text())
    t = pd.read_csv(RES / "affinity_per_ligand.csv")
    sp = pd.read_csv(RES / "assay_spread.csv")
    tests = pd.read_csv(RES / "tests.csv")
    first = pd.read_csv(RES / "tests_first_run_type_units.csv")
    pa = pd.read_csv(RES / "paired_assays.csv")
    for c in ("primary_set", "is_e2_reference", "in_tests_primary", "in_tests_sensitivity"):
        t[c] = t[c].astype(bool)
    rel = t[t.in_tests_primary].sort_values("dpx_median", ascending=False)
    cas = {s["step"]: s for s in n["cleaning_cascade"]}
    g = n["class_gate_primary"]
    cl_l, cl_e = row(tests, "dpX by SV10 ligand group", "ligand"), row(tests, "dpX by SV10 ligand group", "entry")
    fi_l, fi_e = row(tests, "dpX vs filling % (SV14)", "ligand"), row(tests, "dpX vs filling % (SV14)", "entry")
    vo_l, vo_e = (row(tests, "dpX vs ligand volume (SV14)", "ligand"),
                  row(tests, "dpX vs ligand volume (SV14)", "entry"))
    sens = "sensitivity (connectivity_only included)"
    s_cl_l, s_cl_e = row(tests, "dpX by SV10 ligand group", "ligand", sens), row(tests, "dpX by SV10 ligand group",
                                                                                   "entry", sens)
    s_fi_l, s_fi_e = row(tests, "dpX vs filling % (SV14)", "ligand", sens), row(tests, "dpX vs filling % (SV14)",
                                                                                  "entry", sens)
    s_vo_l, s_vo_e = (row(tests, "dpX vs ligand volume (SV14)", "ligand", sens),
                      row(tests, "dpX vs ligand volume (SV14)", "entry", sens))
    f_e = first[(first.comparison == "dpX by SV10 ligand group") & (first.level == "entry")].iloc[0]
    f_l = first[(first.comparison == "dpX by SV10 ligand group") & (first.level == "ligand")].iloc[0]
    f_ve = first[(first.comparison == "dpX vs ligand volume (SV14)") & (first.level == "entry")].iloc[0]
    tp_ = tests[(tests.set == "primary") & tests.tested].set_index(["comparison", "level"])
    ts_ = tests[(tests.set == sens) & tests.tested].set_index(["comparison", "level"])
    common = tp_.index.intersection(ts_.index)
    diff = [k for k in common if (tp_.loc[k, "q_bh"] < 0.05) != (ts_.loc[k, "q_bh"] < 0.05)]
    same_calls = (f"q < 0.05 calls equal to the primary set in all {len(common)} tested rows" if not diff else
                  "q < 0.05 calls differ from the primary set in: " + "; ".join(f"{c} ({lv})" for c, lv in diff))
    bp = n["basis_counts_primary_set"]
    bc = n["basis_counts_connectivity_only"]
    conn_rel = t[(~t.primary_set) & (t.affinity_basis == "relative")].ccd.tolist()
    ent = pd.read_csv(RES / "entry_level.csv")
    top_ent = ent[ent.in_tests_primary.astype(bool)].ccd.value_counts()
    rep = ", ".join(f"{k} {v}" for k, v in top_ent.items() if v >= 3)
    absr = t[t.primary_set & (t.affinity_basis == "absolute")]
    ab_g = absr.sv10_group.value_counts().to_dict()
    ab_c = absr.action_class.value_counts().to_dict()

    line1 = (
        "Not done / substituted: (1) Correction after the first run, disclosed: the first run paired E2 per ligand x "
        "assay x standard type (my declaration), so assay " + ", ".join(n["two_type_assays"]) + " (reports IC50 and "
        f"Ki for both compounds; type-wise dpX differ by <= {n['two_type_max_abs_dpx_difference']:.2f}) counted twice "
        f"for {n['pairs_two_types']} ligands; step 3 asks for one dpX per binding assay, so the reported run uses one "
        "value per ligand x assay (median over its types). This moved the entry-level cl7 vs cl4 test from q "
        f"{f_e.q_bh:.3f} to {cl_e.q_bh:.3f} (ligand level q {f_l.q_bh:.2f} -> {cl_l.q_bh:.2f}; entry-level dpX vs "
        f"ligand volume q {qtxt(f_ve.q_bh)} -> {qtxt(vo_e.q_bh)}); first-run values kept in "
        "affinity_per_ligand.csv dpx_median_type_units and tests_first_run_type_units.csv, not used in the verdict. "
        "(2) Class comparison not tested: coverage gate [DEFAULT] (entries with relative dpX: "
        + ", ".join(f"{c} {v}" for c, v in g["entries_in_sv15_tested_set"].items()) + "; agonist < 5). "
        "(3) Rule 1.4 chain level not run: the step fixes ligand (primary) and entry (secondary) levels and dpX is one "
        "value per ligand. (4) Declared before the outcomes where the run file is open: ligand = CCD (SV10/SV14 "
        "unit) carrying its ChEMBL parent's CHEMBL206 activities; a pair needs the same assay id and the same standard "
        "type; fallback = median over activities, Ki/Kd if any, else IC50; hollow plot position = fallback minus the "
        "E2 median of the same type family (drawing only); spread check on primary-set ligands with >= 2 paired assays, "
        "SD ddof 1, descriptive (no test); entry = SV10 entry ligand and entry group, class = SV12 row (entry, entry "
        "ligand); class test also needs >= 5 entries in the SV15 tested set; Spearman inputs = medians over SV14 "
        "primary-analysis chains whose ligand_ids_used equals the CCD; BH over the primary family, the sensitivity "
        "run (connectivity_only included, run file SV12 step 2) is its own BH family. (5) No network call (every "
        "field was in results/sv/SV12/chembl_activities.csv); nothing installed. (6) Supporting files not named in the "
        "run file: paired_assays.csv, entry_level.csv, tests_first_run_type_units.csv, numbers.json.")

    L = []
    L.append(line1)
    L.append("")
    L.append("# SV15 - Affinities across assay systems, relative to E2 in the same assay")
    L.append("")
    L.append("**Verdict**")
    L.append("")
    verdict = [
        f"- Same-assay E2 reference: {n['paired_assays']} binding assays pair a ligand with E2 ({n['pairs']} "
        f"ligand x assay pairs, {n['paired_parents']} ChEMBL parents); {n['tested_ligands_primary']} primary-set "
        f"ligands (E2 excluded) get a relative dpX, {bp.get('absolute', 0)} only an absolute value; median dpX "
        f"{n['dpx_tested_median']:.2f} (IQR {n['dpx_tested_q1']:.2f} to {fx(n['dpx_tested_q3'])}), "
        f"{n['tested_ligands_dpx_gt0']} of {n['tested_ligands_primary']} bind tighter than E2; "
        f"{n['tested_single_assay_ligands']} rest on one assay.",
        f"- Spread check ({n['spread_ligands']} ligands with >= 2 paired assays): SD across assays median "
        f"{n['spread_median_sd_abs']:.2f} (absolute pX) vs {n['spread_median_sd_dpx']:.2f} (dpX), pooled "
        f"{n['spread_pooled_sd_abs']:.2f} vs {n['spread_pooled_sd_dpx']:.2f}; smaller with the E2 reference in "
        f"{n['spread_n_smaller']} of {n['spread_ligands']}.",
        f"- cl7 vs cl4 (SV10 groups): median dpX {cl_l.median_1:.2f} vs {cl_l.median_2:.2f} ({int(cl_l.n_1)} vs "
        f"{int(cl_l.n_2)} ligands, q {cl_l.q_bh:.2f}); entries {cl_e.median_1:.2f} vs {cl_e.median_2:.2f} "
        f"({int(cl_e.n_1)} vs {int(cl_e.n_2)}, q {cl_e.q_bh:.3f}): not significant at q < 0.05 at either level. "
        f"Classes: not tested "
        f"(agonist {g['entries_in_sv15_tested_set'].get('agonist', 0)} entries with relative dpX, < 5).",
        f"- dpX vs SV14 filling %: rho {fi_l.effect:.2f} ligands (q {fi_l.q_bh:.2f}), {fi_e.effect:.2f} entries "
        f"(q {fi_e.q_bh:.2f}); vs ligand volume: rho {vo_l.effect:.2f} ligands (q {vo_l.q_bh:.2f}, n "
        f"{int(vo_l.n_1)}), {vo_e.effect:.2f} entries (q {qtxt(vo_e.q_bh)}, n {int(vo_e.n_1)}; entries repeat "
        f"ligands: {rep}).",
        "- These are in-vitro binding affinities to the ERa target (CHEMBL206) relative to E2; they say nothing about "
        "cellular activity, receptor degradation or H12 order.",
    ]
    L += verdict
    L.append("")
    L.append("## Inputs (rule 1.10)")
    L.append("")
    L.append("- SV12 card: results/sv/SV12/chembl_activities.csv (1170 CHEMBL206 activities, from the cached ChEMBL "
             "responses), ligand_ids.csv (CCD -> parent, match type; E2 = EST -> CHEMBL135), action_class.csv, "
             "class_coverage.csv (gate classes: " + ", ".join(g["sv12_gate_classes"]) + ").")
    L.append("- SV10 card: results/sv/SV10/ligand_druglikeness.csv (ligand groups, compared = cl7/cl4) and "
             "entry_druglikeness.csv (entry group, entry ligand).")
    L.append("- SV14 card: results/sv/SV14/pocket_per_chain.csv (filling_pct, ligand_volume_A3, in_primary_analysis).")
    L.append("- results/potency_by_cluster.csv is not read or changed; SV15 supersedes it inside the SV section only.")
    L.append("")
    L.append("## Step 1-2 - cleaning and assay systems")
    L.append("")
    L.append(f"- Cleaning cascade (activities kept): raw {cas['raw']['kept']}; target CHEMBL206 "
             f"{cas['target CHEMBL206']['kept']}; relation '=' {cas[REL_KEY]['kept']}; "
             f"pChEMBL present {cas['pChEMBL present']['kept']}; no validity comment "
             f"{cas['no data_validity_comment']['kept']}; not a potential duplicate "
             f"{cas['not a potential duplicate']['kept']} (removed {cas['not a potential duplicate']['removed']}).")
    cs = n["clean_by_system"]
    L.append(f"- Systems: binding (B with Ki/Kd/IC50) {cs.get('binding', 0)} activities "
             f"({', '.join(f'{k} {v}' for k, v in n['binding_act_by_type'].items())}; {n['binding_parents']} parents, "
             f"{n['binding_assays']} assays); functional (F) {cs.get('functional', 0)}, listed only "
             f"({n['functional_ligands']} ligands, column functional_types); dropped {cs.get('dropped', 0)} "
             f"(B with other types: {', '.join(f'{k} {v}' for k, v in n['clean_B_other_types'].items())}; "
             f"A {n['clean_by_assay_type'].get('A', 0)}).")
    L.append(f"- E2 (CHEMBL135): {n['e2_binding_act']} cleaned binding activities in {n['e2_binding_assays']} assays "
             f"({', '.join(f'{k} {v}' for k, v in n['e2_binding_assay_values_by_type'].items())} assay values); "
             f"E2 medians for the hollow plot position: Ki/Kd {n['e2_abs_reference_medians']['Ki/Kd']:.2f}, IC50 "
             f"{n['e2_abs_reference_medians']['IC50']:.2f}.")
    L.append("")
    L.append("## Step 3-4 - relative affinity to E2 and fallback (affinity_per_ligand.csv, paired_assays.csv)")
    L.append("")
    L.append(f"- Pairs: {n['type_level_pairs']} same-assay same-type pairs "
             f"({', '.join(f'{k} {v}' for k, v in n['type_level_pairs_by_type'].items())}) -> {n['pairs']} ligand x "
             f"assay dpX values in {n['paired_assays']} assays; {n['unpaired_type_mismatch_values']} values lost to a "
             f"type mismatch inside an assay; {n['pairs_variant_assays']} pairs in mutant-variant assays.")
    L.append(f"- Per ligand ({n['ligands']} CCDs): " + ", ".join(f"{k} {v}" for k, v in n["basis_counts"].items())
             + f". Primary set: " + ", ".join(f"{k} {v}" for k, v in bp.items()) + "; connectivity_only: "
             + ", ".join(f"{k} {v}" for k, v in bc.items()) + f" (relative: {', '.join(conn_rel)}).")
    ab_t = absr.px_abs_type.value_counts().to_dict()
    L.append(f"- Absolute (hollow, never tested), primary set: {sum(ab_g.values())} ligands (fallback from Ki/Kd "
             f"{ab_t.get('Ki/Kd', 0)}, from IC50 {ab_t.get('IC50', 0)}); "
             f"groups {', '.join(f'{k} {v}' for k, v in ab_g.items())}; classes "
             f"{', '.join(f'{k} {v}' for k, v in ab_c.items())}.")
    L.append("")
    L.append("| CCD | SV10 group | class | dpX median | RBA % | paired assays | filling % (SV14) | entries |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in rel.itertuples():
        L.append(f"| {r.ccd} | {r.sv10_group} | {r.action_class} | {r.dpx_median:+.2f} | {r.rba_median:.3g} | "
                 f"{int(r.n_paired_assays)} | {r.filling_pct_median:.1f} | {int(r.n_entries_sv10_entry_ligand)} |")
    L.append("")
    L.append("## Step 5 - spread check (assay_spread.csv)")
    L.append("")
    L.append("| CCD | paired assays | SD absolute pX | SD dpX | SD E2 pX | ratio |")
    L.append("|---|---|---|---|---|---|")
    for r in sp.itertuples():
        L.append(f"| {r.ccd} | {int(r.n_paired_assays)} | {r.sd_px_abs:.2f} | {r.sd_dpx:.2f} | {r.sd_px_e2:.2f} | "
                 f"{r.sd_ratio_dpx_over_abs:.2f} |")
    L.append("")
    L.append(f"- Median ratio SD(dpX)/SD(pX) {n['spread_median_ratio']:.2f}; median reduction "
             f"{n['spread_median_reduction']:.2f} log units."
             + "".join(f" {r.ccd} ({int(r.n_paired_assays)} assays): SD dpX {r.sd_dpx:.2f} > SD pX {r.sd_px_abs:.2f}, "
                       f"because its E2 values vary more (SD {r.sd_px_e2:.2f})."
                       for r in sp[~sp.dpx_spread_smaller.astype(bool)].itertuples()
                       if r.sd_px_e2 > r.sd_px_abs))
    L.append("")
    L.append("## Step 6 - comparisons (tests.csv; E2 excluded; BH per family)")
    L.append("")
    L.append("| set | comparison | level | n | medians / rho | p | q |")
    L.append("|---|---|---|---|---|---|---|")
    for r in tests.itertuples():
        st = "primary" if r.set == "primary" else "sensitivity"
        if not r.tested:
            L.append(f"| {st} | {r.comparison} | {r.level} | - | not tested: {r.note} | - | - |")
        elif r.test == "Spearman":
            L.append(f"| {st} | {r.comparison} | {r.level} | {int(r.n_1)} | rho {r.effect:+.2f} | {qtxt(r.p)} | "
                     f"{qtxt(r.q_bh)} |")
        else:
            L.append(f"| {st} | {r.comparison} | {r.level} | {int(r.n_1)} vs {int(r.n_2)} | {r.median_1:+.2f} vs "
                     f"{r.median_2:+.2f} (r_rb {r.effect:+.2f}) | {qtxt(r.p)} | {qtxt(r.q_bh)} |")
    L.append("")
    L.append(f"- Class gate: SV12 gate classes {', '.join(g['sv12_gate_classes'])}; entries with relative dpX in the "
             f"SV15 tested set {g['entries_in_sv15_tested_set']}, ligands {g['ligands_in_sv15_tested_set']}; passing "
             f"{g['classes_passing']}. Classes are drawn, not tested.")
    L.append(f"- Sensitivity (connectivity_only included, +{n['tested_ligands_sensitivity'] - n['tested_ligands_primary']} "
             f"ligand): cl7 vs cl4 q {s_cl_l.q_bh:.2f} ligands / {s_cl_e.q_bh:.3f} entries; filling rho "
             f"{s_fi_l.effect:+.2f} / {s_fi_e.effect:+.2f}; ligand volume rho {s_vo_l.effect:+.2f} (q "
             f"{s_vo_l.q_bh:.2f}) / {s_vo_e.effect:+.2f} (q {qtxt(s_vo_e.q_bh)}); {same_calls}.")
    L.append("")
    L.append("## Figures")
    L.append("")
    L.append("- figures/sv/SV15/SV15_dpx_by_group.png/.svg: dpX per ligand by SV10 group and by SV12 class (strip + "
             "box of relative ligands; hollow = absolute, plotted at fallback minus the E2 median of the same type).")
    L.append("- figures/sv/SV15/SV15_dpx_vs_filling.png/.svg: dpX vs SV14 filling % (relative filled and labelled, "
             "absolute hollow).")
    L.append("- Panels: figures/panels/sv/SV15/ (3 PNG).")
    L.append("")
    L.append("## Files")
    L.append("")
    files = [RES / x for x in ("affinity_per_ligand.csv", "assay_spread.csv", "tests.csv", "verdict.md", "SUMMARY.md",
                                "env.txt", "paired_assays.csv", "entry_level.csv", "tests_first_run_type_units.csv",
                                "numbers.json")]
    figs = sorted((ROOT / "figures/sv/SV15").glob("*")) + sorted((ROOT / "figures/panels/sv/SV15").glob("*"))
    scripts = [ROOT / "scripts/sv" / x for x in ("SV15_affinity.py", "SV15_figures.py", "SV15_summary.py")]

    # verdict.md first (<= 5 lines), env.txt, then SUMMARY with sizes
    vtxt = "\n".join(verdict) + "\n"
    (RES / "verdict.md").write_text(vtxt)
    env = [f"# SV15 environment, {time.strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for p in ("numpy", "scipy", "pandas", "statsmodels", "matplotlib"):
        env.append(f"{p} {version(p)}")
    env += ["ChEMBL data: results/sv/SV12/chembl_activities.csv (SV12, ChEMBL_37, cache data/sv/chembl_cache/)",
            "network: none (no field missing from the SV12 activity table)",
            "nothing installed, upgraded or downgraded by SV15"]
    (RES / "env.txt").write_text("\n".join(env) + "\n")
    head = list(L)
    (RES / "SUMMARY.md").write_text("\n".join(head) + "\n")  # first pass so SUMMARY.md has a size
    for _ in range(3):  # converge the listed size of SUMMARY.md itself
        L = head + [f"- {p.relative_to(ROOT)} ({size(p)})" for p in files + figs + scripts]
        L.append("- docs/EXPLORATION.md: appended 'Suggested, not done (SV15, 2026-11-08)' (rule 1.1)")
        (RES / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
