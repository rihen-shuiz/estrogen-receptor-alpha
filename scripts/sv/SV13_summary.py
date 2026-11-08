"""SV13 (RUN_SV1028): verdict.md, SUMMARY.md (rule 1.7), env.txt and numbers.json, from the saved SV13 files only.

Inputs: results/sv/SV13/class_by_cluster.csv, chain_data.csv, gate.csv, tests.csv, numbers_tests.json,
representatives.csv, view.txt; data/sv/SV13/panel_meta.json, model_code_check.json; figure files.
Run: .venv/bin/python -B scripts/sv/SV13_summary.py
"""
from __future__ import annotations

import importlib
import json
import platform
import shutil
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results/sv/SV13"
DAT = ROOT / "data/sv/SV13"


def kb(p: Path) -> str:
    return f"{p.stat().st_size / 1024:.1f} KB"


def main():
    cbc = pd.read_csv(RES / "class_by_cluster.csv")
    d = pd.read_csv(RES / "chain_data.csv")
    gate = pd.read_csv(RES / "gate.csv")
    tests = pd.read_csv(RES / "tests.csv")
    nt = json.loads((RES / "numbers_tests.json").read_text())
    rep = pd.read_csv(RES / "representatives.csv")
    meta = json.loads((DAT / "panel_meta.json").read_text())
    cc = json.loads((DAT / "model_code_check.json").read_text())

    def cell(cls, grp, st="primary"):
        r = cbc[(cbc.set == st) & (cbc.action_class == cls) & (cbc.s08_cluster == grp)]
        if not len(r):
            r = cbc[(cbc.set == st) & (cbc.action_class == cls) & (cbc.cluster_group == grp)
                    & cbc.s08_cluster.str.startswith("minor")]
        r = r.iloc[0]
        return int(r.n_chains), int(r.n_entries), int(r.n_complete_h12)

    g = gate[gate.set == "primary"].set_index(["comparison", "stratum"])
    gs = gate[gate.set == "all"].set_index(["comparison", "stratum"])
    rows = tests[tests.set == "primary"].drop_duplicates(["comparison", "stratum", "outcome"])
    n_rows, n_pass = len(rows), int(rows.gate_pass.astype(str).eq("True").sum())
    cv = nt["comparison_verdicts"]
    r = rep.set_index("action_class")
    agn, serm, serd = cell("agonist", "all"), cell("SERM", "all"), cell("SERD", "all")
    c2b = g.loc[("C2b_agonist_rho_dpX", "cl7")]
    psi = nt["psi537"]
    lchk = meta["ligand_check"]
    n_conn = int((~d.primary_set).sum())
    conn_cls = d.loc[~d.primary_set, "action_class"].value_counts().to_dict()

    # ---------------------------------------------------------------- verdict.md (<= 6 lines)
    verdict = [
        f"- C1 SERD vs SERM (H12 complete, fraction observed, mobility z; all chains and within cl4): "
        f"**{cv['primary|C1_SERD_vs_SERM']}** - SERD {int(g.loc[('C1_SERD_vs_SERM', 'all')].entries_a)} entries "
        f"(all in cl4) vs SERM {int(g.loc[('C1_SERD_vs_SERM', 'all')].entries_b)} "
        f"({int(g.loc[('C1_SERD_vs_SERM', 'cl4')].entries_b)} in cl4); coverage gate needs >= 5 entries per side.",
        f"- C1 antagonist vs SERM: **{cv['primary|C1_antagonist_vs_SERM']}** (0 antagonist entries in SV12). "
        f"C2a partial agonist vs agonist within cl7: **{cv['primary|C2a_partial_agonist_vs_agonist']}** "
        f"(0 partial-agonist entries; {int(g.loc[('C2a_partial_agonist_vs_agonist', 'cl7')].entries_b)} agonist entries in cl7).",
        f"- C2b agonists vs dpX within cl7 (the VDR superagonist point): **{cv['primary|C2b_agonist_rho_dpX']}** "
        f"({int(c2b.entries_a)} agonist entries with an SV15-tested dpX: {c2b.dpx_ligands.replace(':', ' ').replace(';', ', ')} entries; "
        "E2 is the reference and 4OH has an absolute value only).",
        f"- So no class-H12 link can be stated: {n_pass} of {n_rows} comparison x stratum x outcome rows pass the gate "
        f"(sensitivity set with connectivity_only rows: same counts, {n_conn} added chains carry no class).",
        f"- Drawn, not tested (panels, maps): chains with complete H12 538-547 agonist {agn[2]} of {agn[0]}, SERM "
        f"{serm[2]} of {serm[0]}, SERD {serd[2]} of {serd[0]}; classes follow the clusters (agonist "
        f"{cell('agonist', 'cl7')[1]} of {agn[1]} entries in cl7, SERM {cell('SERM', 'cl4')[1]} of {serm[1]} and SERD "
        f"{cell('SERD', 'cl4')[1]} of {serd[1]} in cl4).",
        "- Crystal structures show H12 order and position (H12 disorder, H12 mobility), not receptor degradation or "
        "cellular activity; the classes are provisional ChEMBL-mechanism labels (SV12).",
    ]
    (RES / "verdict.md").write_text("\n".join(verdict) + "\n")

    # ---------------------------------------------------------------- numbers.json
    nums = dict(class_by_cluster_primary={c: dict(zip(("chains", "entries", "complete_h12"), cell(c, "all")))
                                          for c in ("agonist", "SERM", "SERD", "unlabelled")},
                gate_rows_primary=n_rows, gate_rows_pass=n_pass, comparison_verdicts=cv,
                representatives=rep[rep.n_chains > 0][["action_class", "representative", "rep_ligand", "rep_cluster",
                                                         "rep_median_rmsd_A", "rep_resolution_A",
                                                         "rep_h12_complete"]].to_dict(orient="records"),
                view=dict(chain=meta["view"]["view_chain"], n_atoms=meta["view"]["n_atoms"],
                          limits=meta["limits"]),
                ligand_check=lchk, b_z_colour_range=meta["colour_range_b_z"], psi537=psi,
                connectivity_only_chains=n_conn, connectivity_only_classes=conn_cls,
                code_check_rows=len(cc["rows"]), code_check_fitted=sum(x["status"].startswith("fitted")
                                                                       for x in cc["rows"]))
    (RES / "numbers.json").write_text(json.dumps(nums, indent=1, default=str))

    # ---------------------------------------------------------------- env.txt
    tz = timezone(timedelta(hours=5))
    env = [f"# SV13 environment, {datetime.now(tz).strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    for m in ("numpy", "pandas", "scipy", "statsmodels", "matplotlib", "Bio"):
        mod = importlib.import_module(m)
        env.append(f"{'biopython' if m == 'Bio' else m} {mod.__version__}")
    try:
        importlib.import_module("pymol")
        env.append("pymol: importable")
    except ImportError:
        env.append(f"pymol: not importable, executable {shutil.which('pymol')} -> Part C uses the step 12 fallback "
                   "(matplotlib C-alpha traces projected onto the saved view)")
    env += ["mmCIF reading: esr_analysis.read_structure.read_cif (Bio.PDB.MMCIF2Dict), restored files data/raw/",
            "network: none", "nothing installed, upgraded or downgraded by SV13"]
    (RES / "env.txt").write_text("\n".join(env) + "\n")

    # ---------------------------------------------------------------- SUMMARY.md
    line1 = (
        "Not done / substituted: (1) PyMOL not available (PRE addition; env.txt) -> Part C uses the step 12 fallback: "
        "C-alpha traces projected with matplotlib onto the saved view; the representative's cartoon is a thick "
        "C-alpha trace, its ligand of interest sticks with bonds from interatomic distances. (2) Part B: every "
        f"comparison fails the coverage gate ({n_pass} of {n_rows} rows testable), so no model or rank correlation "
        "was fitted and the forest plot shows the gate reason in place of estimates. (3) Representatives follow the "
        f"[DEFAULT] as written: the SERM representative {r.loc['SERM', 'representative']} has no modelled H12 "
        f"(only 528 of 528-550 observed) and the SERD representative {r.loc['SERD', 'representative']} 7 of 10 H12 "
        "residues, so their B z colouring is absent / partial (the existing RMSD matrix compares chains on their "
        "common positions). (4) Declared before any class comparison, where the run file is open: sensitivity set "
        f"= all rows incl. connectivity_only ({n_conn} chains, all unlabelled, so counts are identical); C2b uses the "
        "ligands SV15 itself tests (relative same-assay dpX; E2 reference and absolute-only ligands excluded, SV15 "
        "steps 4 and 6) with a gate of >= 5 agonist entries with dpX; psi 537 outcome = wrapped difference from the "
        f"circular mean of {psi['n_cl7_chains_with_psi537']} cl7 chains ({psi['cl7_circular_mean_psi537_deg']:.1f} deg); "
        "representative ties by resolution, then chain_key; view = principal axes of the 5TM1|A|1 selection (PyMOL "
        "'orient' convention), signs with its H12 centroid at +x/+y. (5) Shared code: optional edge_map parameter "
        "added to scripts/sv/common/h12_map.py (default output unchanged). (6) The model code was exercised once on "
        "seeded random pseudo-labels (data/sv/SV13/model_code_check.json: statuses only, no estimate used) so that "
        "the re-run with the teammate labels has a checked path. (7) Supporting files not named in the run file: "
        "chain_data.csv, gate.csv, numbers_tests.json, numbers.json (results/sv/SV13/), data/sv/SV13/*.")

    def ce(cls, grp):
        c = cell(cls, grp)
        return f"{c[0]} / {c[1]}"

    tab = ["| class | chains / entries all | cl7 | cl4 | minor | complete H12 chains |", "|---|---|---|---|---|---|"]
    for c in ("agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other", "unlabelled"):
        tab.append(f"| {c} | {ce(c, 'all')} | {ce(c, 'cl7')} | {ce(c, 'cl4')} | {ce(c, 'minor')} | "
                   f"{cell(c, 'all')[2]} |")
    gtab = ["| comparison | stratum | side a (chains / entries) | side b (chains / entries) | gate | verdict |",
            "|---|---|---|---|---|---|"]
    for (comp, st), x in g.iterrows():
        sb = "-" if pd.isna(x.chains_b) else f"{x.side_b} {int(x.chains_b)} / {int(x.entries_b)}"
        gtab.append(f"| {comp} | {st} | {x.side_a} {int(x.chains_a)} / {int(x.entries_a)} | {sb} | "
                    f"{'pass' if x.gate_pass else 'fail: ' + str(x.gate_reason)} | {cv['primary|' + comp]} |")
    rtab = ["| class | representative | ligand | cluster | median RMSD to class (A) | resolution (A) | H12 538-547 modelled | runner-up |",
            "|---|---|---|---|---|---|---|---|"]
    tr = pd.read_csv(DAT / "panel_traces.csv")
    for x in rep[rep.n_chains > 0].itertuples():
        nh = int(tr[(tr.chain_key == x.representative) & (tr.role == "class_trace")
                    & tr.position.between(538, 547)].observed.sum())
        rtab.append(f"| {x.action_class} | {x.representative} | {x.rep_ligand} | {x.rep_cluster} | "
                    f"{x.rep_median_rmsd_A:.3f} | {x.rep_resolution_A:.2f} | {nh} of 10 | {x.second_chain} "
                    f"({x.second_median_rmsd_A:.3f}) |")
    figs = sorted((ROOT / "figures/sv/SV13").glob("*")) + sorted((ROOT / "figures/panels/sv/SV13").glob("*"))
    res_files = [RES / f for f in ("class_by_cluster.csv", "tests.csv", "representatives.csv", "view.txt",
                                   "verdict.md", "chain_data.csv", "gate.csv", "numbers_tests.json", "numbers.json",
                                   "env.txt")]
    dat_files = sorted(DAT.glob("*"))
    scripts = [ROOT / "scripts/sv" / f for f in ("SV13_analysis.py", "SV13_panels.py", "SV13_figures.py",
                                                  "SV13_summary.py")]
    view_lines = (RES / "view.txt").read_text().splitlines()
    s = [line1, "",
         "# SV13 - H12 by action class (provisional SV12 classes; H12 window 538-547)", "",
         "**Verdict** (results/sv/SV13/verdict.md)", ""] + verdict + [
        "", "## Inputs (rule 1.10)", "",
        "- SV12 card: results/sv/SV12/chain_class.csv (class per chain, primary_set), class_coverage.csv (cross-check: "
        f"{nt['coverage_crosscheck_vs_SV12']['chains_equal']} of {nt['coverage_crosscheck_vs_SV12']['rows_compared']} "
        "class x cluster cells equal).",
        "- SV00 chain_table.csv; SV01a chain_outcomes.csv (complete, fraction observed); SV01b mobility_per_chain.csv, "
        "mobility_per_residue.csv, lattice_contacts_per_chain.csv, model_data.csv (rebuilt Q/X/lattice/mutation/peptide "
        f"predictors equal model_data.csv: max |diff| cl7 {nt['predictor_rebuild_check']['cl7']['max_abs_diff_rebuilt_vs_model_data']}, "
        f"cl4 {nt['predictor_rebuild_check']['cl4']['max_abs_diff_rebuilt_vs_model_data']}); SV02 phipsi_per_chain.csv; "
        "SV03 geometry_per_chain.csv, core_v3_transforms.csv; SV15 affinity_per_ligand.csv.",
        "- Y03 card: results/ext/embedding_rmsd.csv (mds1/mds2); S04/S05 cards: results/rmsd_entry.npy + "
        "rmsd_entry_index.csv; S12: results/medoids.csv (cl7 medoid 5TM1|A|1); Q18: results/chain_ligand_of_interest.csv; "
        "restored mmCIFs data/raw/.",
        "", "## Part A - classes on the existing maps (class_by_cluster.csv; primary set)", ""] + tab + [
        "", f"- class_by_cluster.csv also lists every S08 cluster id (1-7) and the 'all' set; connectivity_only chains "
        f"({n_conn}) are all unlabelled.",
        "- Maps redrawn with plot_h12_map (SV03 function): fill = S08 cluster, marker shape and edge = SV12 class, "
        "unlabelled small grey-edged circles. MDS map: 788 chains; H12 map: 748 chains with an H12 axis.",
        "", "## Part B - tests fixed in advance (tests.csv, gate.csv)", ""] + gtab + [
        "", f"- tests.csv: {len(tests)} rows (2 sets x {n_rows} comparison x stratum x outcome rows x 2 levels), all "
        "'not testable (coverage gate)'; models and C2b code are in scripts/sv/SV13_analysis.py "
        "(class + Q + X with lattice terms + mutations + peptide; chain level entry-clustered, entry level one chain per "
        "entry; C2b Spearman, chain level bootstrap over entries; BH within each comparison).",
        f"- Sensitivity set (connectivity_only included): identical gate counts; verdicts "
        f"{sorted(set(v for k, v in cv.items() if k.startswith('all|')))}.",
        f"- Code-path check on seeded random pseudo-labels: {nums['code_check_fitted']} of {nums['code_check_rows']} "
        "model / rho fits ran (statuses only, no estimate kept).",
        "", "## Part C - pictures of H12 and the ligand (representatives.csv, view.txt)", ""] + rtab + [
        "", "- Ligand atoms (restored mmCIF, SV14 instance rule) vs data/sv/SV14/atoms_core.npz: " + "; ".join(
            f"{k} {v['n_atoms']} atoms, max {v['max_abs_diff_A']:.1e} A" for k, v in lchk.items()) + ".",
        f"- H12 B z colour range {meta['colour_range_b_z']} (representatives' H12 B z "
        f"{meta['b_z_range_rep_h12'][0]:.2f} to {meta['b_z_range_rep_h12'][1]:.2f}; SV01b per-residue z vs the chain's "
        "core v3 C-alpha B).",
        "- Panels: 6 class panels (partial agonist, antagonist and other say 0 chains) + 1 overlay of the 3 "
        "representatives; captions give chains, entries and chains with complete H12.",
        "- view.txt:"] + [f"  - {v}" for v in view_lines[1:]] + [
        "", "## Figures", "",
        "- figures/sv/SV13/SV13_h12_class_panels.png/.svg (main): per-class H12 panels + overlay, one saved view.",
        "- figures/sv/SV13/SV13_forest.png/.svg (main): C1/C2a class effects and C2b rho; every row shows its gate reason.",
        "- figures/sv/SV13/SV13_mds_map.png/.svg and SV13_h12_map.png/.svg (backup): the two redrawn maps.",
        "", "## Files", ""]
    s += [f"- {p.relative_to(ROOT)} ({kb(p)})" for p in res_files + dat_files + figs + scripts]
    s += ["- results/sv/SV13/SUMMARY.md (this file)",
          "- scripts/sv/common/h12_map.py: optional edge_map parameter (SV13)",
          "- docs/EXPLORATION.md: appended 'Suggested, not done (SV13, 2026-11-08)' (rule 1.1)"]
    (RES / "SUMMARY.md").write_text("\n".join(s) + "\n")
    print("\n".join(s))


if __name__ == "__main__":
    main()
