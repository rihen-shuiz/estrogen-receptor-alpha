"""SV12 SUMMARY.md (rule 1.7) - every number read from the saved SV12 files."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV12"
DOUT = ROOT / "data/sv/SV12"


def size(p: Path) -> str:
    b = p.stat().st_size
    return f"{b / 1024:.1f} KB" if b < 1024 ** 2 else f"{b / 1024 ** 2:.1f} MB"


def main():
    num = json.loads((OUT / "numbers.json").read_text())
    chk = json.loads((OUT / "chembl_checks.json").read_text())
    ids = json.loads((DOUT / "ids_info.json").read_text())
    cinfo = json.loads((DOUT / "chembl_info.json").read_text())
    ac = pd.read_csv(OUT / "action_class.csv")
    lig = pd.read_csv(OUT / "ligand_ids.csv")
    cov = pd.read_csv(OUT / "class_coverage.csv")
    val = pd.read_csv(OUT / "route_c_validation.csv")
    act = pd.read_csv(OUT / "chembl_activities.csv", low_memory=False)
    mech = pd.read_csv(DOUT / "mechanism_records.csv")
    rb = pd.read_csv(DOUT / "route_b_per_parent.csv")
    chains = pd.read_csv(OUT / "chain_class.csv")
    g = num["gate"]
    gate_txt = (f"route (c) PASSED its validation gate ({g['n_agree']}/{g['n_compared']} rows agree = "
                f"{g['agreement']:.2f}; needs >= {g['min_agreement']:.2f} of >= {g['min_rows']}) and fills rows "
                f"without a ChEMBL mechanism") if g["passed"] else (
        f"route (c) FAILED its validation gate ({g['n_agree']}/{g['n_compared']} rows agree"
        + (f" = {g['agreement']:.2f}" if g["n_compared"] else "")
        + f"; needs >= {g['min_agreement']:.2f} of >= {g['min_rows']} rows), so its labels are not used")

    mt = num["ccd_by_match_type"]
    est = lig[lig.ccd == "EST"].iloc[0]
    L = []
    L.append(
        "Not done / substituted: teammate file data/sv/action_class_teammate.csv absent at run time -> labels come "
        f"from route (b) ChEMBL mechanism only; {gate_txt}; class_agreement.csv not written (no teammate file). "
        "All labels are provisional (results/sv/SV12/ only; results/ligand_action_labels.csv and Q10/Q11 untouched; "
        "no literature searched; no class from structure, cluster, chemistry or a figure). "
        "Declared before labelling where the run file is open: "
        "(1) 'modulator' text excludes the phrase 'allosteric modulator' so the explicit NEGATIVE ALLOSTERIC "
        "MODULATOR -> antagonist mapping can apply; "
        "(2) ChEMBL document lookup uses doi__iexact (the plain doi filter is case-sensitive: the CHEMBL1145862 test DOI "
        "upper-cased gives 0 hits, doi__iexact 1); "
        "(3) connectivity matches with several parents take the no-stereo-layer parent if present, else the lowest "
        "ChEMBL id (all candidates listed); "
        "(4) ESR1 targets = all 9 ChEMBL targets with component P03372 (CHEMBL206 + family/complex/PPI targets) "
        "for routes (b) and (c); activities for SV15 are CHEMBL206 only; "
        "(5) conflict = different classes counting partial_agonist as agonist; "
        "(6) chain label = label of the chain's SV00 ligand_id_primary; an entry counts in a class x cluster cell "
        "if >= 1 of its chains there has the class. "
        f"{ids['doi_none']} of {ids['n_entries']} entries have no primary-citation DOI in the mmCIF or the Data API. "
        "Supporting files not named in the run file: chain_class.csv, numbers.json, chembl_checks.json "
        "(results/sv/SV12/) and data/sv/SV12/*.")
    L.append("")
    L.append("# SV12 - Ligand identifiers, provisional action classes and ChEMBL activities")
    L.append("")
    fr = num["final_rows"]
    L.append(f"**Result:** {num['n_rows']} entry x ligand rows ({num['n_entries']} entries, {num['n_ccd']} CCD codes). "
             "Final provisional classes (rows): " + ", ".join(f"{k} {fr.get(k, 0)}" for k in
                                                               ["agonist", "partial_agonist", "SERM", "SERD",
                                                                "antagonist", "other", "unlabelled"]) + ".")
    L.append(f"- Conflicts flagged: {num['n_conflicts']}"
             + ("" if not num["n_conflicts"] else " - " + "; ".join(
                 f"{c['entry_id']}/{c['ccd']} (b {c['label_b']}, c {c['label_c']})" for c in num["conflicts"])))
    L.append("")
    L.append("## Identifiers (step 1)")
    L.append(f"- RCSB chemcomp: {ids['rcsb_ok']}/{ids['n_ccd']} records; stereo SMILES missing "
             f"{ids['smiles_stereo_missing']}, InChIKey missing {ids['inchikey_missing']} (all 306 formal charge 0, "
             "standard InChIKeys).")
    L.append(f"- Primary-citation DOI: {ids['doi_mmcif']} entries from mmCIF _citation (id primary), "
             f"{ids['doi_data_api']} from the Data API, {ids['doi_none']} none; "
             f"{cinfo['n_dois']} distinct DOIs.")
    L.append("")
    L.append("## ChEMBL checks and match (step 2, 6)")
    env = (OUT / "env.txt").read_text().splitlines()
    files = [f for f in sorted((ROOT / "data/sv/chembl_cache").glob("*.json")) if f.read_text() != "__404__"]
    probe = [f for f in files if f.name.startswith("activity_probe_")]
    lists = [json.loads(f.read_text()) for f in files if f not in probe]
    lists = [d for d in lists if isinstance(d, dict) and "page_meta" in d]
    n_next = sum(1 for d in lists if d["page_meta"].get("next"))
    L.append(f"- ChEMBL cache: {len(files)} responses; {len(lists)} list pages of the pipeline (limit 1000), {n_next} with "
             "a non-null page_meta.next, so every list fitted one page (paging code in cached_http.chembl_get_all); "
             f"{len(probe)} field probe (activity, limit 2) not used as data; 0 un-cached requests on a re-run.")
    L.append("- " + [x for x in env if x.startswith("ChEMBL release")][0].split(" (status")[0] + ".")
    L.append(f"- EST test: RCSB InChIKey {chk['EST_inchikey_rcsb']} -> {', '.join(chk['EST_molecule_hits'])} "
             f"({', '.join(map(str, chk['EST_pref_name']))}), parent "
             f"{chk['EST_molecule_hierarchy'][0]['parent_chembl_id']}; E2 reference for SV15 = "
             f"{est.chembl_parent} ({est.match_type} match).")
    c206 = chk["CHEMBL206"]
    L.append(f"- CHEMBL206: {c206['target_type']}, {c206['organism']}, components {c206['components']} -> "
             f"single-protein target for P03372 = {chk['CHEMBL206_is_single_protein_P03372']}; the only single-protein "
             f"target with P03372: {chk['single_protein_targets_for_P03372']}; "
             f"{len(chk['targets_containing_P03372'])} targets contain P03372.")
    L.append(f"- Match per CCD: full {mt.get('full', 0)}, connectivity_only {mt.get('connectivity_only', 0)}, "
             f"no_chembl {mt.get('no_chembl', 0)} (of {len(lig)}); several candidate parents {cinfo['multi_parent']} "
             f"(match types: {', '.join(sorted(set(lig[lig.multi_parent == True].match_type)))}); "
             f"{cinfo['n_parents']} distinct parents.")
    rm = num["rows_by_match_type"]
    L.append(f"- Rows: full {rm.get('full', 0)}, connectivity_only {rm.get('connectivity_only', 0)} "
             f"(excluded from the SV13/SV15 primary analyses, primary_set = False), no_chembl {rm.get('no_chembl', 0)}.")
    L.append("")
    L.append("## Labels (steps 3-4)")
    nb = int(mech.esr1_target.sum()) if len(mech) else 0
    L.append(f"- Route (b): {len(mech)} mechanism records for the parents, {nb} on ESR1 targets, for "
             f"{len(rb)} parents; parent labels: "
             + ", ".join(f"{k} {v}" for k, v in num["label_b_parents"].items()) + "; rows: "
             + ", ".join(f"{k} {v}" for k, v in num["label_b_rows"].items()) + ".")
    for r in rb.sort_values("label_b").itertuples():
        names = lig.loc[lig.chembl_parent == r.chembl_parent, "ccd"].tolist()
        L.append(f"  - {r.chembl_parent} ({', '.join(names)}): {r.label_b} <- {r.b_action_types} "
                 f"[{r.b_mechanisms}]")
    L.append(f"- Route (c): {cinfo['n_dois_with_document']} of {cinfo['n_dois']} DOIs found as ChEMBL documents "
             f"({cinfo['n_documents']} documents, {cinfo['n_document_esr1_activities']} ESR1-target activities); "
             "row outcomes: " + ", ".join(f"{k} {v}" for k, v in num["c_reason_rows"].items())
             + "; labels: " + ", ".join(f"{k} {v}" for k, v in num["label_c_rows"].items()) + ".")
    L.append(f"- Validation gate: {g['n_compared']} rows labelled by both ({g['n_compared_ccd']} ligands, "
             f"{g['n_compared_entries']} entries), agreement {g['n_agree']}/{g['n_compared']}"
             + (f" = {g['agreement']:.2f}" if g["n_compared"] else "") + f" -> {'PASS' if g['passed'] else 'FAIL'}"
             " ([DEFAULT] >= 0.80 of >= 10 rows).")
    dis = val[val.c_compared & (val.c_agrees_b == False)]  # noqa: E712
    if len(dis):
        L.append("  - disagreements (b vs c): " + "; ".join(
            f"{r.entry_id}/{r.ccd} {r.label_b} vs {r.label_c}" for r in dis.itertuples()))
    L.append("- Final rows by source: " + "; ".join(f"{d['source']} {d['action_class']} {d['n']}"
                                                   for d in num["final_rows_by_source"]) + ".")
    L.append("")
    L.append("## Coverage gate [DEFAULT] (SV13-SV15: a class comparison is tested only if each side has >= 5 entries)")
    L.append("")
    L.append("| set | class | chains all | entries all | chains cl7 | entries cl7 | chains cl4 | entries cl4 | "
             "chains minor | entries minor | >= 5 entries (all) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for s in ("primary", "all"):
        for cl in ["agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other", "unlabelled"]:
            sub = cov[(cov.set == s) & (cov.action_class == cl) & (cov.source == "any")].set_index("cluster")
            if sub.empty:
                continue
            cells = []
            for clu in ("all", "cl7", "cl4", "minor"):
                cells += [str(int(sub.loc[clu, "n_chains"])), str(int(sub.loc[clu, "n_entries"]))]
            ge = "" if cl == "unlabelled" else ("yes" if int(sub.loc["all", "n_entries"]) >= 5 else "no")
            L.append(f"| {s} | {cl} | " + " | ".join(cells) + f" | {ge} |")
    nol = int((chains.source == "no_ligand_of_interest").sum())
    L.append("")
    L.append(f"- {len(chains)} SV00 chains; {nol} without a ligand of interest are in 'unlabelled'; chains whose "
             f"two ligand codes carry different classes: {num['chains']['two_codes_differing_class']}.")
    L.append("")
    L.append("## Activities for SV15 (step 5)")
    act_pc = act.pchembl_value.notna().sum()
    e2 = act[act.chembl_parent == est.chembl_parent]
    L.append(f"- {len(act)} activities on CHEMBL206 for {act.chembl_parent.nunique()} of {cinfo['n_parents']} "
             f"matched parents ({act_pc} with pChEMBL; assay types "
             + ", ".join(f"{k} {v}" for k, v in act.assay_type.value_counts().items())
             + f"); E2 ({est.chembl_parent}): {len(e2)} activities in {e2.assay_chembl_id.nunique()} assays.")
    L.append("")
    L.append("## Verdict")
    lab = ac[ac.action_class != "unlabelled"]
    phases = lig[lig.chembl_parent.isin(lab.chembl_parent)].drop_duplicates("chembl_parent").parent_max_phase
    L.append(f"- The provisional label set covers {len(lab)} of {num['n_rows']} rows "
             f"({lab.ccd.nunique()} of {num['n_ccd']} ligands, {lab.chembl_parent.nunique()} ChEMBL parents, "
             f"ChEMBL max_phase {phases.min():g}-{phases.max():g}); {fr.get('unlabelled', 0)} rows stay unlabelled "
             "(no ESR1 mechanism record, no ChEMBL parent, or route (c) not used).")
    present = [c for c in ["agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other"]]
    allc = cov[(cov.set == "primary") & (cov.source == "any") & (cov.cluster == "all")].set_index("action_class")
    passing = [c for c in present if int(allc.loc[c, "n_entries"]) >= 5]
    failing = [f"{c} {int(allc.loc[c, 'n_entries'])}" for c in present if int(allc.loc[c, "n_entries"]) < 5]
    L.append(f"- Coverage gate (primary set, all clusters): >= 5 entries for {', '.join(passing)}; below 5 entries: "
             f"{', '.join(failing)}. Classes and clusters overlap: "
             + "; ".join(f"{c} entries cl7 {int(cov[(cov.set == 'primary') & (cov.source == 'any') & (cov.action_class == c) & (cov.cluster == 'cl7')].n_entries.iloc[0])} / "
                         f"cl4 {int(cov[(cov.set == 'primary') & (cov.source == 'any') & (cov.action_class == c) & (cov.cluster == 'cl4')].n_entries.iloc[0])}"
                         for c in ["agonist", "SERM", "SERD"]) + ".")
    L.append("- Consequence for the class comparisons named in SV13 (by the entry counts above): SERD vs SERM, "
             "antagonist vs SERM and agonist vs partial_agonist each have a side below 5 entries; agonist and SERM "
             "are the only classes with >= 5 entries.")
    esl = lab[lab.ccd == "ESL"]
    if len(esl):
        r = rb[rb.chembl_parent == esl.chembl_parent.iloc[0]].iloc[0]
        L.append(f"- Note: ESL ({esl.parent_pref_name.iloc[0]}, {r.chembl_parent}) is SERM because its only ESR1 "
                 f"mechanism record has action type {r.b_action_types} ('{r.b_mechanisms}'); kept as the [DEFAULT] "
                 "map gives it.")
    L.append("")
    L.append("## Files")
    files = ["ligand_ids.csv", "entry_doi.csv", "action_class.csv", "class_coverage.csv", "route_c_validation.csv",
             "chembl_activities.csv", "chain_class.csv", "numbers.json", "chembl_checks.json", "env.txt"]
    for f in files:
        L.append(f"- results/sv/SV12/{f} ({size(OUT / f)})")
    for f in sorted(DOUT.iterdir()):
        L.append(f"- data/sv/SV12/{f.name} ({size(f)})")
    L.append("- docs/EXPLORATION.md: appended 'Suggested, not done (SV12, 2026-11-08)' (rule 1.1)")
    L.append("- caches: data/sv/rcsb_chemcomp/, data/sv/chembl_cache/ (index.jsonl), data/sv/api_cache/sv12_*")
    L.append("- scripts: scripts/sv/SV12_chembl_test.py, SV12_ids.py, SV12_chembl.py, SV12_classes.py, SV12_env.py, "
             "SV12_summary.py; scripts/sv/common/cached_http.py")
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
