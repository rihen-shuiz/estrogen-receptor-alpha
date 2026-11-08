"""SV12 steps 3-4 and the coverage gate - provisional action classes.

Unit: one row per entry x ligand of interest (distinct entry_id, ligand_id of results/chain_ligand_of_interest.csv).
Sources in order of precedence (run file SV12 step 4):
  (a) data/sv/action_class_teammate.csv if it exists at run time (then ChEMBL labels only cross-check it);
  (b) ChEMBL mechanism records of the parent on ESR1 targets (data/sv/SV12/mechanism_records.csv);
  (c) assay descriptions of the ligand's activities (pChEMBL >= 6, ESR1 targets) in the entry's own paper (DOI ->
      ChEMBL document), used only if it passes the validation gate against (b).
Never from structure, cluster, chemistry or figures (rule 1.9).

Declared before labelling (points the run file leaves open):
- (b) per record: text 'degrad'/'downregulat' or action type DEGRADER -> SERD; action type MODULATOR or text
  'modulator' -> SERM, except that the phrase 'allosteric modulator' (the name of the NEGATIVE/POSITIVE ALLOSTERIC
  MODULATOR action types) does not count as 'modulator' text, so the explicit NAM -> antagonist mapping can apply;
  AGONIST -> agonist; PARTIAL AGONIST -> partial_agonist; ANTAGONIST / INVERSE AGONIST / NEGATIVE ALLOSTERIC
  MODULATOR -> antagonist; a record that gives none of these -> other. Over all ESR1 records of the parent the
  labels are pooled, SERM is added when there are agonist-type (AGONIST, PARTIAL AGONIST) and antagonist-type
  (ANTAGONIST, INVERSE AGONIST, NEGATIVE ALLOSTERIC MODULATOR) records, and the highest of
  SERD > SERM > antagonist > partial_agonist > agonist > other is taken.
- (c) patterns are applied per distinct assay description; the description sets are pooled per row.
- Conflict = two used sources give different classes, counting partial_agonist as agonist (as in the gate).
- Chain label = the label of the chain's ligand_id_primary (SV00: first ligand-of-interest row of the chain);
  chains without a ligand of interest are counted as 'no_ligand_of_interest' in the coverage table.
- Coverage: an entry counts in a class x cluster cell if at least one of its chains in that cluster has the class;
  'primary' set = connectivity_only rows left out (run file SV12 step 2).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV12"
DOUT = ROOT / "data/sv/SV12"
TEAMMATE = ROOT / "data/sv/action_class_teammate.csv"

CLASSES = ["agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other"]
PREC = ["SERD", "SERM", "antagonist", "partial_agonist", "agonist", "other"]
AGO_T = {"AGONIST", "PARTIAL AGONIST"}
ANT_T = {"ANTAGONIST", "INVERSE AGONIST", "NEGATIVE ALLOSTERIC MODULATOR"}

RE_DEG = re.compile(r"degrad|downregulat|down-regulat")
RE_ANT = re.compile(r"antagonis|inhibition of (e2|estradiol|17.?beta-estradiol)[- ]induced")
RE_AGO = re.compile(r"(?<!ant)agonis|transactivation|transcriptional activation")
PCHEMBL_MIN = 6.0
GATE_MIN_ROWS, GATE_MIN_AGREE = 10, 0.80
COVERAGE_MIN_ENTRIES = 5


def record_labels(action_type, moa) -> set[str]:
    at = (action_type or "").strip().upper()
    txt = (moa or "").lower() if isinstance(moa, str) else ""
    labs = set()
    if "degrad" in txt or "downregulat" in txt or at == "DEGRADER":
        labs.add("SERD")
    if at == "MODULATOR" or "modulator" in txt.replace("allosteric modulator", ""):
        labs.add("SERM")
    if at == "AGONIST":
        labs.add("agonist")
    elif at == "PARTIAL AGONIST":
        labs.add("partial_agonist")
    elif at in ANT_T:
        labs.add("antagonist")
    if not labs:
        labs.add("other")
    return labs


def top(labs: set[str]):
    for c in PREC:
        if c in labs:
            return c
    return "unlabelled"


def route_b(mech: pd.DataFrame) -> pd.DataFrame:
    rows = []
    m = mech[mech.esr1_target.astype(bool)]
    for p, g in m.groupby("chembl_parent"):
        labs = set()
        for r in g.itertuples():
            labs |= record_labels(r.action_type, r.mechanism_of_action)
        types = set(g.action_type.fillna("").str.upper())
        both = bool(types & AGO_T) and bool(types & ANT_T)
        if both:
            labs.add("SERM")
        rows.append({"chembl_parent": p, "label_b": top(labs), "b_n_records": len(g),
                     "b_action_types": ";".join(sorted(t for t in types if t)),
                     "b_mechanisms": " | ".join(sorted(set(g.mechanism_of_action.dropna()))),
                     "b_targets": ";".join(sorted(set(g.target_chembl_id))), "b_agonist_and_antagonist": both,
                     "b_record_labels": ";".join(sorted(labs))})
    return pd.DataFrame(rows, columns=["chembl_parent", "label_b", "b_n_records", "b_action_types", "b_mechanisms",
                                       "b_targets", "b_agonist_and_antagonist", "b_record_labels"])


def desc_flags(desc: str):
    d = (desc or "").lower()
    deg = bool(RE_DEG.search(d))
    ant = bool(RE_ANT.search(d))
    ago = (not ant) and bool(RE_AGO.search(d))
    return deg, ant, ago


def c_class(deg, ant, ago):
    if deg:
        return "SERD"
    if ago and ant:
        return "SERM"
    if ant:
        return "antagonist"
    if ago:
        return "agonist"
    return "unlabelled"


def eq(c):
    return "agonist" if c == "partial_agonist" else c


def main():
    loi = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv")
    chains = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv")
    lig = pd.read_csv(OUT / "ligand_ids.csv")
    edoi = pd.read_csv(OUT / "entry_doi.csv")
    mech = pd.read_csv(DOUT / "mechanism_records.csv")
    docs = pd.read_csv(DOUT / "doi_documents.csv")
    dact = pd.read_csv(DOUT / "document_esr1_activities.csv")
    teammate = TEAMMATE.exists()

    rows = (loi.dropna(subset=["ligand_id"]).groupby(["entry_id", "ligand_id"])
            .agg(n_chains_q18=("chain_key", "nunique")).reset_index().rename(columns={"ligand_id": "ccd"}))
    rows = rows.merge(edoi[["entry_id", "doi", "doi_source"]], on="entry_id", how="left")
    rows = rows.merge(lig[["ccd", "chembl_parent", "match_type", "parent_forms", "parent_pref_name"]], on="ccd",
                      how="left")

    # (b)
    rb = route_b(mech)
    rb.to_csv(DOUT / "route_b_per_parent.csv", index=False)
    rows = rows.merge(rb[["chembl_parent", "label_b", "b_action_types", "b_mechanisms"]], on="chembl_parent",
                      how="left")
    rows["label_b"] = rows.label_b.fillna("unlabelled")

    # (c)
    doc_by_doi = docs.dropna(subset=["document_chembl_id"]).groupby("doi").document_chembl_id.apply(set).to_dict()
    dact = dact.copy()
    dact["pchembl_value"] = pd.to_numeric(dact.pchembl_value, errors="coerce")
    ev, lc, reason, ndesc = [], [], [], []
    for r in rows.itertuples():
        if not isinstance(r.chembl_parent, str):
            lc.append("unlabelled"); reason.append("no_chembl_parent"); ndesc.append(0); continue
        if not isinstance(r.doi, str):
            lc.append("unlabelled"); reason.append("no_doi"); ndesc.append(0); continue
        dset = doc_by_doi.get(r.doi, set())
        if not dset:
            lc.append("unlabelled"); reason.append("doi_not_in_chembl"); ndesc.append(0); continue
        forms = set(r.parent_forms.split(";"))
        a = dact[dact.document_chembl_id.isin(dset) & dact.molecule_chembl_id.isin(forms)]
        if a.empty:
            lc.append("unlabelled"); reason.append("ligand_has_no_esr1_activity_in_document"); ndesc.append(0)
            continue
        a6 = a[a.pchembl_value >= PCHEMBL_MIN]
        if a6.empty:
            lc.append("unlabelled"); reason.append("no_activity_with_pchembl_ge6"); ndesc.append(0); continue
        descs = sorted(set(a6.assay_description.fillna("")))
        flags = [desc_flags(d) for d in descs]
        for d, (deg, ant, ago) in zip(descs, flags):
            sub = a6[a6.assay_description.fillna("") == d]
            ev.append({"entry_id": r.entry_id, "ccd": r.ccd, "chembl_parent": r.chembl_parent,
                       "document_chembl_ids": ";".join(sorted(set(sub.document_chembl_id))),
                       "assay_chembl_ids": ";".join(sorted(set(sub.assay_chembl_id))),
                       "n_activities": len(sub), "max_pchembl": sub.pchembl_value.max(),
                       "assay_description": d, "degradation": deg, "antagonism": ant, "agonism": ago})
        deg = any(f[0] for f in flags); ant = any(f[1] for f in flags); ago = any(f[2] for f in flags)
        cl = c_class(deg, ant, ago)
        lc.append(cl); reason.append("pattern_match" if cl != "unlabelled" else "no_pattern_match")
        ndesc.append(len(descs))
    rows["label_c"] = lc
    rows["c_reason"] = reason
    rows["c_n_descriptions"] = ndesc
    pd.DataFrame(ev, columns=["entry_id", "ccd", "chembl_parent", "document_chembl_ids", "assay_chembl_ids",
                              "n_activities", "max_pchembl", "assay_description", "degradation", "antagonism",
                              "agonism"]).to_csv(DOUT / "route_c_evidence.csv", index=False)

    # validation gate (c) vs (b)
    both = (rows.label_b != "unlabelled") & (rows.label_c != "unlabelled")
    rows["c_compared"] = both
    rows["c_agrees_b"] = [eq(b) == eq(c) if k else None for b, c, k in zip(rows.label_b, rows.label_c, both)]
    n_cmp = int(both.sum())
    n_agree = int(sum(1 for v in rows.c_agrees_b if v is True))
    agree = n_agree / n_cmp if n_cmp else float("nan")
    gate = bool(n_cmp >= GATE_MIN_ROWS and agree >= GATE_MIN_AGREE)
    val = rows.loc[:, ["entry_id", "ccd", "chembl_parent", "match_type", "doi", "label_b", "label_c", "c_reason",
                       "c_n_descriptions", "c_compared", "c_agrees_b"]]
    val.to_csv(OUT / "route_c_validation.csv", index=False)

    # final label
    if teammate:
        tm = pd.read_csv(TEAMMATE)
        rows = rows.merge(tm.rename(columns={"entry": "entry_id", "action_class": "label_a", "doi": "doi_a"}),
                          on=["entry_id", "ccd"], how="left")
    else:
        rows["label_a"] = None
    cls, src, conf = [], [], []
    for r in rows.itertuples():
        used = []
        if isinstance(r.label_a, str) and r.label_a:
            used.append(("teammate_file", r.label_a))
        if r.label_b != "unlabelled":
            used.append(("chembl_mechanism", r.label_b))
        if gate and r.label_c != "unlabelled":
            used.append(("entry_paper_assays", r.label_c))
        if not used:
            cls.append("unlabelled"); src.append("none"); conf.append(False); continue
        cls.append(used[0][1]); src.append(used[0][0])
        conf.append(len({eq(c) for _, c in used}) > 1)
    rows["action_class"] = cls
    rows["source"] = src
    rows["class_conflict"] = conf
    rows["route_c_used"] = gate
    rows["primary_set"] = rows.match_type != "connectivity_only"
    out_cols = ["entry_id", "ccd", "action_class", "source", "class_conflict", "doi", "match_type", "primary_set",
                "chembl_parent", "parent_pref_name", "label_b", "label_c", "route_c_used", "c_reason",
                "b_action_types", "b_mechanisms", "n_chains_q18", "doi_source"]
    ac = rows[out_cols].sort_values(["entry_id", "ccd"]).reset_index(drop=True)
    ac.to_csv(OUT / "action_class.csv", index=False)

    if teammate:
        ag = rows[["entry_id", "ccd", "label_a", "label_b", "label_c"]].copy()
        ag["agree_b"] = [eq(a) == eq(b) if isinstance(a, str) and b != "unlabelled" else None
                         for a, b in zip(ag.label_a, ag.label_b)]
        ag.to_csv(OUT / "class_agreement.csv", index=False)

    # chain-level labels and coverage
    ch = chains[["chain_key", "entry_id", "s08_cluster", "ligand_id_primary", "n_ligand_codes", "ligand_codes"]]
    ch = ch.rename(columns={"ligand_id_primary": "ccd"})
    ch["cluster"] = ch.s08_cluster.map({7: "cl7", 4: "cl4"}).fillna("minor")
    ch = ch.merge(ac[["entry_id", "ccd", "action_class", "source", "match_type", "primary_set", "class_conflict"]],
                  on=["entry_id", "ccd"], how="left")
    noli = ch.ccd.isna()
    ch.loc[noli, "action_class"] = "unlabelled"
    ch.loc[noli, "source"] = "no_ligand_of_interest"
    ch.loc[noli, "primary_set"] = True
    # classes of every ligand code of the chain (for chains with two codes)
    acmap = ac.set_index(["entry_id", "ccd"]).action_class.to_dict()
    ch["classes_all_codes"] = [";".join(sorted({acmap.get((e, c), "unlabelled") for c in str(codes).split(";")}))
                               if isinstance(codes, str) else None for e, codes in zip(ch.entry_id, ch.ligand_codes)]
    ch.to_csv(OUT / "chain_class.csv", index=False)

    cov = []
    for set_name, sub in (("all", ch), ("primary", ch[ch.primary_set.astype(bool)])):
        for cl in CLASSES + ["unlabelled"]:
            for srcname in ("any", "teammate_file", "chembl_mechanism", "entry_paper_assays", "none",
                            "no_ligand_of_interest"):
                s = sub[sub.action_class == cl]
                if srcname != "any":
                    s = s[s.source == srcname]
                if srcname != "any" and s.empty:
                    continue
                for clu in ("all", "cl7", "cl4", "minor"):
                    t = s if clu == "all" else s[s.cluster == clu]
                    cov.append({"set": set_name, "action_class": cl, "source": srcname, "cluster": clu,
                                "n_chains": len(t), "n_entries": t.entry_id.nunique(), "n_ligands": t.ccd.nunique(),
                                "entries_ge5": t.entry_id.nunique() >= COVERAGE_MIN_ENTRIES
                                if cl != "unlabelled" else None})
    cov = pd.DataFrame(cov)
    cov.to_csv(OUT / "class_coverage.csv", index=False)

    conflicts = ac[ac.class_conflict]
    num = {
        "teammate_file_exists": teammate,
        "n_rows": len(ac), "n_entries": int(ac.entry_id.nunique()), "n_ccd": int(ac.ccd.nunique()),
        "rows_by_match_type": ac.match_type.value_counts(dropna=False).to_dict(),
        "ccd_by_match_type": lig.match_type.value_counts(dropna=False).to_dict(),
        "label_b_rows": ac.label_b.value_counts().to_dict(),
        "label_b_parents": rb.label_b.value_counts().to_dict(),
        "label_c_rows": ac.label_c.value_counts().to_dict(),
        "c_reason_rows": ac.c_reason.value_counts().to_dict(),
        "gate": {"n_compared": n_cmp, "n_agree": n_agree, "agreement": agree, "passed": gate,
                 "n_compared_ccd": int(rows.loc[both, "ccd"].nunique()),
                 "n_compared_entries": int(rows.loc[both, "entry_id"].nunique()),
                 "min_rows": GATE_MIN_ROWS, "min_agreement": GATE_MIN_AGREE},
        "final_rows": ac.action_class.value_counts().to_dict(),
        "final_rows_by_source": ac.groupby(["source", "action_class"]).size().rename("n").reset_index()
        .to_dict("records"),
        "final_ccd": ac.drop_duplicates("ccd").action_class.value_counts().to_dict(),
        "n_conflicts": len(conflicts), "conflicts": conflicts[["entry_id", "ccd", "label_b", "label_c",
                                                               "action_class"]].to_dict("records"),
        "chains": {"n": len(ch), "by_class": ch.action_class.value_counts().to_dict(),
                   "by_class_cluster": ch.groupby(["action_class", "cluster"]).size().rename("n").reset_index()
                   .to_dict("records"),
                   "two_codes_differing_class": int(ch.classes_all_codes.fillna("").str.contains(";").sum())},
    }
    (OUT / "numbers.json").write_text(json.dumps(num, indent=1, default=str))
    print(json.dumps({k: num[k] for k in ("n_rows", "rows_by_match_type", "label_b_rows", "label_c_rows",
                                          "gate", "final_rows", "n_conflicts")}, indent=1, default=str))


if __name__ == "__main__":
    main()
