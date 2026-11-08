"""SV12 steps 2, 4b/4c (ChEMBL queries), 5 and 6 - ChEMBL match, mechanism records, entry documents, activities.

Run after SV12_chembl_test.py (EST test, CHEMBL206 check, ESR1 target set) and SV12_ids.py (InChIKeys, DOIs).
All ChEMBL responses are cached in data/sv/chembl_cache/ (every page, via page_meta.next).

Declared before the bulk run (rules the run file leaves open):
- Full match: molecule filter molecule_structures__standard_inchi_key = the RCSB InChIKey; parent =
  molecule_hierarchy.parent_chembl_id (the molecule itself when the hierarchy is empty).
- Connectivity match (no full match): filter molecule_structures__standard_inchi_key__startswith = first block.
  Several distinct parents -> the parent of a hit whose second block is UHFFFAOYSA (no stereo layer) if any, else
  the lowest numeric ChEMBL id; all candidates are listed in ligand_ids.csv.
- Several parents for one full InChIKey (not expected) -> lowest numeric ChEMBL id, flagged multi_parent.
- 'Forms' of a parent = molecule_form records with that parent_chembl_id (the parent included).
- Mechanism records: mechanism filter molecule_chembl_id__in = the forms; kept if target_chembl_id is in the ESR1
  target set (every ChEMBL target with component P03372: data/sv/SV12/esr1_targets.csv).
- Documents: document filter doi__iexact = the entry DOI (the plain 'doi' filter is case-sensitive; DOIs are not).
- Route (c) activities: activity filter document_chembl_id = each document, target_chembl_id__in = ESR1 target set;
  the ligand (forms of its parent) and pChEMBL >= 6 are filtered in SV12_classes.py.
- Activities for SV15: activity filter molecule_chembl_id__in = forms of the parent, target_chembl_id = CHEMBL206.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv/common"))
import cached_http as ch  # noqa: E402

C = ROOT / "data/sv/chembl_cache"
OUT = ROOT / "results/sv/SV12"
DOUT = ROOT / "data/sv/SV12"
ACT_KEEP = ["activity_id", "assay_chembl_id", "document_chembl_id", "assay_type", "assay_description",
            "standard_type", "standard_relation", "standard_value", "standard_units", "pchembl_value",
            "data_validity_comment", "data_validity_description", "potential_duplicate", "molecule_chembl_id",
            "parent_molecule_chembl_id", "target_chembl_id", "target_pref_name", "document_year",
            "document_journal", "src_id", "activity_comment", "assay_variant_mutation"]


def cid_num(c: str) -> int:
    return int(re.sub(r"\D", "", c))


def parent_of(m: dict) -> str:
    h = m.get("molecule_hierarchy") or {}
    return h.get("parent_chembl_id") or m["molecule_chembl_id"]


def key_of(m: dict):
    return (m.get("molecule_structures") or {}).get("standard_inchi_key")


def main():
    ids = pd.read_csv(DOUT / "ccd_ids.csv")
    targets = pd.read_csv(DOUT / "esr1_targets.csv")
    tset = set(targets.target_chembl_id)
    edoi = pd.read_csv(OUT / "entry_doi.csv")

    # --- step 2: match -------------------------------------------------------------------------------------------
    rows = []
    for r in ids.itertuples():
        key = r.inchikey
        full = ch.chembl_get_all("molecule", {"molecule_structures__standard_inchi_key": key}, C, "molecule")
        full = [m for m in full if key_of(m) == key]
        rec = {"ccd": r.ccd, "smiles_stereo": r.smiles_stereo, "inchikey": key, "inchikey_block1": key[:14]}
        if full:
            parents = sorted({parent_of(m) for m in full}, key=cid_num)
            rec.update(match_type="full", chembl_molecule_hits=";".join(m["molecule_chembl_id"] for m in full),
                       chembl_parent=parents[0], parent_candidates=";".join(parents),
                       multi_parent=len(parents) > 1)
        else:
            conn = ch.chembl_get_all("molecule", {"molecule_structures__standard_inchi_key__startswith": key[:14]},
                                     C, "molecule_conn")
            conn = [m for m in conn if (key_of(m) or "").startswith(key[:14] + "-")]
            if conn:
                parents = sorted({parent_of(m) for m in conn}, key=cid_num)
                nostereo = sorted({parent_of(m) for m in conn if key_of(m).split("-")[1] == "UHFFFAOYSA"},
                                  key=cid_num)
                chosen = nostereo[0] if nostereo else parents[0]
                rec.update(match_type="connectivity_only",
                           chembl_molecule_hits=";".join(f"{m['molecule_chembl_id']}:{key_of(m)}" for m in conn),
                           chembl_parent=chosen, parent_candidates=";".join(parents), multi_parent=len(parents) > 1,
                           connectivity_choice="no_stereo_layer" if nostereo else "lowest_chembl_id")
            else:
                rec.update(match_type="no_chembl", multi_parent=False)
        rows.append(rec)
    lig = pd.DataFrame(rows)

    parents = sorted(lig.chembl_parent.dropna().unique(), key=cid_num)
    pinfo, forms = {}, {}
    for p in parents:
        m, _ = ch.fetch_json(ch.chembl_url(f"molecule/{p}"), C, "molecule_parent")
        pinfo[p] = {"parent_pref_name": m.get("pref_name"), "parent_max_phase": m.get("max_phase"),
                    "parent_inchikey": key_of(m)}
        fr = ch.chembl_get_all("molecule_form", {"parent_chembl_id": p}, C, "molecule_form")
        forms[p] = sorted({f["molecule_chembl_id"] for f in fr} | {p}, key=cid_num)
    for col in ("parent_pref_name", "parent_max_phase", "parent_inchikey"):
        lig[col] = lig.chembl_parent.map(lambda p, c=col: pinfo.get(p, {}).get(c) if isinstance(p, str) else None)
    lig["parent_forms"] = lig.chembl_parent.map(lambda p: ";".join(forms[p]) if isinstance(p, str) else None)
    lig["is_e2_reference"] = lig.ccd == "EST"
    cols = ["ccd", "smiles_stereo", "inchikey", "chembl_parent", "match_type", "parent_pref_name",
            "parent_max_phase", "parent_inchikey", "parent_candidates", "multi_parent", "connectivity_choice",
            "chembl_molecule_hits", "parent_forms", "inchikey_block1", "is_e2_reference"]
    for c_ in cols:
        if c_ not in lig:
            lig[c_] = None
    lig[cols].to_csv(OUT / "ligand_ids.csv", index=False)

    # --- step 4b: mechanism records ------------------------------------------------------------------------------
    mrows = []
    for p in parents:
        recs = ch.chembl_get_all("mechanism", {"molecule_chembl_id__in": ",".join(forms[p])}, C, "mechanism")
        for m in recs:
            mrows.append({"chembl_parent": p, "molecule_chembl_id": m.get("molecule_chembl_id"),
                          "mec_id": m.get("mec_id"), "target_chembl_id": m.get("target_chembl_id"),
                          "esr1_target": m.get("target_chembl_id") in tset, "action_type": m.get("action_type"),
                          "mechanism_of_action": m.get("mechanism_of_action"),
                          "mechanism_comment": m.get("mechanism_comment"), "max_phase": m.get("max_phase"),
                          "direct_interaction": m.get("direct_interaction"),
                          "mechanism_refs": ";".join(f"{x.get('ref_type')}:{x.get('ref_id')}"
                                                     for x in (m.get("mechanism_refs") or []))})
    mech = pd.DataFrame(mrows, columns=["chembl_parent", "molecule_chembl_id", "mec_id", "target_chembl_id",
                                        "esr1_target", "action_type", "mechanism_of_action", "mechanism_comment",
                                        "max_phase", "direct_interaction", "mechanism_refs"])
    mech.to_csv(DOUT / "mechanism_records.csv", index=False)

    # --- step 4c: documents of the entry DOIs and their ESR1-target activities -----------------------------------
    drows = []
    for doi in sorted(edoi.doi.dropna().unique()):
        docs = ch.chembl_get_all("document", {"doi__iexact": doi}, C, "document")
        if not docs:
            drows.append({"doi": doi, "document_chembl_id": None})
        for d in docs:
            drows.append({"doi": doi, "document_chembl_id": d["document_chembl_id"], "chembl_doi": d.get("doi"),
                          "doc_type": d.get("doc_type"), "journal": d.get("journal"), "year": d.get("year"),
                          "pubmed_id": d.get("pubmed_id"), "src_id": d.get("src_id")})
    docs = pd.DataFrame(drows)
    docs.to_csv(DOUT / "doi_documents.csv", index=False)
    arows = []
    for dcid in sorted(docs.document_chembl_id.dropna().unique(), key=cid_num):
        acts = ch.chembl_get_all("activity", {"document_chembl_id": dcid,
                                              "target_chembl_id__in": ",".join(sorted(tset))}, C, "activity_doc")
        for a in acts:
            arows.append({k: a.get(k) for k in ACT_KEEP})
    dact = pd.DataFrame(arows, columns=ACT_KEEP)
    dact.to_csv(DOUT / "document_esr1_activities.csv", index=False)

    # --- step 5: all activities on CHEMBL206 for every matched parent --------------------------------------------
    p2ccd = lig.dropna(subset=["chembl_parent"]).groupby("chembl_parent").ccd.apply(lambda s: ";".join(sorted(s)))
    p2mt = lig.dropna(subset=["chembl_parent"]).groupby("chembl_parent").match_type.apply(
        lambda s: ";".join(sorted(set(s))))
    act = []
    for p in parents:
        recs = ch.chembl_get_all("activity", {"molecule_chembl_id__in": ",".join(forms[p]),
                                              "target_chembl_id": "CHEMBL206"}, C, "activity_206")
        for a in recs:
            row = {"chembl_parent": p, "ccd_codes": p2ccd[p], "match_types": p2mt[p]}
            row.update({k: a.get(k) for k in ACT_KEEP})
            act.append(row)
    adf = pd.DataFrame(act, columns=["chembl_parent", "ccd_codes", "match_types"] + ACT_KEEP)
    adf.to_csv(OUT / "chembl_activities.csv", index=False)

    info = {"n_ccd": len(lig), "match_type": lig.match_type.value_counts().to_dict(),
            "multi_parent": int(lig.multi_parent.fillna(False).astype(bool).sum()),
            "n_parents": len(parents), "n_mechanism_records": len(mech),
            "n_mechanism_records_esr1": int(mech.esr1_target.sum()) if len(mech) else 0,
            "n_dois": int(edoi.doi.nunique()), "n_dois_with_document": int(docs.dropna(
                subset=["document_chembl_id"]).doi.nunique()),
            "n_documents": int(docs.document_chembl_id.nunique()), "n_document_esr1_activities": len(dact),
            "n_activities_chembl206": len(adf), "n_parents_with_activities": int(adf.chembl_parent.nunique()),
            "EST": lig.loc[lig.ccd == "EST", ["chembl_parent", "match_type", "parent_pref_name"]].to_dict("records"),
            "http": dict(ch.STATS)}
    (DOUT / "chembl_info.json").write_text(json.dumps(info, indent=1, default=str))
    print(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
