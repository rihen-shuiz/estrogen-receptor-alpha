"""SV12 step 1 - identifiers.

- Per CCD code of a ligand of interest (results/chain_ligand_of_interest.csv, Q18): stereo SMILES and InChIKey
  from the RCSB chemcomp record (rcsb_chem_comp_descriptor.SMILES_stereo / .InChIKey), cached in
  data/sv/rcsb_chemcomp/.
- Per entry: primary-citation DOI from the restored mmCIF (data/raw/<ID>.cif, _citation.pdbx_database_id_DOI where
  _citation.id = 'primary'); the Data API rcsb_primary_citation only for entries whose mmCIF has none (GraphQL,
  tested on 5TM1 first, cached in data/sv/api_cache/).
Writes data/sv/SV12/ccd_ids.csv (ChEMBL columns are added by SV12_chembl.py into results/sv/SV12/ligand_ids.csv)
and results/sv/SV12/entry_doi.csv.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import gemmi
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv/common"))
import cached_http as ch  # noqa: E402

RCSB_CC = ROOT / "data/sv/rcsb_chemcomp"
API = ROOT / "data/sv/api_cache"
OUT = ROOT / "results/sv/SV12"
DOUT = ROOT / "data/sv/SV12"
GQL = "https://data.rcsb.org/graphql"
QUERY = ("query ($ids: [String!]!) { entries(entry_ids: $ids) { rcsb_id "
         "rcsb_primary_citation { id pdbx_database_id_DOI pdbx_database_id_PubMed title year } } }")


def nan(v):
    return None if v in (None, "?", ".", "") else v


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOUT.mkdir(parents=True, exist_ok=True)
    loi = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv")
    ccds = sorted(loi.ligand_id.dropna().unique())

    rows = []
    for c in ccds:
        d, path = ch.fetch_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{c}", RCSB_CC, f"chemcomp_{c}",
                                allow_404=True)
        if d is None:
            rows.append({"ccd": c, "rcsb_status": "404"})
            continue
        desc = d.get("rcsb_chem_comp_descriptor") or {}
        cc = d.get("chem_comp") or {}
        rows.append({"ccd": c, "rcsb_status": "ok", "name": cc.get("name"), "formula": cc.get("formula"),
                     "formal_charge": cc.get("pdbx_formal_charge"), "smiles_stereo": desc.get("SMILES_stereo"),
                     "smiles": desc.get("SMILES"), "inchi": desc.get("InChI"), "inchikey": desc.get("InChIKey"),
                     "rcsb_cache": str(path.relative_to(ROOT))})
    ids = pd.DataFrame(rows)
    ids.to_csv(DOUT / "ccd_ids.csv", index=False)

    # primary-citation DOI per entry
    entries = sorted(loi.entry_id.unique())
    erows = []
    for e in entries:
        f = ROOT / "data/raw" / f"{e}.cif"
        doi = pmid = title = year = None
        n_primary = 0
        if f.exists():
            b = gemmi.cif.read(str(f)).sole_block()
            for r in b.find("_citation.", ["id", "?pdbx_database_id_DOI", "?pdbx_database_id_PubMed", "?title",
                                            "?year"]):
                if r.str(0) == "primary":
                    n_primary += 1
                    doi = nan(r.str(1)) if r.has(1) else None
                    pmid = nan(r.str(2)) if r.has(2) else None
                    title = nan(r.str(3)) if r.has(3) else None
                    year = nan(r.str(4)) if r.has(4) else None
        erows.append({"entry_id": e, "mmcif_found": f.exists(), "n_primary_citation_rows": n_primary,
                      "doi": doi.strip() if doi else None, "pubmed_id": pmid,
                      "citation_title": " ".join(title.split()) if title else None, "citation_year": year,
                      "doi_source": "mmcif_citation_primary" if doi else None})
    edf = pd.DataFrame(erows)

    # Data API fallback (rcsb_primary_citation), test on 5TM1 first
    test, _ = ch.fetch_json(GQL, API, "sv12_citation_test", body={"query": QUERY, "variables": {"ids": ["5TM1"]}})
    if "errors" in test:
        raise SystemExit(f"GraphQL test on 5TM1 failed: {test['errors']}")
    miss = edf.loc[edf.doi.isna(), "entry_id"].tolist()
    api = {}
    for i in range(0, len(miss), 50):
        batch = miss[i:i + 50]
        d, _ = ch.fetch_json(GQL, API, "sv12_citation", body={"query": QUERY, "variables": {"ids": batch}})
        if "errors" in d:
            raise SystemExit(f"GraphQL error: {d['errors']}")
        for ent in d["data"]["entries"]:
            api[ent["rcsb_id"]] = ent.get("rcsb_primary_citation") or {}
    for e in miss:
        pc = api.get(e, {})
        doi = nan(pc.get("pdbx_database_id_DOI"))
        m = edf.entry_id == e
        if doi:
            edf.loc[m, "doi"] = doi.strip()
            edf.loc[m, "doi_source"] = "data_api_rcsb_primary_citation"
        else:
            edf.loc[m, "doi_source"] = "none"
        for k_out, k_in in (("pubmed_id", "pdbx_database_id_PubMed"), ("citation_title", "title"),
                            ("citation_year", "year")):
            if pd.isna(edf.loc[m, k_out]).all() and pc.get(k_in) is not None:
                edf.loc[m, k_out] = str(pc.get(k_in))
    edf.to_csv(OUT / "entry_doi.csv", index=False)
    info = {"n_ccd": len(ccds), "rcsb_ok": int((ids.rcsb_status == "ok").sum()),
            "inchikey_missing": int(ids.inchikey.isna().sum()) if "inchikey" in ids else len(ccds),
            "smiles_stereo_missing": int(ids.smiles_stereo.isna().sum()) if "smiles_stereo" in ids else len(ccds),
            "n_entries": len(entries), "doi_mmcif": int((edf.doi_source == "mmcif_citation_primary").sum()),
            "doi_data_api": int((edf.doi_source == "data_api_rcsb_primary_citation").sum()),
            "doi_none": int((edf.doi_source == "none").sum()),
            "graphql_test_5TM1": test["data"]["entries"], "http": dict(ch.STATS)}
    (DOUT / "ids_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main()
