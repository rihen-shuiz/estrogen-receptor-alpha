"""SV00 step 6 (metadata part): entry metadata from the restored mmCIFs (gemmi), RCSB Data API gaps/cross-check and
the SIFTS UniProt-author numbering offset per kept chain. Rules: docs/sv/RUN_SV1028.md section 3.

Entries: the 351 entries holding the 788 kept chains (results/rmsd_entry_index.csv).
mmCIF (gemmi.cif): _exptl_crystal_grow.pH/.temp (crystal_id 1), _diffrn.ambient_temp (id 1), _exptl.method,
_symmetry.space_group_name_H-M, _refine.ls_d_res_high/.ls_R_factor_R_free/.ls_R_factor_R_work (first row),
TLS refined = category _pdbx_refine_tls exists, _pdbx_audit_revision_history (initial release, latest revision).
Cleaning: '?' and '.' -> NaN; several rows -> crystal / diffraction id 1 (else the first row) and the entry is flagged
multi_valued (the categories that caused it are listed in multi_valued_fields); [DEFAULT] plausibility: growth pH 3-11,
growth temperature 270-315 K, ambient temperature 50-320 K, values outside -> NaN and counted.
Data API: the Appendix B GraphQL query (tested once on 5TM1 by SV00_api_test.py), batches of <= 50 entries, raw JSON
cached in data/sv/api_cache/; used only to fill a field the mmCIF lacks (same id-1 rule and plausibility ranges) and as
a cross-check of the fields both report.
Offsets: _pdbx_sifts_unp_segments is absent from the RCSB files, so the SIFTS fallback of section 3 is used:
https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{id} (JSON cached in data/sv/sifts_cache/). Per kept chain
(struct_asym_id = label_asym_id) and P03372 segment overlapping UniProt 305-550: offset = unp_start - author start
residue number (checked against unp_end - author end). Cross-checks: the deposited _struct_ref_seq segment
(db_align_beg - pdbx_auth_seq_align_beg) and the existing residue mapping (esr1_position - auth_seq_id).
Politeness: <= 4 requests/s, retries with backoff, every response cached and reused.
Outputs: results/sv/SV00/entry_metadata_mmcif.csv, results/sv/SV00/chain_offsets.csv,
results/sv/SV00/metadata_cleaning.csv, results/sv/SV00/api_crosscheck.csv.
Run: .venv/bin/python -I scripts/sv/SV00_metadata.py
"""
import json
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV00"
API_CACHE = ROOT / "data/sv/api_cache"
SIFTS_CACHE = ROOT / "data/sv/sifts_cache"
GRAPHQL = "https://data.rcsb.org/graphql"
SIFTS_URL = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{id}"
UNIPROT = "P03372"
LBD = (305, 550)
PLAUSIBLE = {"ph": (3.0, 11.0), "growth_temp_K": (270.0, 315.0), "ambient_temp_K": (50.0, 320.0)}  # [DEFAULT]
MIN_INTERVAL_S = 0.25
_last = [0.0]

API_QUERY = """query ($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    exptl { method }
    exptl_crystal_grow { crystal_id pH temp method }
    diffrn { id crystal_id ambient_temp }
    symmetry { space_group_name_H_M }
    cell { length_a length_b length_c angle_alpha angle_beta angle_gamma }
    refine { ls_d_res_high ls_R_factor_R_work ls_R_factor_R_free }
    rcsb_accession_info { deposit_date initial_release_date }
    pdbx_audit_revision_history { ordinal revision_date major_revision minor_revision }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids uniprot_ids }
      entity_poly { rcsb_mutation_count pdbx_seq_one_letter_code_can }
      rcsb_polymer_entity_align {
        reference_database_name reference_database_accession
        aligned_regions { entity_beg_seq_id ref_beg_seq_id length }
      }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          auth_asym_id auth_to_entity_poly_seq_mapping
        }
      }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name formula_weight } }
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers { auth_asym_id auth_seq_id }
        rcsb_nonpolymer_instance_validation_score { RSCC is_subject_of_investigation }
      }
    }
  }
}"""


def polite(method, url, **kw):
    for attempt in range(1, 5):
        wait = MIN_INTERVAL_S - (time.monotonic() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.monotonic()
        try:
            r = requests.request(method, url, timeout=90, **kw)
            if r.status_code == 404:
                return r
            r.raise_for_status()
            return r
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(2 * attempt)


def post_graphql(ids):
    r = polite("POST", GRAPHQL, json={"query": API_QUERY, "variables": {"ids": list(ids)}})
    return r.json()


def num(value):
    if value is None or value in ("?", "."):
        return np.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return "unparseable"


def pick(rows, id_field):
    """Row with id '1' (else the first) and whether there were several rows."""
    if not rows:
        return {}, False
    chosen = next((r for r in rows if r.get(id_field) == "1"), rows[0])
    return chosen, len(rows) > 1


def table_opt(block, cat, fields):
    present = []
    for f in fields:
        if block.find_values(cat + f):
            present.append(f)
    if not present:
        return []
    t = block.find(cat, present)
    return [{f: gemmi.cif.as_string(row[i]) if row[i] not in ("?", ".") else None for i, f in enumerate(present)}
            for row in t]


def mmcif_entry(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    grow, m1 = pick(table_opt(block, "_exptl_crystal_grow.", ["crystal_id", "pH", "temp"]), "crystal_id")
    diff, m2 = pick(table_opt(block, "_diffrn.", ["id", "crystal_id", "ambient_temp"]), "id")
    refs = table_opt(block, "_refine.", ["pdbx_refine_id", "ls_d_res_high", "ls_R_factor_R_free", "ls_R_factor_R_work"])
    ref, m3 = (refs[0], len(refs) > 1) if refs else ({}, False)
    sg = block.find_value("_symmetry.space_group_name_H-M")
    method = block.find_values("_exptl.method")
    rev = table_opt(block, "_pdbx_audit_revision_history.", ["ordinal", "revision_date"])
    dates = sorted(r["revision_date"] for r in rev if r.get("revision_date"))
    initial = next((r["revision_date"] for r in rev if r.get("ordinal") == "1"), dates[0] if dates else None)
    multi = [n for n, m in (("exptl_crystal_grow", m1), ("diffrn", m2), ("refine", m3)) if m]
    return {
        "exp_method": ";".join(gemmi.cif.as_string(v) for v in method) if method else None,
        "ph_raw": grow.get("pH"), "growth_temp_raw": grow.get("temp"), "ambient_temp_raw": diff.get("ambient_temp"),
        "space_group_mmcif": gemmi.cif.as_string(sg) if sg not in (None, "?", ".") else None,
        "resolution_mmcif": num(ref.get("ls_d_res_high")), "r_free_mmcif": num(ref.get("ls_R_factor_R_free")),
        "r_work_mmcif": num(ref.get("ls_R_factor_R_work")),
        "tls_refined": bool(block.find_mmcif_category("_pdbx_refine_tls.")),
        "initial_release_date": initial, "latest_revision_date": dates[-1] if dates else None,
        "n_revisions": len(rev), "n_crystal_grow_rows": len(table_opt(block, "_exptl_crystal_grow.", ["crystal_id"])),
        "n_diffrn_rows": len(table_opt(block, "_diffrn.", ["id"])), "n_refine_rows": len(refs),
        "multi_valued": bool(multi), "multi_valued_fields": ";".join(multi),
        "poly_seq_scheme": {(r["asym_id"], int(r["seq_id"])): (r.get("pdb_seq_num"), r.get("pdb_ins_code"))
                            for r in table_opt(block, "_pdbx_poly_seq_scheme.",
                                               ["asym_id", "seq_id", "pdb_seq_num", "pdb_ins_code"])},
        "struct_ref_seq": table_opt(block, "_struct_ref_seq.", ["pdbx_strand_id", "pdbx_db_accession", "db_align_beg",
                                                              "db_align_end", "pdbx_auth_seq_align_beg",
                                                              "pdbx_auth_seq_align_end"]),
    }


def api_entries(ids):
    API_CACHE.mkdir(parents=True, exist_ok=True)
    out = {}
    for b in range(0, len(ids), 50):
        batch = ids[b:b + 50]
        path = API_CACHE / f"batch_{b // 50:02d}.json"
        if path.exists() and json.loads(path.read_text()).get("ids") == batch:
            res = json.loads(path.read_text())["response"]
        else:
            res = post_graphql(batch)
            path.write_text(json.dumps({"ids": batch, "query": API_QUERY, "response": res}))
        if res.get("errors"):
            raise RuntimeError(f"Data API errors in batch {b // 50}: {res['errors'][:2]}")
        for e in res["data"]["entries"]:
            out[e["rcsb_id"].upper()] = e
    return out


def api_values(e):
    grow, _ = pick([{k: (str(v) if v is not None else None) for k, v in r.items()} for r in (e.get("exptl_crystal_grow") or [])], "crystal_id")
    diff, _ = pick([{k: (str(v) if v is not None else None) for k, v in r.items()} for r in (e.get("diffrn") or [])], "id")
    ref = (e.get("refine") or [{}])[0]
    return {"ph": num(grow.get("pH")), "growth_temp_K": num(grow.get("temp")), "ambient_temp_K": num(diff.get("ambient_temp")),
            "space_group": (e.get("symmetry") or {}).get("space_group_name_H_M"),
            "resolution": num(ref.get("ls_d_res_high")), "r_free": num(ref.get("ls_R_factor_R_free")),
            "r_work": num(ref.get("ls_R_factor_R_work"))}


def sifts(pdb):
    SIFTS_CACHE.mkdir(parents=True, exist_ok=True)
    path = SIFTS_CACHE / f"{pdb}.json"
    if not path.exists():
        r = polite("GET", SIFTS_URL.format(id=pdb.lower()))
        path.write_text(r.text if r.status_code == 200 else json.dumps({"_status": r.status_code}))
    return json.loads(path.read_text())


def sifts_offsets(pdb, data, scheme):
    """{label_asym_id: per-residue offsets UniProt - author number over UniProt 305-550}.

    SIFTS gives label_seq_id -> UniProt per contiguous segment (start/end residue_number, unp_start/unp_end); the
    author number of each label_seq_id comes from the file's _pdbx_poly_seq_scheme.pdb_seq_num (residues with an
    insertion code are skipped). Also returns the number of segments whose label and UniProt lengths disagree."""
    res, bad = {}, 0
    maps = data.get(pdb.lower(), {}).get("UniProt", {}).get(UNIPROT, {}).get("mappings", [])
    for m in maps:
        if m["unp_end"] < LBD[0] or m["unp_start"] > LBD[1]:
            continue
        s_lab, e_lab = m["start"]["residue_number"], m["end"]["residue_number"]
        if e_lab - s_lab != m["unp_end"] - m["unp_start"]:
            bad += 1
            continue
        asym = m["struct_asym_id"]
        for lab in range(s_lab, e_lab + 1):
            unp = m["unp_start"] + lab - s_lab
            auth, ins = scheme.get((asym, lab), (None, None))
            if not (LBD[0] <= unp <= LBD[1]) or auth is None or ins or not str(auth).lstrip("-").isdigit():
                continue
            res.setdefault(asym, []).append(unp - int(auth))
    return res, bad


def summarise_offsets(values):
    vals = [v for v in values if v is not None]
    if not vals:
        return np.nan, None
    uniq = sorted(set(vals))
    return float(max(uniq, key=vals.count)), len(uniq) == 1


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    entries = sorted(idx.entry_id.unique())
    lbd = pd.read_csv(ROOT / "data/sv/lbd_chains.csv", keep_default_na=False, na_values=[""], low_memory=False)
    lbd = lbd.set_index("chain_key")

    rows, srs, schemes = [], {}, {}
    for pdb in entries:
        r = mmcif_entry(ROOT / "data/raw" / f"{pdb}.cif")
        srs[pdb] = r.pop("struct_ref_seq")
        schemes[pdb] = r.pop("poly_seq_scheme")
        r["entry_id"] = pdb
        rows.append(r)
    meta = pd.DataFrame(rows)

    api = api_entries(entries)
    cleaning, cross = [], []
    for field, raw in (("ph", "ph_raw"), ("growth_temp_K", "growth_temp_raw"), ("ambient_temp_K", "ambient_temp_raw")):
        vals = meta[raw].map(num)
        n_missing_raw = int(vals.map(lambda v: isinstance(v, float) and np.isnan(v)).sum())
        unparse = vals.map(lambda v: v == "unparseable")
        vals = vals.where(~unparse, np.nan).astype(float)
        lo, hi = PLAUSIBLE[field]
        out_rng = vals.notna() & ~vals.between(lo, hi)
        out_ids = meta.entry_id[out_rng].tolist()
        out_vals = vals[out_rng].tolist()
        vals = vals.where(~out_rng, np.nan)
        api_v = meta.entry_id.map(lambda e: api_values(api[e])[field] if e in api else np.nan)
        api_v = pd.to_numeric(api_v.where(api_v != "unparseable", np.nan), errors="coerce")
        api_ok = api_v.between(lo, hi)
        fill = vals.isna() & api_ok & ~out_rng & ~unparse
        both = vals.notna() & api_v.notna()
        cross.append({"field": field, "n_both": int(both.sum()),
                      "n_agree": int((np.abs(vals[both] - api_v[both]) < 1e-6).sum())})
        meta[field] = vals.where(~fill, api_v)
        meta[field + "_source"] = np.where(fill, "data_api", np.where(vals.notna(), "mmcif", "missing"))
        cleaning.append({"field": field, "plausible_range": f"{lo}-{hi}", "n_entries": len(meta),
                         "n_missing_in_mmcif": n_missing_raw, "n_unparseable": int(unparse.sum()),
                         "n_out_of_range_set_nan": int(out_rng.sum()),
                         "out_of_range_entries": ";".join(f"{e}={v:g}" for e, v in zip(out_ids, out_vals)),
                         "n_filled_from_api": int(fill.sum()), "n_missing_final": int(meta[field].isna().sum())})
    for field, col in (("resolution", "resolution_mmcif"), ("r_free", "r_free_mmcif"), ("r_work", "r_work_mmcif"),
                       ("space_group", "space_group_mmcif")):
        a = meta[col]
        b = meta.entry_id.map(lambda e: api_values(api[e])[field] if e in api else np.nan)
        if field == "space_group":
            both = a.notna() & b.notna()
            agree = (a[both].astype(str).str.replace(" ", "") == b[both].astype(str).str.replace(" ", "")).sum()
        else:
            a = pd.to_numeric(a, errors="coerce")
            b = pd.to_numeric(b, errors="coerce")
            both = a.notna() & b.notna()
            agree = (np.abs(a[both] - b[both]) < 1e-6).sum()
        cross.append({"field": field, "n_both": int(both.sum()), "n_agree": int(agree)})
        if field == "r_work":
            fill = a.isna() & b.notna()
            meta["r_work"] = a.where(~fill, b)
            meta["r_work_source"] = np.where(fill, "data_api", np.where(a.notna(), "mmcif", "missing"))
            cleaning.append({"field": "r_work", "plausible_range": "", "n_entries": len(meta),
                             "n_missing_in_mmcif": int(a.isna().sum()), "n_unparseable": 0,
                             "n_out_of_range_set_nan": 0, "out_of_range_entries": "",
                             "n_filled_from_api": int(fill.sum()), "n_missing_final": int(meta.r_work.isna().sum())})
    cleaning.append({"field": "multi_valued", "plausible_range": "", "n_entries": len(meta),
                     "n_missing_in_mmcif": 0, "n_unparseable": 0, "n_out_of_range_set_nan": 0,
                     "out_of_range_entries": ";".join(f"{e}:{f}" for e, f in zip(meta.entry_id[meta.multi_valued],
                                                                                 meta.multi_valued_fields[meta.multi_valued])),
                     "n_filled_from_api": 0, "n_missing_final": int(meta.multi_valued.sum())})
    meta.to_csv(OUT / "entry_metadata_mmcif.csv", index=False)
    pd.DataFrame(cleaning).to_csv(OUT / "metadata_cleaning.csv", index=False)

    # offsets per kept chain
    mp = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_residue", "auth_seq_id"])
    mp = mp[mp.chain_key.isin(idx.chain_key) & mp.has_residue.astype(bool) & mp.esr1_position.between(*LBD)]
    mp = mp[pd.to_numeric(mp.auth_seq_id, errors="coerce").notna()]
    mp_off = (mp.esr1_position - mp.auth_seq_id.astype(float)).groupby(mp.chain_key).agg(
        lambda s: (float(s.mode().iloc[0]), s.nunique() == 1))
    orows = []
    for key in idx.chain_key:
        pdb, label, _ = key.split("|")
        auth = lbd.loc[key, "auth_asym_id"]
        data = sifts(pdb)
        so_all, n_bad = sifts_offsets(pdb, data, schemes[pdb])
        so = so_all.get(label, [])
        off, const = summarise_offsets(so)
        segs = [s for s in srs[pdb] if s.get("pdbx_strand_id") == auth and s.get("pdbx_db_accession") in (UNIPROT, "ESR1_HUMAN")]
        srv = []
        for s in segs:
            try:
                b, e = int(s["db_align_beg"]), int(s["db_align_end"])
                if e < LBD[0] or b > LBD[1]:
                    continue
                srv.append(b - int(s["pdbx_auth_seq_align_beg"]))
            except (TypeError, ValueError, KeyError):
                continue
        sr_off, sr_const = summarise_offsets(srv)
        m_off, m_const = mp_off.get(key, (np.nan, None))
        source = "sifts_pdbe_api" if np.isfinite(off) else ("struct_ref_seq" if np.isfinite(sr_off) else "missing")
        orows.append({"chain_key": key, "entry_id": pdb, "label_asym_id": label, "auth_asym_id": auth,
                      "offset_sifts": off, "offset_sifts_constant": const, "n_sifts_residues": len(so),
                      "n_sifts_segments_skipped_entry": n_bad,
                      "offset_struct_ref_seq": sr_off, "offset_struct_ref_seq_constant": sr_const,
                      "offset_existing_mapping": m_off, "offset_existing_mapping_constant": m_const,
                      "uniprot_author_offset": off if np.isfinite(off) else sr_off, "offset_source": source})
    offs = pd.DataFrame(orows)
    offs.to_csv(OUT / "chain_offsets.csv", index=False)
    f = offs.offset_sifts.notna()
    cross.append({"field": "offset_sifts_vs_struct_ref_seq", "n_both": int((f & offs.offset_struct_ref_seq.notna()).sum()),
                  "n_agree": int((offs.offset_sifts == offs.offset_struct_ref_seq).sum())})
    cross.append({"field": "offset_sifts_vs_existing_mapping", "n_both": int((f & offs.offset_existing_mapping.notna()).sum()),
                  "n_agree": int((offs.offset_sifts == offs.offset_existing_mapping).sum())})
    pd.DataFrame(cross).to_csv(OUT / "api_crosscheck.csv", index=False)
    print(pd.DataFrame(cleaning).drop(columns=["out_of_range_entries"]).to_string(index=False))
    print(pd.DataFrame(cross).to_string(index=False))
    print(offs.offset_source.value_counts().to_dict(), offs.uniprot_author_offset.value_counts(dropna=False).to_dict())


if __name__ == "__main__":
    main()
