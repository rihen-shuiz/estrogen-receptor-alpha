"""SV12 §3 checks before any bulk ChEMBL use.

1. EST test: InChIKey from the RCSB chemcomp record of EST -> ChEMBL molecule (filter
   molecule_structures__standard_inchi_key) -> molecule_hierarchy.parent_chembl_id.
2. CHEMBL206: confirm it is the SINGLE PROTEIN target whose component is UniProt P03372.
3. Every ChEMBL target containing P03372 (target_components__accession=P03372): the 'ESR1 target set' used by
   route (b) (CHEMBL206 or a family/complex target containing it) and route (c).
Writes results/sv/SV12/chembl_checks.json and data/sv/SV12/esr1_targets.csv. Responses cached (§3).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv/common"))
import cached_http as ch  # noqa: E402

RCSB_CC = ROOT / "data/sv/rcsb_chemcomp"
CHEMBL = ROOT / "data/sv/chembl_cache"
OUT = ROOT / "results/sv/SV12"
DOUT = ROOT / "data/sv/SV12"


def chemcomp(ccd: str):
    return ch.fetch_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}", RCSB_CC, f"chemcomp_{ccd}",
                         allow_404=True)[0]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DOUT.mkdir(parents=True, exist_ok=True)
    checks = {}
    est = chemcomp("EST")
    key = est["rcsb_chem_comp_descriptor"]["InChIKey"]
    checks["EST_inchikey_rcsb"] = key
    data, path = ch.fetch_json(ch.chembl_url("molecule", {"molecule_structures__standard_inchi_key": key}),
                               CHEMBL, "molecule")
    mols = data["molecules"]
    checks["EST_molecule_hits"] = [m["molecule_chembl_id"] for m in mols]
    checks["EST_molecule_hierarchy"] = [m.get("molecule_hierarchy") for m in mols]
    checks["EST_pref_name"] = [m.get("pref_name") for m in mols]
    checks["EST_standard_inchi_key_returned"] = [(m.get("molecule_structures") or {}).get("standard_inchi_key")
                                                for m in mols]
    checks["EST_response_cache"] = str(path.relative_to(ROOT))

    t, path = ch.fetch_json(ch.chembl_url("target/CHEMBL206"), CHEMBL, "target")
    comps = [(c.get("accession"), c.get("component_type")) for c in t.get("target_components", [])]
    checks["CHEMBL206"] = {"pref_name": t.get("pref_name"), "target_type": t.get("target_type"),
                           "organism": t.get("organism"), "components": comps,
                           "cache": str(path.relative_to(ROOT))}
    checks["CHEMBL206_is_single_protein_P03372"] = bool(
        t.get("target_type") == "SINGLE PROTEIN" and comps == [("P03372", "PROTEIN")])

    targets = ch.chembl_get_all("target", {"target_components__accession": "P03372"}, CHEMBL, "target_list")
    rows = []
    for tg in targets:
        acc = sorted({c.get("accession") for c in tg.get("target_components", []) if c.get("accession")})
        rows.append({"target_chembl_id": tg["target_chembl_id"], "pref_name": tg.get("pref_name"),
                     "target_type": tg.get("target_type"), "organism": tg.get("organism"),
                     "n_components": len(tg.get("target_components", [])), "accessions": ";".join(acc)})
    tdf = pd.DataFrame(rows).sort_values("target_chembl_id")
    tdf.to_csv(DOUT / "esr1_targets.csv", index=False)
    checks["targets_containing_P03372"] = tdf.to_dict("records")
    single = tdf[(tdf.target_type == "SINGLE PROTEIN") & (tdf.accessions == "P03372")]
    checks["single_protein_targets_for_P03372"] = single.target_chembl_id.tolist()
    checks["http"] = dict(ch.STATS)
    (OUT / "chembl_checks.json").write_text(json.dumps(checks, indent=1))
    print(json.dumps(checks, indent=1))


if __name__ == "__main__":
    main()
