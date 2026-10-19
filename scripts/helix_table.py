"""Tracker K02: helix table for the ERalpha LBD, from the mmCIF helix annotation of the reference chain.

Reference chain: read from docs/reference.md (Chain key). Source: struct_conf rows with conf_type_id HELX_P
in data/raw/7NFB.cif for that label chain. DSSP is used only if a binary is on PATH (none was found when
this was written; the source column says which one was used).
Canonical ESR1 numbering: residue_mapping.csv (esr1_position) for the reference chain.
Names: the source gives helix IDs only (HELX_Pnn). A helix is named H12 only when it covers at least half of the
settings h12_range (analysis_settings.json); every other helix is name "unknown". Both H12 ranges are reported.
The settings are never changed here.
Output: results/helix_table.csv
Run from the repo root: python scripts/helix_table.py
"""
import json
import re
import shutil
from pathlib import Path

import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[1]
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
LBD = (SETTINGS["esr1_lbd_start"], SETTINGS["esr1_lbd_end"])
H12_SETTINGS = tuple(SETTINGS["h12_range"])
OUT = ROOT / "results/helix_table.csv"


def reference_chain_key():
    text = (ROOT / "docs/reference.md").read_text()
    match = re.search(r"Chain key `([^`]+)`", text)
    if not match:
        raise ValueError("docs/reference.md has no 'Chain key' line")
    return match.group(1)


def main():
    key = reference_chain_key()
    pdb_id, label_chain, _model = key.split("|")
    cif = MMCIF2Dict(str(ROOT / f"data/raw/{pdb_id}.cif"))
    fields = [name.split(".")[1] for name in cif if name.startswith("_struct_conf.")]
    count = len(cif["_struct_conf.id"])
    rows = [{f: cif[f"_struct_conf.{f}"][i] for f in fields} for i in range(count)]
    helices = [r for r in rows if r["conf_type_id"] == "HELX_P" and r["beg_label_asym_id"] == label_chain]

    mapping = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", dtype=str, keep_default_na=False)
    mapping = mapping[mapping.chain_key == key]
    by_label = {int(r.label_seq_id): (int(float(r.esr1_position)), r.has_ca == "True") for r in mapping.itertuples()
                if r.esr1_position != "" and r.label_seq_id != ""}

    dssp = shutil.which("mkdssp") or shutil.which("dssp")
    source = f"mmCIF struct_conf HELX_P (data/raw/{pdb_id}.cif, chain {label_chain}); " + (
        f"DSSP {dssp}" if dssp else "DSSP binary not found")

    h12_overlap = {}
    records = []
    for helix in helices:
        beg, end = int(helix["beg_label_seq_id"]), int(helix["end_label_seq_id"])
        start, end_e1 = by_label[beg][0], by_label[end][0]
        if end_e1 < LBD[0] or start > LBD[1]:
            continue
        span = range(start, end_e1 + 1)
        ca = sum(1 for label in range(beg, end + 1) if label in by_label and by_label[label][1])
        overlap = len(set(span) & set(range(H12_SETTINGS[0], H12_SETTINGS[1] + 1)))
        h12_overlap[helix["pdbx_PDB_helix_id"]] = overlap
        records.append({
            "helix": helix["pdbx_PDB_helix_id"], "name": "unknown", "name_status": "unknown_source_gives_no_name",
            "start": start, "end": end_e1, "length": end_e1 - start + 1, "ca_observed": ca,
            "source_label_seq_start": beg, "source_label_seq_end": end,
            "source_auth_start": helix["beg_auth_seq_id"], "source_auth_end": helix["end_auth_seq_id"],
            "helix_class": helix["pdbx_PDB_helix_class"], "source": source,
        })
    table = pd.DataFrame(records)
    named = max(h12_overlap, key=h12_overlap.get)
    if h12_overlap[named] * 2 < H12_SETTINGS[1] - H12_SETTINGS[0] + 1:
        raise ValueError("No helix covers half of the settings H12 range")
    row = table.helix == named
    table.loc[row, "name"] = "H12"
    table.loc[row, "name_status"] = "named_from_settings_h12_range"
    table["settings_h12_range"] = f"{H12_SETTINGS[0]}-{H12_SETTINGS[1]}"
    source_h12 = table.loc[row].iloc[0]
    table["h12_source_vs_settings"] = ""
    table.loc[row, "h12_source_vs_settings"] = (
        "agrees" if (source_h12.start, source_h12.end) == H12_SETTINGS
        else f"disagrees: source {source_h12.start}-{source_h12.end}, settings {H12_SETTINGS[0]}-{H12_SETTINGS[1]}; settings not changed")
    table.to_csv(OUT, index=False)
    print(table[["helix", "name", "start", "end", "length", "ca_observed", "source_auth_start", "source_auth_end", "helix_class"]].to_string(index=False))
    print(f"K02: {len(table)} helices, reference {key}; {source}")
    print(table.loc[row, "h12_source_vs_settings"].iloc[0])


if __name__ == "__main__":
    main()
