"""Rebuild tracker Q08 (results/ligands.csv) and C01 (results/smiles_parse_report.csv).

Inputs: results/rebuild/entries.csv (Q06 rebuild; ligand_ids per entry), CCD files in
esr_analysis/resources/data/ccd (local) or downloaded from https://files.rcsb.org/ligands/view/<code>.cif
into results/rebuild/ccd/ when no local copy exists.
Method (Q08): one row per distinct non-water CCD code over all entries, sorted by code;
esr_analysis.build_ligand_table.read_ccd picks the CCD SMILES (SMILES_CANONICAL, CACTVS first) and writes
the RDKit canonical isomeric SMILES; ccd_source says 'local esr_analysis/resources/data/ccd' or
'RCSB <url>'; n_entries/entry_ids = carrying entries (sorted, ';'-joined).
Method (C01): parse every CCD SMILES with RDKit; parsed, heavy-atom count, error text; failures are kept
as rows so no ligand silently drops out.
Outputs: results/rebuild/ligands.csv, results/rebuild/smiles_parse_report.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_ligands.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd
import requests
from rdkit import Chem, RDLogger

from esr_analysis.build_ligand_table import read_ccd
from rebuild_common import CCD_CACHE, CCD_LOCAL, OUT, read_csv

CCD_URL = "https://files.rcsb.org/ligands/view/{code}.cif"
COLUMNS = ["ligand_id", "ligand_name", "smiles", "canonical_isomeric_smiles", "smiles_source", "ccd_source",
           "ccd_sha256", "chemistry_status", "chemistry_error", "n_entries", "entry_ids"]


def ccd_file(code):
    local = CCD_LOCAL / f"{code}.cif"
    if local.exists():
        return local, "local esr_analysis/resources/data/ccd"
    cached = CCD_CACHE / f"{code}.cif"
    if not cached.exists():
        CCD_CACHE.mkdir(parents=True, exist_ok=True)
        response = requests.get(CCD_URL.format(code=code), timeout=60)
        response.raise_for_status()
        cached.write_bytes(response.content)
    return cached, f"RCSB {CCD_URL.format(code=code)}"


def main():
    RDLogger.DisableLog("rdApp.*")
    entries = read_csv(OUT / "entries.csv")
    carriers = {}
    for entry in entries.itertuples():
        for code in str(entry.ligand_ids if isinstance(entry.ligand_ids, str) else "").split(";"):
            if code:
                carriers.setdefault(code, []).append(entry.entry_id)
    rows, parse = [], []
    for code in sorted(carriers):
        path, source = ccd_file(code)
        chem, _ = read_ccd(path)
        rows.append({"ligand_id": code, "ligand_name": chem["ligand_name"], "smiles": chem["smiles"],
                     "canonical_isomeric_smiles": chem["canonical_isomeric_smiles"],
                     "smiles_source": chem["smiles_source"], "ccd_source": source, "ccd_sha256": chem["ccd_sha256"],
                     "chemistry_status": chem["chemistry_status"], "chemistry_error": chem["chemistry_error"],
                     "n_entries": len(carriers[code]), "entry_ids": ";".join(sorted(carriers[code]))})
        mol = Chem.MolFromSmiles(chem["smiles"]) if chem["smiles"] else None
        parse.append({"ligand_id": code, "parsed": mol is not None, "n_atoms": mol.GetNumAtoms() if mol else 0,
                      "error": "" if mol is not None else (chem["chemistry_error"] or "RDKit could not parse SMILES")})
    pd.DataFrame(rows, columns=COLUMNS).to_csv(OUT / "ligands.csv", index=False)
    report = pd.DataFrame(parse)
    report.to_csv(OUT / "smiles_parse_report.csv", index=False)
    print(f"Q08/C01 rebuild: {len(rows)} ligands; parsed {int(report.parsed.sum())}/{len(report)}")


if __name__ == "__main__":
    main()
