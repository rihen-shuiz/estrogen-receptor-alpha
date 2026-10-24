"""Rebuild tracker Q05/Q06: master entry table with QC flags (results/entries.csv).

Inputs: results/rebuild/lbd_extract/entries.csv (from scripts/rebuild_lbd_extract.py: resolution, R-free,
mean B and their sources, space group, deposit date, qc_status, declared_protein_chains, model_count, sha256)
and the mmCIFs (data/raw, else results/rebuild/raw) for depositing authors and ligand IDs.
Method: depositing_authors = _audit_author.name joined by '; '; ligand_ids = sorted distinct
_pdbx_entity_nonpoly.comp_id without HOH, joined by ';'. Flags per value (settings qc_resolution_A 2.5,
qc_r_free 0.30, qc_mean_B_A2 50): 'fail' if above the cut, 'unknown' if missing, else 'pass'.
qc_pass = combined status of esr_analysis.flag_quality.quality_status (any fail -> fail, else any missing ->
unknown, else pass), so a missing value is never a pass. No entry is dropped.
Output: results/rebuild/entries.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_entries_qc.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from esr_analysis.flag_quality import quality_status
from rebuild_common import OUT, SETTINGS, raw_path, read_csv

COLUMNS = ["entry_id", "resolution", "r_free", "mean_b", "space_group", "deposit_date", "depositing_authors",
           "ligand_ids", "resolution_source", "r_free_source", "B_source", "qc_status", "declared_protein_chains",
           "model_count", "source_sha256", "qc_res", "qc_rfree", "qc_b", "qc_pass"]


def flag(value, limit):
    if not np.isfinite(value):
        return "unknown"
    return "fail" if value > limit else "pass"


def main():
    lbd = read_csv(OUT / "lbd_extract/entries.csv")
    rows = []
    for entry in lbd.itertuples():
        data = MMCIF2Dict(str(raw_path(entry.pdb_id)))
        authors = "; ".join(data.get("_audit_author.name", []))
        ligands = sorted({code for code in data.get("_pdbx_entity_nonpoly.comp_id", []) if code != "HOH"})
        res, rfree, b = float(entry.resolution_A), float(entry.r_free), float(entry.mean_B_factor_A2)
        rows.append({"entry_id": entry.pdb_id, "resolution": res, "r_free": rfree, "mean_b": b,
                     "space_group": entry.space_group, "deposit_date": entry.deposit_date,
                     "depositing_authors": authors, "ligand_ids": ";".join(ligands),
                     "resolution_source": entry.resolution_source, "r_free_source": entry.r_free_source,
                     "B_source": entry.B_source, "qc_status": entry.qc_status,
                     "declared_protein_chains": entry.declared_protein_chains, "model_count": entry.model_count,
                     "source_sha256": entry.source_sha256,
                     "qc_res": flag(res, SETTINGS["qc_resolution_A"]), "qc_rfree": flag(rfree, SETTINGS["qc_r_free"]),
                     "qc_b": flag(b, SETTINGS["qc_mean_B_A2"]), "qc_pass": quality_status(res, rfree, b, SETTINGS)})
    table = pd.DataFrame(rows, columns=COLUMNS)
    table.to_csv(OUT / "entries.csv", index=False)
    print("Q06 rebuild:", len(table), "entries;", table.qc_pass.value_counts().to_dict())


if __name__ == "__main__":
    main()
