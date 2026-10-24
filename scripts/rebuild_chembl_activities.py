"""Rebuild tracker X13: full ChEMBL activity download for ESR1 (target CHEMBL206).

Inputs: ChEMBL REST API (https://www.ebi.ac.uk/chembl/api/data/activity.json), network.
Method: filter target_chembl_id=CHEMBL206, standard_type in {IC50, Ki, Kd}, standard_relation '=';
follow page_meta.next across ALL pages (limit 1000 per page) and keep rows whose relation is also '='
(the frozen file has relation == standard_relation == '=' for every row). Rows are sorted by activity_id.
Blank canonical_smiles are kept blank. The download stops after a 15 minute budget and then writes
nothing final (a partial file would look complete), only a note in the log.
Output: results/rebuild/chembl_esr1_activities.csv with the frozen column set; frozen file is
data/chembl_esr1_activities.csv (never overwritten). The live ChEMBL release may differ from the frozen one.
Run from the repo root: PYTHONPATH=. python scripts/rebuild_chembl_activities.py
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import pandas as pd
import requests

from rebuild_common import OUT

BASE = "https://www.ebi.ac.uk"
FIRST = ("/chembl/api/data/activity.json?target_chembl_id=CHEMBL206&standard_type__in=IC50,Ki,Kd"
         "&standard_relation__exact=%3D&limit=1000&offset=0")
COLUMNS = ["activity_id", "assay_chembl_id", "molecule_chembl_id", "pchembl_value", "relation",
           "standard_relation", "standard_type", "standard_units", "standard_value", "target_chembl_id",
           "type", "units", "value", "canonical_smiles"]
BUDGET_S = 15 * 60


def get(url):
    for attempt in range(5):
        try:
            response = requests.get(url, timeout=120)
            response.raise_for_status()
            return response.json()
        except Exception:
            if attempt == 4:
                raise
            time.sleep(3 * (attempt + 1))


def main():
    start, url, rows, total = time.time(), BASE + FIRST, [], None
    while url:
        if time.time() - start > BUDGET_S:
            print(f"ChEMBL: budget exceeded after {len(rows)} of {total} rows; no output written")
            (OUT / "chembl_download_incomplete.txt").write_text(f"stopped after {len(rows)} of {total} rows\n")
            sys.exit(2)
        page = get(url)
        total = page["page_meta"]["total_count"]
        rows.extend(page["activities"])
        nxt = page["page_meta"]["next"]
        url = BASE + nxt if nxt else None
        print(f"  {len(rows)}/{total}", flush=True)
    table = pd.DataFrame(rows)
    for column in COLUMNS:
        if column not in table:
            table[column] = ""
    table = table[COLUMNS]
    table = table[(table.relation == "=") & (table.standard_relation == "=")]
    table = table.drop_duplicates("activity_id").sort_values("activity_id")
    table.to_csv(OUT / "chembl_esr1_activities.csv", index=False)
    (OUT / "chembl_download_incomplete.txt").unlink(missing_ok=True)
    print(f"ChEMBL: {len(table)} activities kept (API total_count {total}); "
          f"{table.molecule_chembl_id.nunique()} molecules, {int(table.canonical_smiles.isna().sum())} rows without SMILES")


if __name__ == "__main__":
    main()
