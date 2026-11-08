"""SV00 step 1-2: restore the frozen mmCIF entry list into data/raw and write data/sv/manifest.csv.

Entry list: the frozen list of the original run, read with scripts/fetch_data.read_ids (data/pdb_ids_esr1.txt then
data/pdb_ids_esr2.txt; 487 + 39 = 526 IDs; tracker Q01/Q02 out, Q03 log ev data/raw/download_manifest.csv and
data/raw/sha256sums.txt). No RCSB search is run.
Download: esr_analysis.download_structures.download_one, the function scripts/fetch_data.py calls (same URL
https://files.rcsb.org/download/{ID}.cif, same plain .cif file format and validation, valid files already present are
kept). scripts/fetch_data.py itself is not run as a whole because its download_structures() call rewrites
data/raw/download_manifest.csv and its later steps write data/ccd and data/lbd; rule 1.2 allows SV00 to add only
missing files to data/raw and data/clean. The manifest of this step goes to data/sv/manifest.csv instead.
Required entries (those holding the 788 kept chains, results/entry_clusters.csv) are fetched first, then the rest.
Politeness: one request at a time, >= 0.25 s between requests (<= 4 req/s), download_one retries with backoff.
Each file is checked against data/raw/sha256sums.txt; the latest revision date comes from the file's
_pdbx_audit_revision_history (gemmi).
Outputs: data/raw/<ID>.cif (missing files only), data/sv/manifest.csv, data/sv/fetch.log.
Run: .venv/bin/python -I scripts/sv/SV00_fetch.py
"""
import datetime as dt
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import gemmi  # noqa: E402
import pandas as pd  # noqa: E402

from esr_analysis.download_structures import RCSB_CIF_URL, download_one  # noqa: E402
from esr_analysis.run_io import sha256  # noqa: E402
from fetch_data import read_ids  # noqa: E402

RAW = ROOT / "data/raw"
OUT = ROOT / "data/sv"
ORIGINAL_DOWNLOAD_DATE = "2026-10-14"  # Q03 DONE log line (tracker/esr1.log) recording data/raw/sha256sums.txt
MIN_INTERVAL_S = 0.25


def revision_info(path):
    block = gemmi.cif.read_file(str(path)).sole_block()
    table = block.find("_pdbx_audit_revision_history.", ["ordinal", "revision_date", "major_revision", "minor_revision"])
    rows = [(int(r[0]), r[1], r[2], r[3]) for r in table]
    if not rows:
        return {"latest_revision_date": "", "latest_revision": "", "n_revisions": 0}
    last = max(rows, key=lambda r: (r[1], r[0]))
    return {"latest_revision_date": last[1], "latest_revision": f"{last[2]}.{last[3]}", "n_revisions": len(rows)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids = read_ids()
    esr1 = set(x.strip().upper() for x in (ROOT / "data/pdb_ids_esr1.txt").read_text().split())
    required = set(pd.read_csv(ROOT / "results/entry_clusters.csv").entry_id.str.upper())
    frozen = {}
    for line in (RAW / "sha256sums.txt").read_text().splitlines():
        if line.strip():
            digest, name = line.split()
            frozen[name.replace(".cif", "").upper()] = digest
    order = [i for i in ids if i.upper() in required] + [i for i in ids if i.upper() not in required]
    log = open(OUT / "fetch.log", "a")
    log.write(f"# SV00 fetch {dt.datetime.now().isoformat(timespec='seconds')}: {len(order)} IDs, "
              f"{len(required)} required\n")
    rows, last_request = [], 0.0
    for n, pdb in enumerate(order, 1):
        pdb = pdb.upper()
        path = RAW / f"{pdb}.cif"
        existed = path.exists()
        if not existed:
            wait = MIN_INTERVAL_S - (time.monotonic() - last_request)
            if wait > 0:
                time.sleep(wait)
            last_request = time.monotonic()
        ok, status = download_one(pdb, path)
        if not existed and ok:
            status = "downloaded"
        row = {"entry": pdb, "receptor": "esr1" if pdb in esr1 else "esr2", "required": pdb in required,
               "path": path.relative_to(ROOT).as_posix(), "source_url": RCSB_CIF_URL.format(pdb_id=pdb),
               "download_ok": ok, "download_status": status,
               "download_date": dt.date.today().isoformat() if not existed else "",
               "sha256": "", "sha256_frozen": frozen.get(pdb, ""), "sha_matches_frozen": False}
        if ok:
            row["sha256"] = sha256(path)
            row["sha_matches_frozen"] = row["sha256"] == row["sha256_frozen"]
            row.update(revision_info(path))
        rows.append(row)
        log.write(f"{pdb},{ok},{status}\n")
        if n % 25 == 0:
            log.flush()
            print(f"  {n}/{len(order)}", flush=True)
    log.close()
    man = pd.DataFrame(rows)
    man["original_download_date"] = ORIGINAL_DOWNLOAD_DATE
    man["revised_after_original"] = man.latest_revision_date.fillna("").astype(str) > ORIGINAL_DOWNLOAD_DATE
    # an entry downloaded earlier in this step keeps its first download date
    old = OUT / "manifest.csv"
    if old.exists():
        prev = pd.read_csv(old, dtype=str, keep_default_na=False).set_index("entry")["download_date"]
        man["download_date"] = [d or prev.get(e, "") for e, d in zip(man.entry, man.download_date)]
    man.to_csv(OUT / "manifest.csv", index=False)
    print(f"SV00 fetch: {int(man.download_ok.sum())}/{len(man)} ok; required ok "
          f"{int(man[man.required].download_ok.sum())}/{int(man.required.sum())}; sha256 equal to frozen "
          f"{int(man.sha_matches_frozen.sum())}; revised after {ORIGINAL_DOWNLOAD_DATE}: "
          f"{int(man.revised_after_original.sum())}")


if __name__ == "__main__":
    main()
