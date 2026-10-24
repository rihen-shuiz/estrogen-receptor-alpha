"""Rebuild step 0 (input for tracker Q04/Q06/X01): make every listed mmCIF available.

Inputs: data/pdb_ids_esr1.txt + data/pdb_ids_esr2.txt (526 IDs), data/raw (local copies, git-ignored, may be
partial), data/raw/sha256sums.txt (hashes of the files the frozen results were built from).
Method: IDs with a local data/raw file are used as is; the rest are downloaded from RCSB with
esr_analysis.download_structures.download_one into results/rebuild/raw/ (data/raw is never written).
Every file is hashed and compared with data/raw/sha256sums.txt; RCSB re-releases can change a file.
Outputs: results/rebuild/raw/*.cif, results/rebuild/raw_manifest.csv (pdb_id, path, status, sha256, sha_matches_frozen).
Run from the repo root: PYTHONPATH=. python scripts/rebuild_fetch_raw.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import pandas as pd

from esr_analysis.download_structures import download_one
from esr_analysis.run_io import sha256
from rebuild_common import OUT, RAW_CACHE, RAW_LOCAL, ROOT, read_ids, rel


def main():
    frozen_sha = {}
    for line in (RAW_LOCAL / "sha256sums.txt").read_text().splitlines():
        if line.strip():
            digest, name = line.split()
            frozen_sha[name.replace(".cif", "").upper()] = digest
    RAW_CACHE.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb in read_ids():
        local = RAW_LOCAL / f"{pdb}.cif"
        if local.exists():
            path, status = local, "local"
        else:
            path = RAW_CACHE / f"{pdb}.cif"
            ok, status = download_one(pdb, path)
            if not ok:
                rows.append({"pdb_id": pdb, "path": "", "status": f"failed: {status}", "sha256": "", "sha_matches_frozen": False})
                continue
        digest = sha256(path)
        rows.append({"pdb_id": pdb, "path": rel(path), "status": status, "sha256": digest,
                     "sha_matches_frozen": digest == frozen_sha.get(pdb)})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(OUT / "raw_manifest.csv", index=False)
    print(manifest.status.value_counts().to_dict(), "sha matches:", int(manifest.sha_matches_frozen.sum()), "/", len(manifest))


if __name__ == "__main__":
    main()
