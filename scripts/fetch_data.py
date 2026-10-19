"""Rebuild git-ignored inputs from the tracked ID lists.

Steps, in order:
1. Download every ESR1 and ESR2 mmCIF listed in data/pdb_ids_*.txt into data/raw.
   Valid files already present are kept.
2. Download CCD ligand files for every non-water ligand in results/entries.csv that
   has no local copy in esr_analysis/resources/data/ccd. They go to data/ccd.
3. Rebuild data/lbd with esr_analysis.extract_lbd.

Run from the repo root:
    python scripts/fetch_data.py          # download and rebuild
    python scripts/fetch_data.py --plan   # print what would be fetched, no network
"""
import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pandas as pd
import requests

from esr_analysis.download_structures import download_structures

CCD_URL = "https://files.rcsb.org/ligands/view/{code}.cif"
LOCAL_CCD = ROOT / "esr_analysis/resources/data/ccd"
CCD_DIR = ROOT / "data/ccd"


def read_ids():
    ids = []
    for name in ["data/pdb_ids_esr1.txt", "data/pdb_ids_esr2.txt"]:
        ids += [line.strip() for line in (ROOT / name).read_text().splitlines() if line.strip()]
    return ids


def ligand_codes():
    entries = pd.read_csv(ROOT / "results/entries.csv")
    codes = set()
    for value in entries.ligand_ids.fillna(""):
        codes.update(code for code in value.split(";") if code)
    return sorted(codes)


def missing_ccd(codes):
    return [code for code in codes
            if not (LOCAL_CCD / f"{code}.cif").exists() and not (CCD_DIR / f"{code}.cif").exists()]


def fetch_ccd(codes):
    CCD_DIR.mkdir(parents=True, exist_ok=True)
    for code in codes:
        response = requests.get(CCD_URL.format(code=code), timeout=60)
        response.raise_for_status()
        (CCD_DIR / f"{code}.cif").write_bytes(response.content)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--plan", action="store_true", help="print the steps and counts without network calls")
    args = parser.parse_args()

    ids = read_ids()
    codes = ligand_codes()
    todo_ccd = missing_ccd(codes)
    print(f"structures: {len(ids)} IDs -> data/raw")
    print(f"ligands: {len(codes)} codes, {len(todo_ccd)} CCD files missing -> data/ccd")
    print("rebuild: data/lbd via esr_analysis.extract_lbd")
    if args.plan:
        return

    manifest = download_structures(ids, ROOT / "data/raw", root=ROOT)
    if not manifest.download_ok.all():
        raise RuntimeError("Some structure downloads failed; see data/raw/download_manifest.csv")

    fetch_ccd(todo_ccd)

    cmd = [sys.executable, "-m", "esr_analysis.extract_lbd", "--ids", *ids,
           "--raw-folder", str(ROOT / "data/raw"), "--output", str(ROOT / "data/lbd")]
    subprocess.run(cmd, cwd=ROOT, check=True)
    print("done")


if __name__ == "__main__":
    main()
