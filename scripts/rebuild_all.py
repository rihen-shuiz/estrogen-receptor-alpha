"""Run every upstream rebuild script in dependency order, then the comparison.

Order: fetch raw mmCIFs (step 0) -> Q04 LBD extraction -> Q06 entries/QC -> X01 ESR2 entries ->
Q08/C01 ligands -> C02-C04 fingerprints -> C07/C10 K-means -> C05/C06 embeddings -> E03/E04 ->
X13 ChEMBL download (network, up to ~15 min; a failure is recorded as not_rebuilt, not fatal) ->
X14 ChEMBL-PDB matches -> comparison report.
All outputs go to results/rebuild/; frozen files are never written. Final report:
results/rebuild/rebuild_report.csv.
Run from the repo root: PYTHONPATH=. python scripts/rebuild_all.py [--skip-chembl]
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STEPS = [
    "rebuild_fetch_raw.py",
    "rebuild_lbd_extract.py",
    "rebuild_entries_qc.py",
    "rebuild_esr2_entries.py",
    "rebuild_ligands.py",
    "rebuild_fingerprints.py",
    "rebuild_kmeans.py",
    "rebuild_embeddings.py",
    "rebuild_esr1_esr2_maps.py",
    "rebuild_chembl_activities.py",
    "rebuild_chembl_pdb_matches.py",
    "rebuild_compare.py",
]
OPTIONAL = {"rebuild_chembl_activities.py"}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--skip-chembl", action="store_true", help="do not download ChEMBL (X14 then uses data/ copy)")
    args = parser.parse_args()
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    for step in STEPS:
        if args.skip_chembl and step == "rebuild_chembl_activities.py":
            continue
        print(f"== {step}", flush=True)
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / step)], cwd=ROOT, env=env)
        if result.returncode != 0:
            if step in OPTIONAL:
                print(f"   {step} failed (exit {result.returncode}); continuing, comparison records not_rebuilt")
                continue
            sys.exit(f"{step} failed with exit code {result.returncode}")


if __name__ == "__main__":
    main()
