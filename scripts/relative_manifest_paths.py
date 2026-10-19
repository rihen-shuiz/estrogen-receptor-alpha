"""Rewrite data/raw/download_manifest.csv with repo-relative paths (tracker task Q17).

The manifest stored absolute paths with the local user name. Each path becomes
data/raw/<PDB_ID>.cif, pointing at the file already in data/raw. No download runs.
Nothing is written unless every row names its own file and that file exists.

Run from anywhere: python scripts/relative_manifest_paths.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data/raw/download_manifest.csv"


def main():
    manifest = pd.read_csv(MANIFEST, dtype=str, keep_default_na=False)
    paths = []
    for pdb_id, old_path in zip(manifest.pdb_id, manifest.path):
        name = Path(old_path).name
        if name != f"{pdb_id}.cif":
            raise ValueError(f"{pdb_id}: manifest path {old_path!r} does not name {pdb_id}.cif")
        target = ROOT / "data/raw" / name
        if not target.is_file():
            raise FileNotFoundError(f"{pdb_id}: {target} missing; manifest not rewritten")
        paths.append(target.relative_to(ROOT).as_posix())
    manifest["path"] = paths
    manifest.to_csv(MANIFEST, index=False)
    print(f"rewrote {len(paths)} paths in {MANIFEST.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
