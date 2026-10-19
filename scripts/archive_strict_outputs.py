"""Tracker K00: copy the strict (351-chain, v2) outputs into results/archive_v2_strict/ before the rewrite.

Copies only: the sources stay where they are. Top-level files keep their file names; the alignment
folder goes into archive/alignment/. MANIFEST.csv lists every copy with its sha256 and source path.
The script refuses to run when the archive already exists, so an archive is never overwritten.
Run from the repo root: python scripts/archive_strict_outputs.py
"""
import hashlib
import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/archive_v2_strict"
TOP_LEVEL = [
    "results/rmsd_entry.npy",
    "results/rmsd_entry_index.csv",
    "results/entry_clusters.csv",
    "results/linkage_entry.npy",
    "results/cluster_stats.csv",
    "results/medoids.csv",
    "results/qc_clustering_delta.csv",
    "results/h12_metric.csv",
    "results/alignment_core_rmsd.csv",
    "figures/dendrogram.png",
    "figures/structural_sensitivity.png",
    "docs/limitations.md",
]
ALIGNMENT = [
    "results/alignment/aligned_coordinates.npz",
    "results/alignment/alignment_transforms.json",
    "results/alignment/global_alignment.csv",
    "results/alignment/chain_keys.txt",
    "results/alignment/rmsd_entry_index.txt",
    "results/alignment/overlap_entry.npy",
    "results/alignment/distance_exclusions.json",
]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if ARCHIVE.exists():
        raise FileExistsError(f"{ARCHIVE} already exists; refusing to overwrite an archive")
    sources = [ROOT / name for name in TOP_LEVEL + ALIGNMENT]
    missing = [str(path) for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Cannot archive, missing sources: {missing}")

    rows = []
    for name in TOP_LEVEL + ALIGNMENT:
        source = ROOT / name
        relative = Path(Path(name).name) if name in TOP_LEVEL else Path("alignment") / Path(name).name
        target = ARCHIVE / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if sha256(source) != sha256(target):
            raise ValueError(f"Copy differs from source: {name}")
        rows.append({"file": relative.as_posix(), "sha256": sha256(target),
                     "bytes": target.stat().st_size, "source": name})
    manifest = pd.DataFrame(rows)
    manifest.to_csv(ARCHIVE / "MANIFEST.csv", index=False)
    print(f"K00: {len(manifest)} files copied to {ARCHIVE.relative_to(ROOT)}; sources untouched")


if __name__ == "__main__":
    main()
