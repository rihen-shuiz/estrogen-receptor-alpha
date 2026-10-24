"""Rebuild tracker C05 (UMAP) and C06 (t-SNE) of the ligand fingerprint space.

Inputs: results/rebuild/fp_ecfp4.npy + fp_index.csv (C02 rebuild).
Method: C05 umap.UMAP(metric='jaccard', random_state=42, defaults n_neighbors=15, min_dist=0.1) on boolean
ECFP4 bits; rows UMAP leaves as disconnected vertices come back NaN and are kept with finite=False.
C06 sklearn TSNE(n_components=2, metric='jaccard', random_state=42) on boolean ECFP4 bits.
Both embeddings are stochastic/implementation dependent: compare by geometry, not coordinates.
Outputs: results/rebuild/umap_coords.csv (ligand_id, umap1, umap2, finite), results/rebuild/tsne_coords.csv
Run from the repo root: PYTHONPATH=. python scripts/rebuild_embeddings.py
"""
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
import pandas as pd

from rebuild_common import OUT, SEED, read_csv


def main():
    warnings.filterwarnings("ignore")
    import umap
    from sklearn.manifold import TSNE
    X = np.load(OUT / "fp_ecfp4.npy").astype(bool)
    ids = read_csv(OUT / "fp_index.csv").ligand_id.astype(str)
    u = umap.UMAP(metric="jaccard", random_state=SEED).fit_transform(X)
    finite = np.isfinite(u).all(1)
    pd.DataFrame({"ligand_id": ids, "umap1": u[:, 0], "umap2": u[:, 1], "finite": finite}).to_csv(OUT / "umap_coords.csv", index=False)
    print(f"C05 rebuild: {int(finite.sum())} of {len(ids)} UMAP rows finite")
    t = TSNE(n_components=2, metric="jaccard", random_state=SEED).fit_transform(X)
    pd.DataFrame({"ligand_id": ids, "tsne1": t[:, 0], "tsne2": t[:, 1], "finite": np.isfinite(t).all(1)}).to_csv(OUT / "tsne_coords.csv", index=False)
    print("C06 rebuild: t-SNE done")


if __name__ == "__main__":
    main()
