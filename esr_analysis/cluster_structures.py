"""Block 11: original average-linkage partitions only; no medoid/PCA task."""
from __future__ import annotations

from collections.abc import Sequence

import numpy as np

import pandas as pd

from scipy.cluster.hierarchy import fcluster, linkage

from scipy.spatial.distance import squareform

from .distance_matrix import DISTANCE_TOLERANCE_A, validate_distance_matrix

def structural_groups(ids, matrix, thresholds, primary):
    """Return cut assignments and average linkage for explicit indexed distances."""
    ids = list(ids)
    matrix = validate_distance_matrix(matrix, expected_size=len(ids))
    if not ids:
        raise ValueError("At least one chain is required for clustering")
    if len(set(ids)) != len(ids):
        raise ValueError("Chain keys must be unique")
    thresholds = [float(threshold) for threshold in thresholds]
    if not thresholds or not np.isfinite(thresholds).all() or any(t < 0 for t in thresholds):
        raise ValueError("Thresholds must be finite nonnegative distances")
    if primary not in thresholds:
        raise ValueError(f"Primary threshold {primary} not present")
    columns = [f"cluster_{threshold:.1f}A" for threshold in thresholds]
    if len(set(columns)) != len(columns):
        raise ValueError("Thresholds must have distinct column names at one decimal place")

    # SciPy linkage expects a condensed distance vector, not distance-row features.
    if len(ids) == 1:
        z = np.empty((0, 4), dtype=float)
    else:
        z = linkage(squareform(matrix, checks=False), method="average")
    assignments = pd.DataFrame({"chain_key": ids})
    for threshold, column in zip(thresholds, columns):
        labels = np.ones(1, dtype=int) if len(ids) == 1 else fcluster(z, t=threshold, criterion="distance")
        assignments[column] = labels

    return assignments, z


def main():
    from ._task_io import parser, prepare
    args = parser('Block 11: average-linkage report baseline only; no representatives or QC reclustering').parse_args()
    root, source, out, settings = prepare(args)
    ids = (source / 'rmsd_entry_index.txt').read_text().splitlines()
    matrix = np.load(source / 'rmsd_entry.npy')
    primary = settings['primary_structural_threshold_A']
    assignments, z = structural_groups(ids, matrix, [primary], primary)
    assignments.to_csv(out / 'entry_clusters.csv', index=False)
    np.save(out / 'linkage_entry.npy', z)
    print(f"Block 11: {assignments.iloc[:, 1].nunique()} baseline groups")


if __name__ == "__main__":
    main()
