"""Block 10: shared-position CA RMSD in the saved reference frame.

Input: ordered chain keys and aligned canonical mappings, minimum shared CA.
Method: no pair refit; preserve historical greedy missing-pair pruning.
Output: retained index, RMSD array (A), overlap counts, explicit drop reasons."""
from __future__ import annotations
from collections.abc import Mapping, Sequence
import numpy as np
from .align_structures import rmsd


def compute_distances(
    ids: Sequence[str],
    caches: Mapping[str, Mapping[int, np.ndarray]],
    min_common: int = 150,
) -> tuple[list[str], np.ndarray, np.ndarray, list[dict]]:
    """Compute historical distances and greedily remove missing-pair endpoints.

    Pruning drops the first chain with the most missing pairs, preserving input
    order for ties. This is greedy, not a guaranteed minimum removal set.
    Historical pruning stops at two chains: unresolved NaNs can remain and must
    be rejected before clustering. Each dropped record names incompatible
    remaining chains and their observed overlap counts at the removal step.
    """
    ids = list(ids)
    if len(set(ids)) != len(ids):
        raise ValueError("Chain keys must be unique")
    if min_common < 1:
        raise ValueError("min_common must be positive")
    n = len(ids)
    matrix = np.zeros((n, n), dtype=float)
    common_counts = np.zeros((n, n), dtype=int)
    for i in range(n):
        common_counts[i, i] = len(caches[ids[i]])
        for j in range(i + 1, n):
            common = sorted(set(caches[ids[i]]) & set(caches[ids[j]]))
            common_counts[i, j] = common_counts[j, i] = len(common)
            if len(common) < min_common:
                matrix[i, j] = matrix[j, i] = np.nan
            else:
                a = np.array([caches[ids[i]][p] for p in common])
                b = np.array([caches[ids[j]][p] for p in common])
                matrix[i, j] = matrix[j, i] = rmsd(a, b)

    keep = np.arange(n)
    work = matrix.copy()
    drop_reasons = {}
    while np.isnan(work).any() and len(keep) > 2:
        nan_counts = np.isnan(work).sum(axis=1)
        drop_local = int(np.argmax(nan_counts))
        drop_index = int(keep[drop_local])
        incompatible = keep[np.isnan(work[drop_local])]
        drop_reasons[ids[drop_index]] = {
            "chain_key": ids[drop_index],
            "reason": "insufficient_common_ca",
            "pruning_step": len(drop_reasons) + 1,
            "min_common": min_common,
            "incompatible_chain_keys": [ids[j] for j in incompatible],
            "common_ca_counts": [int(common_counts[drop_index, j]) for j in incompatible],
        }
        keep = np.delete(keep, drop_local)
        work = matrix[np.ix_(keep, keep)]

    kept_ids = [ids[i] for i in keep]
    dropped = [drop_reasons[key] for key in ids if key in drop_reasons]
    kept_counts = common_counts[np.ix_(keep, keep)]
    return kept_ids, work, kept_counts, dropped



def main():
    from ._task_io import parser, prepare, aligned_caches, json_write
    args = parser('Block 10: ordered shared-position RMSD from bundled saved aligned observations').parse_args()
    root, source, out, settings = prepare(args)
    index = (source / 'rmsd_entry_index.txt').read_text().splitlines()
    caches = aligned_caches(source / 'aligned_coordinates.npz', index)
    ids, matrix, overlap, excluded = compute_distances(index, caches, settings['min_pairwise_common_ca'])
    np.save(out / 'rmsd_entry.npy', matrix)
    np.save(out / 'overlap_entry.npy', overlap)
    (out / 'rmsd_entry_index.txt').write_text('\n'.join(ids) + '\n')
    json_write(out / 'distance_exclusions.json', excluded)
    print(f'Block 10: {len(ids)} entries, {len(ids)*(len(ids)-1)//2} distinct pairs')


if __name__ == "__main__":
    main()
