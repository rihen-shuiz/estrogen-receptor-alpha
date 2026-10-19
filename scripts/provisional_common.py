"""Shared helpers for the provisional September analyses (tracker tasks P01-P04).

Every input is read, not recomputed: the September RMSD matrix and its entry index
(results/task_outputs/block_10) and the combined QC flag per entry (results/entries.csv).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[1]
RMSD_PATH = ROOT / "results/task_outputs/block_10/rmsd_entry.npy"
INDEX_PATH = ROOT / "results/task_outputs/block_10/rmsd_entry_index.txt"
PROVISIONAL = ROOT / "results/provisional"
CUTS = (1.0, 1.5, 2.0, 2.5)
SEED = 42          # analysis_settings.json, random_seed
N_ORDERINGS = 100


def load_september():
    """Chain keys (one per entry, matrix order) and the 351 x 351 September RMSD matrix in Å."""
    chain_keys = [line.strip() for line in INDEX_PATH.read_text().splitlines() if line.strip()]
    return chain_keys, np.load(RMSD_PATH)


def entry_id(chain_key):
    return chain_key.split("|")[0]


def combined_qc():
    entries = pd.read_csv(ROOT / "results/entries.csv", dtype=str, keep_default_na=False)
    return dict(zip(entries.entry_id, entries.qc_pass))


def cut_labels(rmsd, cut):
    """Average-linkage clusters at a distance cut in Å; labels as scipy numbers them (from 1)."""
    tree = linkage(squareform(rmsd, checks=False), method="average")
    return fcluster(tree, t=cut, criterion="distance")


def cluster_sizes(labels):
    """Cluster sizes, largest first."""
    return pd.Series(labels).value_counts()


def entropy_bits(labels):
    proportions = cluster_sizes(labels).to_numpy() / len(labels)
    return float(-(proportions * np.log2(proportions)).sum())


def top3_share(labels):
    return float(cluster_sizes(labels).to_numpy()[:3].sum() / len(labels))


def orderings(n):
    """N_ORDERINGS random orderings of n entries, from one seeded generator."""
    rng = np.random.default_rng(SEED)
    return np.array([rng.permutation(n) for _ in range(N_ORDERINGS)])


def clusters_seen(labels, orders):
    """Distinct clusters seen after the first N entries, one row per ordering."""
    codes = pd.factorize(labels)[0]
    seen_curve = np.zeros((len(orders), len(labels)))
    for row, order in enumerate(orders):
        seen = np.zeros(codes.max() + 1, dtype=bool)
        count = 0
        for position, entry in enumerate(order):
            if not seen[codes[entry]]:
                seen[codes[entry]] = True
                count += 1
            seen_curve[row, position] = count
    return seen_curve


def saturation(labels):
    """(mean, sd, total clusters, n_90) over the orderings; n_90 is the smallest N whose mean reaches 90% of the total."""
    curves = clusters_seen(labels, orderings(len(labels)))
    mean, sd = curves.mean(axis=0), curves.std(axis=0, ddof=1)
    total = len(np.unique(labels))
    n_90 = int(np.argmax(mean >= 0.9 * total)) + 1
    return mean, sd, total, n_90
