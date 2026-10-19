"""QC-conditioned clustering delta of the September clustering (tracker task P03).

Recluster the September entries without those that fail combined QC (results/entries.csv, qc_pass
'fail'). Average linkage runs on the loaded September RMSD sub-matrix at the four P01 cuts. For each
cut and set, reports clusters, top-3 share, entropy in bits, n_90 (P02 method, seed 42), and the
full-set clusters that disappear because every member fails combined QC. Writes
results/provisional/qc_delta.csv. Failing entries are kept in P01's file, only excluded here.

Run from anywhere: python scripts/provisional_qc_delta.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provisional_common import (CUTS, PROVISIONAL, combined_qc, cut_labels, entropy_bits,  # noqa: E402
                                entry_id, load_september, saturation, top3_share)


def metrics(cut, set_name, labels):
    _, _, total, n_90 = saturation(labels)
    return {"cut": f"{cut:.1f}", "set": set_name, "n_entries": len(labels), "n_clusters": total,
            "top3_fraction": top3_share(labels), "entropy_bits": entropy_bits(labels), "n_90": n_90}


def main():
    chain_keys, rmsd = load_september()
    qc = combined_qc()
    passing = np.array([qc.get(entry_id(key)) == "pass" for key in chain_keys])
    clusters = pd.read_csv(PROVISIONAL / "clusters.csv", dtype=str, keep_default_na=False)
    keep = np.flatnonzero(passing)
    rows = []
    for cut in CUTS:
        full = clusters[f"cut_{cut:.1f}"].astype(int).to_numpy()
        rows.append(metrics(cut, "all_351", full) | {"n_clusters_absent": 0, "clusters_absent": ""})
        disappeared = [label for label in np.unique(full) if not passing[full == label].any()]
        absent = ";".join(f"{label}:{int((full == label).sum())}" for label in disappeared)
        pass_labels = cut_labels(rmsd[np.ix_(keep, keep)], cut)
        rows.append(metrics(cut, f"qc_pass_{len(keep)}", pass_labels)
                    | {"n_clusters_absent": len(disappeared), "clusters_absent": absent})
    table = pd.DataFrame(rows)
    table.to_csv(PROVISIONAL / "qc_delta.csv", index=False)
    print(table.to_string(index=False))
    print(f"wrote {PROVISIONAL / 'qc_delta.csv'}")


if __name__ == "__main__":
    main()
