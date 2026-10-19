"""Provisional threshold sensitivity of the September clustering (tracker task P01).

Average linkage on the September RMSD matrix (loaded, not recomputed), cut at 1.0, 1.5, 2.0 and
2.5 Å. Writes results/provisional/clusters.csv (one label column per cut),
results/provisional/entropy_sensitivity.csv (clusters, top-3 share, entropy in bits per cut) and
figures/provisional_threshold_sensitivity.png. The 1.5 Å cut is checked against the September
groups (212, 124, 6, 2, 2, 2, 1, 1, 1) and partition in results/task_outputs/block_11/entry_clusters.csv.

Run from anywhere: python scripts/provisional_clusters.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import adjusted_rand_score  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from provisional_common import (CUTS, PROVISIONAL, ROOT, cluster_sizes, cut_labels,  # noqa: E402
                                combined_qc, entropy_bits, entry_id, load_september, top3_share)

SEP_SIZES = [212, 124, 6, 2, 2, 2, 1, 1, 1]
SEP_CUT = 1.5
FIGURE = ROOT / "figures/provisional_threshold_sensitivity.png"
INK, MUTED, GRID, AXIS, SURFACE, SERIES = "#0b0b0b", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb", "#2a78d6"


def sensitivity_table(labels):
    rows = []
    for cut in CUTS:
        sizes = cluster_sizes(labels[cut])
        rows.append({"cut": f"{cut:.1f}", "n_clusters": len(sizes), "top3_fraction": top3_share(labels[cut]),
                     "entropy_bits": entropy_bits(labels[cut]), "largest_cluster": int(sizes.iloc[0]),
                     "n_entries": len(labels[cut])})
    return pd.DataFrame(rows)


def reproduction_check(chain_keys, labels):
    sep = pd.read_csv(ROOT / "results/task_outputs/block_11/entry_clusters.csv", dtype=str, keep_default_na=False)
    sep_label = dict(zip(sep.chain_key, sep[f"cluster_{SEP_CUT}A"].astype(int)))
    sizes = sorted(cluster_sizes(labels[SEP_CUT]).tolist(), reverse=True)
    agreement = adjusted_rand_score([sep_label[key] for key in chain_keys], labels[SEP_CUT])
    matches = sizes == SEP_SIZES and agreement == 1.0
    print(f"1.5 Å sizes {sizes}; adjusted Rand index vs September partition {agreement:.6f}")
    print("REPRODUCES September 1.5 Å groups" if matches else "PART: does not reproduce September 1.5 Å groups")
    return matches


def figure(table):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), dpi=150, facecolor=SURFACE)
    panels = [("n_clusters", "Clusters (count)"), ("top3_fraction", "Top-3 share of entries"),
              ("entropy_bits", "Entropy (bits)")]
    for ax, (column, title) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        cuts = [float(c) for c in table.cut]
        ax.plot(cuts, table[column], color=SERIES, lw=1.2, marker="o", ms=6, mec=SURFACE, mew=1)
        ax.axvline(SEP_CUT, color=AXIS, lw=1)
        point = table.loc[table.cut == f"{SEP_CUT:.1f}", column].iloc[0]
        ax.annotate(f"{point:.3g}", (SEP_CUT, point), textcoords="offset points", xytext=(7, 6),
                    color=INK, fontsize=9)
        ax.set_title(title, loc="left", fontsize=10, color=INK)
        ax.set_xlabel("Distance cut (Å)", color=MUTED, fontsize=9)
        ax.set_xticks(cuts)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(colors=MUTED, labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
    fig.tight_layout()
    fig.savefig(FIGURE, facecolor=SURFACE)
    plt.close(fig)


def main():
    chain_keys, rmsd = load_september()
    labels = {cut: cut_labels(rmsd, cut) for cut in CUTS}
    clusters = pd.DataFrame({"entry_id": [entry_id(key) for key in chain_keys], "chain_key": chain_keys})
    clusters["qc_combined"] = clusters.entry_id.map(combined_qc())
    for cut in CUTS:
        clusters[f"cut_{cut:.1f}"] = labels[cut]
    PROVISIONAL.mkdir(parents=True, exist_ok=True)
    clusters.to_csv(PROVISIONAL / "clusters.csv", index=False)
    table = sensitivity_table(labels)
    table.to_csv(PROVISIONAL / "entropy_sensitivity.csv", index=False)
    figure(table)
    print(table.to_string(index=False))
    reproduction_check(chain_keys, labels)
    print(f"wrote {PROVISIONAL / 'clusters.csv'}, {PROVISIONAL / 'entropy_sensitivity.csv'} and {FIGURE}")


if __name__ == "__main__":
    main()
