"""Tracker Q14: do conformational singletons / rare clusters have worse resolution or R-free?

Inputs: results/entry_clusters.csv (S08, all 788 eligible chains), results/qc_tiers.csv (K05: resolution, R-free,
mean B and qc_label per entry; read, not recomputed).
Rare cluster: fewer than RARE_BELOW chains at that cut (a singleton is the extreme case). An entry counts as rare
when any of its chains sits in a rare cluster, otherwise common. The test is per entry, because resolution and
R-free belong to the entry and repeating them per chain would inflate n.
Test: two-sided Mann-Whitney U, rare vs common entries, for resolution and R-free, at each cut. The rare groups are
small, so the p-values are descriptive; nothing is filtered or dropped.
Outputs: results/qc_outliers.csv (one row per chain: cluster, cluster size, rare flag at each cut, entry QC values),
results/qc_outliers_stats.csv (one row per cut and metric), docs/qc_outliers.md, figures/qc_outliers.png.
Run from the repo root: python scripts/qc_outliers.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUTS = ["1.0", "1.5", "2.0", "2.5"]
PRIMARY = f'{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
RARE_BELOW = 10
METRICS = {"resolution_A": "Resolution (Å)", "r_free": "R-free"}
INK, INK_SECONDARY, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
COMMON, RARE = "#2a78d6", "#eb6834"


def main():
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    qc = pd.read_csv(ROOT / "results/qc_tiers.csv", usecols=["entry_id", "resolution_A", "r_free", "mean_b_A2", "qc_label"])
    chains = clusters.drop(columns="qc_label").merge(qc, on="entry_id", how="left", validate="many_to_one")
    for cut in CUTS:
        size = chains[f"cut_{cut}"].map(chains[f"cut_{cut}"].value_counts())
        chains[f"size_{cut}"] = size
        chains[f"rare_{cut}"] = size < RARE_BELOW
    chains.to_csv(ROOT / "results/qc_outliers.csv", index=False)

    entries = chains.groupby("entry_id").agg(
        resolution_A=("resolution_A", "first"), r_free=("r_free", "first"),
        **{f"rare_{cut}": (f"rare_{cut}", "any") for cut in CUTS})
    stats = []
    for cut in CUTS:
        rare = entries[entries[f"rare_{cut}"]]
        common = entries[~entries[f"rare_{cut}"]]
        for metric in METRICS:
            a, b = rare[metric].dropna(), common[metric].dropna()
            p = float(mannwhitneyu(a, b, alternative="two-sided").pvalue) if len(a) and len(b) else np.nan
            stats.append({"cut": cut, "metric": metric, "rare_below_chains": RARE_BELOW,
                          "n_rare_entries": len(a), "n_common_entries": len(b),
                          "median_rare": a.median() if len(a) else np.nan, "median_common": b.median(),
                          "mannwhitney_p": p})
    stats = pd.DataFrame(stats)
    stats.to_csv(ROOT / "results/qc_outliers_stats.csv", index=False)

    rare_chains = chains[chains[f"rare_{PRIMARY}"]].sort_values([f"cut_{PRIMARY}", "chain_key"])
    figure(entries, rare_chains)
    write_doc(stats, rare_chains, entries)
    print(stats.to_string(index=False))


def figure(entries, rare_chains):
    fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=150, facecolor=SURFACE)
    rng = np.random.default_rng(0)
    groups = [("common", ~entries[f"rare_{PRIMARY}"], COMMON), ("rare", entries[f"rare_{PRIMARY}"], RARE)]
    for ax, (metric, label) in zip(axes, METRICS.items()):
        ax.set_facecolor(SURFACE)
        for x, (name, mask, colour) in enumerate(groups):
            values = entries.loc[mask, metric].dropna()
            ax.scatter(x + rng.uniform(-0.18, 0.18, len(values)), values, s=10, color=colour,
                       alpha=0.5 if name == "common" else 0.9, lw=0)
            ax.hlines(values.median(), x - 0.28, x + 0.28, color=INK, lw=1.5)
        ax.set_xticks([0, 1], [f"common\n(n={int(groups[0][1].sum())})", f"rare\n(n={int(groups[1][1].sum())})"])
        ax.set_ylabel(label, color=INK_SECONDARY, fontsize=9)
        ax.tick_params(colors=MUTED, labelsize=8)
        ax.grid(axis="y", color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
    fig.suptitle(f"Entries in rare clusters (< {RARE_BELOW} chains at {PRIMARY} Å) vs the rest; line = median",
                 x=0.02, ha="left", fontsize=10, color=INK, fontweight="bold")
    fig.tight_layout()
    fig.savefig(ROOT / "figures/qc_outliers.png", facecolor=SURFACE)
    plt.close(fig)


def write_doc(stats, rare_chains, entries):
    primary = stats[stats.cut == PRIMARY].set_index("metric")
    res, rfree = primary.loc["resolution_A"], primary.loc["r_free"]
    lines = [
        "# QC question 1: do conformational singletons have worse resolution?",
        "",
        "Generated by `scripts/qc_outliers.py` from `results/entry_clusters.csv` (S08, 788 chains) and "
        "`results/qc_tiers.csv` (K05). Tables: `results/qc_outliers.csv`, `results/qc_outliers_stats.csv`; "
        "plot: `figures/qc_outliers.png`.",
        "",
        f"A cluster is rare when it holds fewer than {RARE_BELOW} chains. An entry is rare when any of its chains "
        "is in a rare cluster. Resolution and R-free are entry properties, so the comparison is per entry "
        f"({len(entries)} entries). Two-sided Mann-Whitney U.",
        "",
        f"**Stats line ({PRIMARY} Å cut):** {int(res.n_rare_entries)} rare vs {int(res.n_common_entries)} common entries; "
        f"median resolution {res.median_rare:.2f} vs {res.median_common:.2f} Å (p = {res.mannwhitney_p:.3g}); "
        f"median R-free {rfree.median_rare:.3f} vs {rfree.median_common:.3f} (p = {rfree.mannwhitney_p:.3g}).",
        "",
        "| cut (Å) | metric | rare entries | common entries | median rare | median common | p |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in stats.itertuples():
        lines.append(f"| {row.cut} | {row.metric} | {row.n_rare_entries} | {row.n_common_entries} | "
                     f"{row.median_rare:.3f} | {row.median_common:.3f} | {row.mannwhitney_p:.3g} |")
    lines += ["", f"Chains in rare clusters at {PRIMARY} Å:", "",
              "| chain | cluster | cluster size | resolution (Å) | R-free | qc_label |", "|---|---|---|---|---|---|"]
    for row in rare_chains.to_dict("records"):
        lines.append(f"| {row['chain_key']} | {row['cut_' + PRIMARY]} | {row['size_' + PRIMARY]} | "
                     f"{row['resolution_A']:.2f} | {row['r_free']:.3f} | {row['qc_label']} |")
    lines += ["", "The rare groups are small, so these p-values describe this dataset and are not corrected for the "
              "eight tests in the table. A rare cluster can reflect a real conformational state, a crystal-packing "
              "effect or poor data; resolution and R-free alone cannot separate these.", ""]
    (ROOT / "docs/qc_outliers.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
