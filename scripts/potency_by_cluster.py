"""Tracker X15 / X16: ChEMBL potency of PDB ligands grouped by structural cluster, and the tests.

Inputs (read only): results/chain_ligand_of_interest.csv (chain -> ligand of interest), results/chembl_pdb_matches.csv
(X14, PDB ligand -> ChEMBL molecule, Tanimoto >= 0.9 or exact), data/chembl_esr1_activities.csv (ESR1 activities),
results/entry_clusters.csv (S08 labels at the primary cut), results/h12_metric.csv (S16).
Potency: pChEMBL (-log10 molar) of IC50, Ki and Kd rows with standard_relation "=", median per ChEMBL molecule;
a PDB ligand matched to several molecules takes the median over them.
Unit = one PDB ligand (chains of one ligand are not independent). Its cluster is the cluster holding most of its
chains (ties: the larger cluster, then the lower number); n_chains and the share in that cluster are kept.
Tests (X16, results/potency_tests.csv): Kruskal-Wallis across clusters with >= MIN_GROUP ligands; Mann-Whitney between
the two largest clusters; Spearman between potency and the ligand's median H12 RMSD (continuous structural reading).
Figure (X15, figures/potency_by_cluster.png): potency as nM on a log axis per cluster, median marked.
Also writes results/potency_by_cluster.csv (one row per ligand).
Run from the repo root: python scripts/potency_by_cluster.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import kruskal, mannwhitneyu, spearmanr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CUT = f'cut_{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
TYPES = ["IC50", "Ki", "Kd"]
MIN_GROUP = 3
SURFACE, INK, INK_SECONDARY, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SERIES = "#2a78d6"


def ligand_table():
    activities = pd.read_csv(ROOT / "data/chembl_esr1_activities.csv", low_memory=False)
    activities = activities[activities.standard_type.isin(TYPES) & activities.standard_relation.eq("=")
                            & activities.pchembl_value.notna()]
    per_molecule = activities.groupby("molecule_chembl_id").agg(pchembl=("pchembl_value", "median"),
                                                                n_activities=("activity_id", "size"))
    matches = pd.read_csv(ROOT / "results/chembl_pdb_matches.csv", dtype={"ligand_id": str})
    matches = matches.join(per_molecule, on="molecule_chembl_id").dropna(subset=["pchembl"])
    per_ligand = matches.groupby("ligand_id").agg(pchembl=("pchembl", "median"), n_molecules=("molecule_chembl_id", "nunique"),
                                                  n_activities=("n_activities", "sum"), min_tanimoto=("tanimoto", "min"))

    loi = pd.read_csv(ROOT / "results/chain_ligand_of_interest.csv", dtype={"ligand_id": str})
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")[["chain_key", CUT]]
    h12 = pd.read_csv(ROOT / "results/h12_metric.csv", usecols=["chain_key", "h12_rmsd"])
    chains = loi.merge(clusters, on="chain_key").merge(h12, on="chain_key", how="left")
    chains = chains[chains.ligand_id.isin(per_ligand.index)]
    size = clusters[CUT].value_counts()
    rows = []
    for ligand, group in chains.groupby("ligand_id"):
        counts = group[CUT].value_counts()
        best = sorted(counts.index, key=lambda c: (-counts[c], -size[c], c))[0]
        rows.append({"ligand_id": ligand, "cluster": int(best), "n_chains": len(group),
                     "n_entries": group.entry_id.nunique(), "share_in_cluster": counts[best] / len(group),
                     "median_h12_rmsd_A": group.h12_rmsd.median()})
    table = pd.DataFrame(rows).join(per_ligand, on="ligand_id")
    table["potency_nM"] = 10 ** (9 - table.pchembl)
    return table.sort_values(["cluster", "ligand_id"]).reset_index(drop=True)


def tests(table):
    counts = table.cluster.value_counts()
    groups = [g.pchembl.to_numpy() for c, g in table.groupby("cluster") if counts[c] >= MIN_GROUP]
    rows = []
    stat, p = kruskal(*groups)
    rows.append({"test": "kruskal_wallis", "groups": ",".join(str(c) for c in sorted(counts[counts >= MIN_GROUP].index)),
                 "n": int(sum(len(g) for g in groups)), "statistic": stat, "p_value": p})
    a, b = counts.index[:2]
    x, y = table.pchembl[table.cluster == a], table.pchembl[table.cluster == b]
    stat, p = mannwhitneyu(x, y, alternative="two-sided")
    rows.append({"test": "mann_whitney", "groups": f"{a} vs {b}", "n": len(x) + len(y), "statistic": stat, "p_value": p,
                 "median_a": x.median(), "median_b": y.median()})
    sub = table.dropna(subset=["median_h12_rmsd_A"])
    rho, p = spearmanr(sub.pchembl, sub.median_h12_rmsd_A)
    rows.append({"test": "spearman_pchembl_vs_h12_rmsd", "groups": "all", "n": len(sub), "statistic": rho, "p_value": p})
    out = pd.DataFrame(rows)
    out["cut"] = CUT.removeprefix("cut_")
    out["unit"] = "PDB ligand (median pChEMBL over IC50/Ki/Kd, relation '=')"
    return out


def figure(table):
    clusters = sorted(table.cluster.unique(), key=lambda c: -(table.cluster == c).sum())
    fig, ax = plt.subplots(figsize=(7.5, 4.5), dpi=150, facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    rng = np.random.default_rng(0)
    for x, cluster in enumerate(clusters):
        values = table.potency_nM[table.cluster == cluster]
        ax.scatter(x + rng.uniform(-0.2, 0.2, len(values)), values, s=16, color=SERIES, alpha=0.6, lw=0)
        ax.hlines(values.median(), x - 0.3, x + 0.3, color=INK, lw=1.5)
    ax.set_yscale("log")
    ticks = [10.0 ** k for k in range(-2, 6) if table.potency_nM.min() / 3 <= 10.0 ** k <= table.potency_nM.max() * 3]
    ax.set_yticks(ticks, [f"{t:g}" for t in ticks])
    ax.minorticks_off()
    ax.invert_yaxis()
    ax.set_xticks(range(len(clusters)), [f"cluster {c}\n(n={(table.cluster == c).sum()})" for c in clusters])
    ax.set_ylabel("Potency, nM (IC50/Ki/Kd median; log scale, stronger up)", color=INK_SECONDARY, fontsize=9)
    ax.set_title(f"ChEMBL potency of PDB ligands by structural cluster ({CUT[4:]} Å); line = median",
                 loc="left", fontsize=10, color=INK, fontweight="bold")
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    fig.tight_layout()
    fig.savefig(ROOT / "figures/potency_by_cluster.png", facecolor=SURFACE)
    plt.close(fig)


def main():
    table = ligand_table()
    table.to_csv(ROOT / "results/potency_by_cluster.csv", index=False)
    figure(table)
    result = tests(table)
    result.to_csv(ROOT / "results/potency_tests.csv", index=False)
    print(f"{len(table)} ligands with potency; per cluster {table.cluster.value_counts().to_dict()}")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
