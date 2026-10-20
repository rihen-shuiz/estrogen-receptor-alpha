"""Tracker Q15: do high-B entries cluster together or scatter across the structural clusters?

Inputs: results/entry_clusters.csv (S08, 788 chains, cut at primary_structural_threshold_A),
results/helix_bfactor_by_chain.csv (K03: b_group = tertile of the chain's mean C-alpha B),
results/qc_tiers.csv (K05: entry mean B and tier_b). All read, not recomputed; nothing is dropped.
Chain level: cross-table of b_group (low / mid / high) against cluster, with the share of high-B chains per cluster
and the enrichment (share in cluster / share overall). Chi-square on the clusters with at least MIN_TEST chains
(small clusters are listed but would break the expected-count rule).
Entry level (chains of one entry share crystal and refinement, so chain counts overstate n): entries whose chains all
sit in one tested cluster, tier_b (A = mean B within the tier-A limit, B = above) against that cluster, Fisher exact
when two clusters are tested, chi-square otherwise.
Outputs: results/qc_bfactor_clusters.csv (one row per cluster: counts per b_group, high share, enrichment, entry
tier_b counts) and results/qc_bfactor_clusters_tests.json.
Run from the repo root: python scripts/qc_bfactor_clusters.py
"""
import json
from pathlib import Path

import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

ROOT = Path(__file__).resolve().parents[1]
CUT = f'cut_{json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())["primary_structural_threshold_A"]:.1f}'
MIN_TEST = 20
GROUPS = ["low", "mid", "high"]


def main():
    clusters = pd.read_csv(ROOT / "results/entry_clusters.csv")
    bfac = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv", usecols=["chain_key", "mean_b_chain", "b_group"])
    qc = pd.read_csv(ROOT / "results/qc_tiers.csv", usecols=["entry_id", "mean_b_A2", "tier_b"])
    chains = clusters.merge(bfac, on="chain_key", how="left", validate="one_to_one")
    if chains.b_group.isna().any():
        raise ValueError("chains without a K03 b_group")

    table = pd.crosstab(chains[CUT], chains.b_group).reindex(columns=GROUPS, fill_value=0)
    table.columns = [f"n_{g}" for g in GROUPS]
    table["n_chains"] = table.sum(axis=1)
    table["n_entries"] = chains.groupby(CUT).entry_id.nunique()
    table["median_mean_b_chain"] = chains.groupby(CUT).mean_b_chain.median()
    overall_high = (chains.b_group == "high").mean()
    table["high_share"] = table.n_high / table.n_chains
    table["high_enrichment"] = table.high_share / overall_high
    table["tested"] = table.n_chains >= MIN_TEST

    tested = table[table.tested].index.tolist()
    chi2, p_chain, dof, _ = chi2_contingency(table.loc[tested, [f"n_{g}" for g in GROUPS]])

    per_entry = chains.groupby("entry_id")[CUT].agg(["nunique", "first"])
    single = per_entry[(per_entry["nunique"] == 1) & per_entry["first"].isin(tested)].rename(columns={"first": "cluster"})
    single = single.join(qc.set_index("entry_id"))
    entry_table = pd.crosstab(single.cluster, single.tier_b)
    for tier in entry_table.columns:
        table[f"entries_tier_b_{tier}"] = entry_table[tier]
    if len(tested) == 2 and entry_table.shape[1] == 2:
        odds, p_entry = fisher_exact(entry_table.to_numpy())
        entry_test = {"test": "fisher_exact", "odds_ratio": float(odds), "p": float(p_entry)}
    else:
        stat, p_entry, _, _ = chi2_contingency(entry_table)
        entry_test = {"test": "chi2", "statistic": float(stat), "p": float(p_entry)}

    table.index.name = "cluster"
    table.reset_index().assign(cut=CUT.removeprefix("cut_")).to_csv(ROOT / "results/qc_bfactor_clusters.csv", index=False)
    tests = {
        "task": "Q15", "cut": CUT.removeprefix("cut_"), "n_chains": len(chains), "overall_high_share": float(overall_high),
        "chain_level": {"clusters_tested": [int(c) for c in tested], "min_chains_to_test": MIN_TEST,
                        "chi2": float(chi2), "dof": int(dof), "p": float(p_chain)},
        "entry_level": {"n_entries": len(single), "entries_split_or_untested": int(len(per_entry) - len(single)),
                        "tier_b_by_cluster": {str(k): v for k, v in entry_table.to_dict(orient="index").items()},
                        **entry_test},
        "note": "chain-level counts repeat entries with several chains; the entry-level test is the safer reading",
    }
    (ROOT / "results/qc_bfactor_clusters_tests.json").write_text(json.dumps(tests, indent=2) + "\n")
    print(table.to_string())
    print(json.dumps({k: tests[k] for k in ("chain_level", "entry_level")}, indent=1))


if __name__ == "__main__":
    main()
