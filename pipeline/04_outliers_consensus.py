"""Stage 04 - consensus RMSD outliers: where they sit among the lid groups and why they are unusual.

Input: data/outliers/outlier_flags.csv (five detectors run in October on the RMSD matrix: k-nearest-neighbour
distance, distance to the group medoid, local outlier factor, isolation forest, Mahalanobis; consensus = flagged by
>= 3) and outlier_verdicts.csv (quality reason vs conformation), plus the stage-02 groups.
Test: hypergeometric (two-sided Fisher) for cl7, cl4 and the minor groups against all chains, and cl7 vs cl4 within the
main groups, at chain level and entry level ((entry, group) pairs); BH over the 8 tests.
Outputs (results/04_outliers_consensus/): enrichment_tests.csv, outliers.csv, numbers.json.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, hypergeom

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("04_outliers_consensus")
t = c.chains()
t["outlier"] = t.consensus_outlier.astype(str).eq("True")
rows = []
for level in ("chain", "entry"):
    u = t if level == "chain" else t.groupby(["entry_id", "group"]).outlier.max().reset_index()
    for test, univ, grp in (("cl7 vs all", None, "cl7"), ("cl4 vs all", None, "cl4"), ("minor vs all", None, "minor"),
                            ("cl7 vs cl4 (main only)", ("cl7", "cl4"), "cl7")):
        v = u if univ is None else u[u.group.isin(univ)]
        n, k, m = len(v), int(v.outlier.sum()), int((v.group == grp).sum())
        obs = int((v.outlier & (v.group == grp)).sum())
        tab = [[obs, m - obs], [k - obs, n - m - k + obs]]
        rows.append(dict(level=level, test=test, n_units=n, n_group=m, n_outliers=k, observed=obs,
                         expected=m * k / n, p_two_sided=fisher_exact(tab)[1],
                         p_depleted=hypergeom.cdf(obs, n, m, k), p_enriched=hypergeom.sf(obs - 1, n, m, k)))
tests = pd.DataFrame(rows)
tests["q_bh"] = c.bh(tests.p_two_sided)
tests.to_csv(OUT / "enrichment_tests.csv", index=False)
ol = t[t.outlier][["chain_key", "entry_id", "group", "n_methods", "outlier_verdict", "outlier_reasons", "resolution_A",
                   "peptide_bound", "mutation_status", "space_group", "ligand_id_primary", "action_class"]]
ol.to_csv(OUT / "outliers.csv", index=False)
tt = tests.set_index(["level", "test"])
c.write_numbers(OUT, {
    "n_consensus_outliers": int(t.outlier.sum()), "n_entries": int(ol.entry_id.nunique()),
    "by_group": ol.group.value_counts().to_dict(),
    "verdicts": ol.outlier_verdict.value_counts().to_dict(),
    "cl7_vs_cl4_chain": tt.loc[("chain", "cl7 vs cl4 (main only)"), ["observed", "expected", "q_bh"]].to_dict(),
    "cl7_vs_cl4_entry": tt.loc[("entry", "cl7 vs cl4 (main only)"), ["observed", "expected", "q_bh"]].to_dict(),
})
print(tests.round(4).to_string())
