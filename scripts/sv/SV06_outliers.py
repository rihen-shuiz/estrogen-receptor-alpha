"""SV06: the Y16 consensus outliers on the clustering map (RUN_SV1028). No new outlier method.

Inputs (found through tracker/esr1.tasks card outs and the DONE lines of tracker/esr1.log, rule 1.10):
  Y15/Y16  results/ext/outlier_flags.csv (Y15 out: 22 chains with n_methods >= 3 = consensus_outlier) and
           results/ext/outlier_verdicts.csv (listed as evidence of Y16 DONE in tracker/esr1.log: verdict
           quality / conformation and the reason text, one row per consensus outlier).
  Y03      results/ext/embedding_rmsd.csv, columns mds1/mds2 (the shared map = metric MDS, Y03 DONE note).
  Y07      results/ext/leaf_order.csv (dendrogram leaf rank of every chain).
  S08/S12  results/medoids.csv (medoid chain of every S08 cluster at the 1.5 A cut).
  Y10      results/ext/cluster_displacement_profile.csv, columns medoid_cl7/medoid_cl4 (complete-H12 medoids of the
           two main clusters; the same two medoids the Y15 medoid detector used).
  S04/S05  results/rmsd_entry.npy + results/rmsd_entry_index.csv (788 x 788 C-alpha RMSD, loaded, not recomputed).
  SV00     results/sv/SV00/chain_table.csv (s08_cluster, auth_asym_id; same row order as the matrix).

Rules fixed before any result was looked at:
  Cluster groups: cl7, cl4 (the two main S08 clusters) and 'minor' = S08 clusters 1, 2, 3, 5, 6 pooled (rule 1.8).
  Labels: 'PDB|author chain' (auth_asym_id); chain_key keeps the label_asym_id, the CSV carries both.
  Panel B distances: own medoid = medoid of the chain's main cluster; for a minor-cluster chain the main cluster
    is the one whose medoid is nearer (same rule as Y15 'home_cluster'). Other medoid = the medoid of the other
    main cluster. The RMSD to the medoid of the chain's own S08 cluster (results/medoids.csv) is also written
    (rmsd_own_s08_medoid). own_other_ratio = rmsd_own / rmsd_other (1 = on the diagonal, equidistant from the two
    main medoids). The run file sets no threshold for 'near the diagonal', so none is applied.
    SUBSTITUTION made after the first run (recorded in SUMMARY line 1): the rule was first applied with the S12
    medoids of results/medoids.csv. Their cl4 medoid 2QE4|B|1 has no modelled H12 (h12_label not_modelled), and
    the S04 RMSD is taken over the C-alpha positions two chains share, so 2QE4|B|1 sits 0.36 A from the cl7 medoid
    (median cl7-cl4 RMSD 3.49 A) and cannot separate the clusters: a degenerate reference. Panel B therefore uses
    the complete-H12 medoids of Y10 (cl7 5TM1|A|1 = the S12 cl7 medoid, cl4 7UJ7|C|1). The S12-medoid distances
    stay in the CSV (columns *_s12) and the degeneracy numbers in numbers.json (medoid_check).
  Hypergeometric tests (two-sided, minimum-likelihood = Fisher exact on the 2 x 2 table group x outlier), at two
    unit levels (rule 1.4):
      chain level  units = the 788 matrix chains; outlier unit = consensus outlier chain.
      entry level  units = (entry, cluster group) pairs; an entry with chains in two groups counts once per group;
                   outlier unit = the pair holds at least one consensus outlier chain.
    Tests per level: cl7 vs all, cl4 vs all, minor vs all, and cl7 within the main clusters only (cl7 vs cl4, the
    minor units removed; cl4 is the complement, so this one test covers both). BH over all 8 tests (rule 1.4).
    Per S08 cluster the observed and expected counts are listed descriptively (rows with tested = False).
  Verdict rule: a group deviates from expectation when q < 0.05 at chain AND entry level ('pattern'); at one level
    only = 'chain level only' / 'entry level only'; otherwise 'no pattern'. The cl7-vs-cl4 test answers whether
    outliers sit in one main cluster. Minor-cluster chains are far from both main clusters by definition, which is
    part of what the Y15 detectors measure (medoid distance, kNN, Mahalanobis to the nearer main cluster), so a minor
    enrichment is reported as expected by construction, not as a finding.
Outputs: results/sv/SV06/outliers_on_map.csv, hypergeom_tests.csv, numbers.json, verdict.md, env.txt.
Run: .venv/bin/python -B -I scripts/sv/SV06_outliers.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, hypergeom

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV06"
MAIN = (7, 4)
ALPHA = 0.05


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1)
    return q


def group_of(c: int) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(c), "minor")


def hyper_row(level, test, universe, n_univ, k_group, n_out, x_obs):
    """Two-sided hypergeometric (Fisher exact, minlike) for x_obs outliers among k_group units of the group."""
    table = [[x_obs, k_group - x_obs], [n_out - x_obs, (n_univ - k_group) - (n_out - x_obs)]]
    p2 = fisher_exact(table, alternative="two-sided")[1]
    rv = hypergeom(n_univ, k_group, n_out)
    exp = n_out * k_group / n_univ
    return dict(level=level, test=test, universe=universe, tested=True, n_units=n_univ, n_group=k_group,
                n_outliers=n_out, observed=x_obs, expected=exp, obs_over_exp=x_obs / exp if exp else np.nan,
                p_enriched=float(rv.sf(x_obs - 1)), p_depleted=float(rv.cdf(x_obs)), p_two_sided=float(p2))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = np.load(ROOT / "results/rmsd_entry.npy")
    idx = pd.read_csv(ROOT / "results/rmsd_entry_index.csv")
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    assert d.shape == (788, 788) and len(idx) == 788
    assert (idx.chain_key.values == ct.chain_key.values).all()
    of = pd.read_csv(ROOT / "results/ext/outlier_flags.csv")
    ov = pd.read_csv(ROOT / "results/ext/outlier_verdicts.csv")
    emb = pd.read_csv(ROOT / "results/ext/embedding_rmsd.csv")
    leaf = pd.read_csv(ROOT / "results/ext/leaf_order.csv")
    med = pd.read_csv(ROOT / "results/medoids.csv")
    ec = pd.read_csv(ROOT / "results/entry_clusters.csv")

    # consistency checks of the frozen inputs
    chk = ct[["chain_key", "s08_cluster"]].merge(ec[["chain_key", "cut_1.5"]], on="chain_key", validate="one_to_one")
    assert (chk.s08_cluster == chk["cut_1.5"]).all(), "SV00 s08_cluster differs from S08 cut_1.5"
    cons = of[of.consensus_outlier].chain_key
    assert set(cons) == set(ov.chain_key) and len(cons) == 22, "Y15 consensus set != Y16 verdict set"
    y16 = json.loads((ROOT / "results/ext/numbers_Y16.json").read_text())
    y16v = {v["chain_key"]: (v["verdict"], v["reasons"]) for v in y16["verdicts"]}
    assert y16v == {r.chain_key: (r.verdict, r.reasons) for r in ov.itertuples()}, "verdict CSV != numbers_Y16.json"
    assert (emb.chain_key.values == idx.chain_key.values).all()
    assert (of.chain_key.values == idx.chain_key.values).all()
    assert len(leaf) == 788 and leaf.chain_key.is_unique

    prof = pd.read_csv(ROOT / "results/ext/cluster_displacement_profile.csv")

    row = {k: i for i, k in enumerate(idx.chain_key)}
    med_row = {int(r.cluster): row[r.chain_key] for r in med.itertuples()}      # S12 medoids, every S08 cluster
    med_key = {int(r.cluster): r.chain_key for r in med.itertuples()}
    s7, s4 = med_row[7], med_row[4]
    ref_key = {7: prof.medoid_cl7.iloc[0], 4: prof.medoid_cl4.iloc[0]}          # Y10 complete-H12 medoids
    assert prof.medoid_cl7.nunique() == 1 and prof.medoid_cl4.nunique() == 1
    m7, m4 = row[ref_key[7]], row[ref_key[4]]
    cl = ct.s08_cluster.values
    medoid_check = dict(
        s12_cl7=med_key[7], s12_cl4=med_key[4],
        s12_cl7_h12_label=str(ct.h12_label.iloc[s7]), s12_cl4_h12_label=str(ct.h12_label.iloc[s4]),
        rmsd_s12_cl7_cl4=float(d[s7, s4]),
        y10_cl7=ref_key[7], y10_cl4=ref_key[4],
        y10_cl7_h12_label=str(ct.h12_label.iloc[m7]), y10_cl4_h12_label=str(ct.h12_label.iloc[m4]),
        rmsd_y10_cl7_cl4=float(d[m7, m4]),
        median_rmsd_cl7_cl4_chains=float(np.median(d[np.ix_(cl == 7, cl == 4)])),
        median_rmsd_s12_cl4_medoid_to_cl7=float(np.median(d[s4, cl == 7])),
        median_rmsd_s12_cl4_medoid_to_cl4=float(np.median(d[s4, cl == 4])),
        median_rmsd_y10_cl4_medoid_to_cl7=float(np.median(d[m4, cl == 7])),
        median_rmsd_y10_cl4_medoid_to_cl4=float(np.median(d[m4, cl == 4])))

    # ------------------------------------------------------------------ per-outlier table
    t = ov[["chain_key", "entry_id", "verdict", "reasons", "n_methods"]].rename(columns={"verdict": "reason"})
    t = t.merge(ct[["chain_key", "auth_asym_id", "s08_cluster"]], on="chain_key", validate="one_to_one")
    t = t.merge(emb[["chain_key", "mds1", "mds2"]], on="chain_key", validate="one_to_one")
    t = t.merge(leaf[["chain_key", "leaf_rank"]], on="chain_key", validate="one_to_one")
    t["label"] = t.entry_id + "|" + t.auth_asym_id.astype(str)
    t["cluster_group"] = t.s08_cluster.map(group_of)
    r = t.chain_key.map(row).values
    t["rmsd_cl7_medoid"] = d[r, m7]
    t["rmsd_cl4_medoid"] = d[r, m4]
    main_c = np.where(t.s08_cluster.isin(MAIN), t.s08_cluster, np.where(d[r, m7] < d[r, m4], 7, 4))
    t["own_main_cluster"] = main_c
    home = of.set_index("chain_key").home_cluster
    assert (t.chain_key.map(home).values == main_c).all(), "own main cluster != Y15 home_cluster"
    t["own_main_rule"] = np.where(t.s08_cluster.isin(MAIN), "S08 cluster", "nearer main medoid")
    t["own_medoid"] = [ref_key[c] for c in main_c]
    t["other_medoid"] = [ref_key[4 if c == 7 else 7] for c in main_c]
    t["rmsd_own_medoid"] = np.where(main_c == 7, t.rmsd_cl7_medoid, t.rmsd_cl4_medoid)
    t["rmsd_other_medoid"] = np.where(main_c == 7, t.rmsd_cl4_medoid, t.rmsd_cl7_medoid)
    t["own_other_ratio"] = t.rmsd_own_medoid / t.rmsd_other_medoid
    t["closer_to_other_medoid"] = t.rmsd_other_medoid < t.rmsd_own_medoid
    t["own_s08_medoid"] = t.s08_cluster.map(med_key)
    t["rmsd_own_s08_medoid"] = [d[i, med_row[int(c)]] for i, c in zip(r, t.s08_cluster)]
    t["rmsd_cl7_medoid_s12"] = d[r, s7]
    t["rmsd_cl4_medoid_s12"] = d[r, s4]
    t["rmsd_own_medoid_s12"] = np.where(main_c == 7, t.rmsd_cl7_medoid_s12, t.rmsd_cl4_medoid_s12)
    t["rmsd_other_medoid_s12"] = np.where(main_c == 7, t.rmsd_cl4_medoid_s12, t.rmsd_cl7_medoid_s12)
    cols = ["chain_key", "label", "entry_id", "auth_asym_id", "s08_cluster", "cluster_group", "n_methods", "reason",
            "reasons", "mds1", "mds2", "leaf_rank", "own_main_cluster", "own_main_rule", "own_medoid", "other_medoid",
            "rmsd_own_medoid", "rmsd_other_medoid", "own_other_ratio", "closer_to_other_medoid", "rmsd_cl7_medoid",
            "rmsd_cl4_medoid", "own_s08_medoid", "rmsd_own_s08_medoid", "rmsd_own_medoid_s12",
            "rmsd_other_medoid_s12", "rmsd_cl7_medoid_s12", "rmsd_cl4_medoid_s12"]
    t = t[cols].sort_values(["cluster_group", "leaf_rank"],
                            key=lambda s: s.map({"cl7": 0, "cl4": 1, "minor": 2}) if s.name == "cluster_group" else s)
    t.to_csv(OUT / "outliers_on_map.csv", index=False)

    # ------------------------------------------------------------------ hypergeometric tests
    u = ct[["chain_key", "entry_id", "s08_cluster"]].copy()
    u["group"] = u.s08_cluster.map(group_of)
    u["outlier"] = u.chain_key.isin(cons)
    ent = u.groupby(["entry_id", "group"], as_index=False).outlier.any()
    rows = []
    for level, df in (("chain", u), ("entry", ent)):
        n_all, n_out = len(df), int(df.outlier.sum())
        for g in ("cl7", "cl4", "minor"):
            m = df.group == g
            rows.append(hyper_row(level, f"{g} vs all", "all units", n_all, int(m.sum()), n_out,
                                  int((m & df.outlier).sum())))
        mm = df[df.group.isin(("cl7", "cl4"))]
        m = mm.group == "cl7"
        rows.append(hyper_row(level, "cl7 vs cl4 (main only)", "cl7 + cl4 units", len(mm), int(m.sum()),
                              int(mm.outlier.sum()), int((m & mm.outlier).sum())))
    ht = pd.DataFrame(rows)
    ht["q_bh"] = bh(ht.p_two_sided.values)
    ht["direction"] = np.where(ht.observed > ht.expected, "more than expected",
                               np.where(ht.observed < ht.expected, "fewer than expected", "as expected"))
    # descriptive per S08 cluster (chain level and entry level), not tested
    desc = []
    for level, df, key in (("chain", u, "s08_cluster"),
                           ("entry", u.groupby(["entry_id", "s08_cluster"], as_index=False).outlier.any(),
                            "s08_cluster")):
        n_all, n_out = len(df), int(df.outlier.sum())
        for c in sorted(df[key].unique()):
            m = df[key] == c
            k = int(m.sum())
            desc.append(dict(level=level, test=f"S08 cluster {c} ({group_of(c)})", universe="all units",
                             tested=False, n_units=n_all, n_group=k, n_outliers=n_out,
                             observed=int((m & df.outlier).sum()), expected=n_out * k / n_all,
                             obs_over_exp=(int((m & df.outlier).sum()) / (n_out * k / n_all))))
    ht = pd.concat([ht, pd.DataFrame(desc)], ignore_index=True)
    ht.to_csv(OUT / "hypergeom_tests.csv", index=False)

    # ------------------------------------------------------------------ verdict
    tested = ht[ht.tested]

    def level_call(test):
        s = tested[tested.test == test].set_index("level")
        sig = {lv: bool(s.loc[lv, "q_bh"] < ALPHA) for lv in ("chain", "entry")}
        if sig["chain"] and sig["entry"]:
            return "pattern", s
        if sig["chain"] or sig["entry"]:
            return ("chain level only" if sig["chain"] else "entry level only"), s
        return "no pattern", s

    calls = {}
    for test in ("cl7 vs cl4 (main only)", "cl7 vs all", "cl4 vs all", "minor vs all"):
        c, s = level_call(test)
        calls[test] = dict(call=c, **{f"{lv}_{k}": (float(s.loc[lv, k]) if k not in ("observed", "n_group") else
                                                    int(s.loc[lv, k]))
                                      for lv in ("chain", "entry")
                                      for k in ("observed", "expected", "n_group", "p_two_sided", "q_bh")})

    nmain = t[t.cluster_group != "minor"]
    nums = dict(
        n_chains=788, n_outliers=int(len(t)), n_outlier_entries=int(t.entry_id.nunique()),
        n_entry_units=int(len(ent)), n_entry_units_outlier=int(ent.outlier.sum()),
        outliers_by_group={g: int((t.cluster_group == g).sum()) for g in ("cl7", "cl4", "minor")},
        reason_by_group={g: {rr: int(((t.cluster_group == g) & (t.reason == rr)).sum())
                             for rr in ("quality", "conformation")} for g in ("cl7", "cl4", "minor")},
        chains_by_group={g: int((u.group == g).sum()) for g in ("cl7", "cl4", "minor")},
        medoid_cl7=ref_key[7], medoid_cl4=ref_key[4], rmsd_medoid_cl7_cl4=float(d[m7, m4]),
        medoid_check=medoid_check,
        n_main_closer_to_other_s12=int((nmain.rmsd_other_medoid_s12 < nmain.rmsd_own_medoid_s12).sum()),
        n_main_closer_to_other=int(nmain.closer_to_other_medoid.sum()), n_main=int(len(nmain)),
        ratio_main_min=float(nmain.own_other_ratio.min()), ratio_main_max=float(nmain.own_other_ratio.max()),
        ratio_main_median=float(nmain.own_other_ratio.median()),
        ratio_minor_min=float(t[t.cluster_group == "minor"].own_other_ratio.min()),
        ratio_minor_max=float(t[t.cluster_group == "minor"].own_other_ratio.max()),
        calls=calls, alpha=ALPHA, n_tests_bh=int(len(tested)),
    )
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1))

    def fmt(test):
        c = calls[test]
        return (f"chain {c['chain_observed']} vs {c['chain_expected']:.1f} expected, q {c['chain_q_bh']:.2g}; "
                f"entry {c['entry_observed']} vs {c['entry_expected']:.1f}, q {c['entry_q_bh']:.2g}")

    top = nmain.sort_values("own_other_ratio", ascending=False).iloc[0]
    main_call = calls["cl7 vs cl4 (main only)"]["call"]
    ob = nums["outliers_by_group"]
    lines = [
        f"cl7 vs cl4: {main_call} - of the {nums['n_main']} consensus outliers inside the main clusters, "
        f"{ob['cl7']} are in cl7 and {ob['cl4']} in cl4 ({fmt('cl7 vs cl4 (main only)')}; hypergeometric, BH over "
        f"{nums['n_tests_bh']} tests).",
        f"Minor clusters hold {ob['minor']} of {nums['n_outliers']} outliers from {nums['chains_by_group']['minor']} "
        f"chains ({calls['minor vs all']['call']}: {fmt('minor vs all')}); expected by construction, since the Y15 "
        f"detectors measure distance from the main clusters.",
        f"Panel B (complete-H12 medoids {ref_key[7]} / {ref_key[4]}, {nums['rmsd_medoid_cl7_cl4']:.2f} A apart): "
        f"{nums['n_main_closer_to_other']} of {nums['n_main']} main-cluster outliers are closer to the other medoid "
        f"than to their own (own/other RMSD ratio {nums['ratio_main_min']:.2f}-{nums['ratio_main_max']:.2f}, median "
        f"{nums['ratio_main_median']:.2f}; highest {top.label} {top.own_other_ratio:.2f}, Y16 reason: {top.reasons}); "
        f"the S12 cl4 medoid {med_key[4]} has no modelled H12 and is not used.",
    ]
    (OUT / "verdict.md").write_text("\n".join(lines) + "\n")

    import matplotlib
    tz = timezone(timedelta(hours=5))
    env = [f"# SV06 environment, {datetime.now(tz).strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}",
           f"platform {platform.platform()}",
           f"numpy {np.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}",
           f"matplotlib {matplotlib.__version__}",
           "nothing installed, upgraded or downgraded by SV06; no network used"]
    (OUT / "env.txt").write_text("\n".join(env) + "\n")
    print(t[["label", "s08_cluster", "reason", "rmsd_own_medoid", "rmsd_other_medoid", "own_other_ratio",
             "rmsd_own_s08_medoid", "rmsd_own_medoid_s12", "rmsd_other_medoid_s12"]].to_string())
    print(json.dumps(medoid_check, indent=1))
    print(ht.to_string())
    print("\n".join(lines))


if __name__ == "__main__":
    main()
