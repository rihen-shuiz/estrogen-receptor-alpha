"""SV03 steps 4-5 (RUN_SV1028): cl7 vs cl4 per descriptor, within-cluster associations, "more sensitive than RMSD".

Input: results/sv/SV03/geometry_per_chain.csv (SV03_geometry.py; covariates copied from results/sv/SV00/chain_table.csv).
Rules fixed before any test was run:
- Descriptors (16): uncensored lengths of the 8 helices (H12 = 538-547, run file section 2), angles H12-H3,
  H12-H11, H12-H4/5, H11-H3, distances H12 centroid - H3 axis and - H11 axis, H11 kink, H11-H12 hinge.
  "involves_H12" marks the 7 descriptors that use H12; the other 9 describe the other helices.
- Rank tests at two levels (rule 1.4; the SV11 / SV10 conventions of this run):
  chain level = GEE (Gaussian, identity, independence working correlation, clusters = entries) of the average rank
    of the descriptor (ranked among the compared chains) on the group indicator (binary) or on the average rank of
    the covariate (continuous); cluster-robust Wald p. Effect size: rank-biserial r = 2U / (n1 n0) - 1 from the
    Mann-Whitney U of the chains (positive = higher in the second group: cl4 in (a), covariate True in (b)); for the
    continuous covariate (ligand heavy atoms) Spearman rho of the chains.
  entry level = two-sided Mann-Whitney U (binary) or Spearman (continuous) on entry values; same effect sizes.
- (a) cl7 vs cl4: chains with s08_cluster 7 or 4. Entry level: entry group by the SV10 rule (>= 80 % of the entry's
  kept chains in cl7 or in cl4; other entries dropped); entry value = median over the entry's chains of that cluster.
- (b) within cl7 and within cl4: covariates ligand heavy atoms (continuous), basic amine, carboxylic acid (SV00
  flags of the primary ligand of interest), peptide (chain-level peptide_bound), construct mutations coded as in
  SV01a: Y537S, L536S, D538G, other (any other construct mutation). Chains without a ligand of interest are left out
  of the three ligand tests. Entry level: unit = entry x cluster; value = median over the unit's chains that have the
  descriptor; covariate over the same chains = median (heavy atoms) or True when >= 50 % of them are True.
- Testable: binary - each group >= 2 chains from >= 2 entries (chain level), >= 2 entries (entry level); continuous -
  >= 3 units from >= 2 entries; the descriptor and the covariate vary among the units. Otherwise not tested (no p,
  left out of BH).
- BH (rule 1.4) within each family: (a) chain, (a) entry, (b) chain, (b) entry.
- Step 5 [DEFAULT] "more sensitive than RMSD": a (b) comparison (descriptor x covariate x cluster) with |r| >= 0.3
  (|rho| >= 0.3 for heavy atoms, the continuous analogue) and q < 0.05 at both the chain and the entry level.
  None -> "no added sensitivity".
Outputs: results/sv/SV03/tests.csv (one row per test), sensitivity_check.csv (step 5, one row per (b) comparison with
  both levels), cl7_vs_cl4_wide.csv (step 4a, both levels side by side), numbers_tests.json.
Run: .venv/bin/python -B -I scripts/sv/SV03_tests.py
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import mannwhitneyu, rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV03"
R_MIN, Q_MAX = 0.3, 0.05       # step 5 [DEFAULT]
ENTRY_SHARE_MIN = 0.8          # SV10 entry-group rule
BIN_SHARE = 0.5                # entry-level binary covariate: >= 50 % of the unit's chains

DESCRIPTORS = {
    "ulen_H1": "uncensored length H1 (res)", "ulen_H3": "uncensored length H3 (res)",
    "ulen_H45": "uncensored length H4/5 (res)", "ulen_H7": "uncensored length H7 (res)",
    "ulen_H8": "uncensored length H8 (res)", "ulen_H9": "uncensored length H9 (res)",
    "ulen_H1011": "uncensored length H10/11 (res)", "ulen_H12": "uncensored length H12 (res)",
    "angle_H12_H3_deg": "angle H12-H3 (deg)", "angle_H12_H11_deg": "angle H12-H11 (deg)",
    "angle_H12_H45_deg": "angle H12-H4/5 (deg)", "angle_H11_H3_deg": "angle H11-H3 (deg)",
    "dist_H12c_H3axis_A": "H12 centroid to H3 axis (A)", "dist_H12c_H11axis_A": "H12 centroid to H11 axis (A)",
    "kink_H11_deg": "H11 kink angle (deg)", "hinge_H11C_H12_deg": "H11-H12 hinge angle (deg)"}
INVOLVES_H12 = {d: ("H12" in d) for d in DESCRIPTORS}
COVARIATES = {"ligand_heavy_atoms": ("continuous", "ligand heavy atoms"),
              "ligand_basic_amine": ("binary", "basic amine"),
              "ligand_carboxylic_acid": ("binary", "carboxylic acid"),
              "peptide_bound": ("binary", "peptide bound"),
              "mut_Y537S": ("binary", "Y537S"), "mut_L536S": ("binary", "L536S"),
              "mut_D538G": ("binary", "D538G"), "mut_other": ("binary", "other construct mutation")}
LIGAND_COVS = ("ligand_heavy_atoms", "ligand_basic_amine", "ligand_carboxylic_acid")


def to_bool(s):
    return s.map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0, 1: 1.0, 0: 0.0})


def gee_p(y, x, groups):
    X = sm.add_constant(np.asarray(x, float))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = sm.GEE(np.asarray(y, float), X, groups=np.asarray(groups), family=sm.families.Gaussian(),
                     cov_struct=sm.cov_struct.Independence()).fit()
    lo, hi = res.conf_int()[1]
    return float(res.params[1]), float(lo), float(hi), float(res.params[1] / res.bse[1]), float(res.pvalues[1])


def rank_biserial(a, b):
    """r > 0: values in b tend to be higher than in a."""
    u = mannwhitneyu(b, a, alternative="two-sided")
    return 2 * float(u.statistic) / (len(a) * len(b)) - 1, float(u.statistic), float(u.pvalue)


def iqr_text(x):
    x = pd.Series(x).dropna()
    return f"{x.quantile(0.25):.2f}-{x.quantile(0.75):.2f}" if len(x) else ""


def binary_test(d, val, grp, entry, level):
    """d: units with val (descriptor), grp (0/1), entry. Returns dict."""
    a, b = d[d[grp] == 0], d[d[grp] == 1]
    r = {"n0": len(a), "n1": len(b), "entries0": a[entry].nunique(), "entries1": b[entry].nunique(),
         "median0": a[val].median(), "median1": b[val].median(), "iqr0": iqr_text(a[val]), "iqr1": iqr_text(b[val])}
    need = 2
    if min(len(a), len(b)) < 2 or min(r["entries0"], r["entries1"]) < need:
        return {**r, "tested": False, "note": "fewer than 2 units or 2 entries in a group"}
    if d[val].nunique() < 2:
        return {**r, "tested": False, "note": "no variation in the descriptor"}
    rb, u, p_mw = rank_biserial(a[val].to_numpy(), b[val].to_numpy())
    r.update(effect=rb, effect_type="rank-biserial r", U=u)
    if level == "chain":
        coef, lo, hi, z, p = gee_p(rankdata(d[val]), d[grp], d[entry])
        r.update(tested=True, test="GEE on average ranks, entry clusters (robust Wald)", gee_coef_rank=coef,
                 gee_ci_low=lo, gee_ci_high=hi, statistic=z, p=p, p_mannwhitney_chains=p_mw)
    else:
        r.update(tested=True, test="Mann-Whitney U (two-sided) on entry values", statistic=u, p=p_mw)
    return r


def continuous_test(d, val, cov, entry, level):
    r = {"n0": np.nan, "n1": len(d), "entries1": d[entry].nunique()}
    if len(d) < 3 or d[entry].nunique() < 2:
        return {**r, "tested": False, "note": "fewer than 3 units or 2 entries"}
    if d[val].nunique() < 2 or d[cov].nunique() < 2:
        return {**r, "tested": False, "note": "no variation"}
    rho, p_sp = spearmanr(d[val], d[cov])
    r.update(effect=float(rho), effect_type="Spearman rho", median1=d[val].median())
    if level == "chain":
        coef, lo, hi, z, p = gee_p(rankdata(d[val]), rankdata(d[cov]), d[entry])
        r.update(tested=True, test="GEE rank(descriptor) ~ rank(covariate), entry clusters (robust Wald)",
                 gee_coef_rank=coef, gee_ci_low=lo, gee_ci_high=hi, statistic=z, p=p, p_spearman_chains=float(p_sp))
    else:
        r.update(tested=True, test="Spearman on entry values", statistic=float(rho), p=float(p_sp))
    return r


def entry_groups(ct):
    eg = {}
    for e, g in ct.groupby("entry_id"):
        c = g.s08_cluster.to_numpy()
        eg[e] = "cl7" if np.mean(c == 7) >= ENTRY_SHARE_MIN else ("cl4" if np.mean(c == 4) >= ENTRY_SHARE_MIN else "mixed")
    return pd.Series(eg, name="entry_group")


def part_a(w, eg):
    rows = []
    ch = w[w.group.isin(["cl7", "cl4"])].copy()
    ch["is_cl4"] = (ch.group == "cl4").astype(float)
    for dname in DESCRIPTORS:
        base = {"part": "a", "comparison": "cl4 vs cl7", "cluster": "cl7+cl4", "descriptor": dname,
                "descriptor_label": DESCRIPTORS[dname], "involves_H12": INVOLVES_H12[dname],
                "covariate": "cluster (cl4 = 1)", "group0": "cl7", "group1": "cl4"}
        d = ch[ch[dname].notna()]
        rows.append({**base, "level": "chain", **binary_test(d, dname, "is_cl4", "entry_id", "chain")})
        e = d.merge(eg.rename("entry_group"), left_on="entry_id", right_index=True)
        e = e[e.group == e.entry_group]
        ev = e.groupby(["entry_id", "entry_group"])[dname].median().reset_index()
        ev["is_cl4"] = (ev.entry_group == "cl4").astype(float)
        rows.append({**base, "level": "entry", **binary_test(ev, dname, "is_cl4", "entry_id", "entry")})
    return rows


def part_b(w):
    rows = []
    for cl in ("cl7", "cl4"):
        c = w[w.group == cl]
        for dname in DESCRIPTORS:
            for cov, (kind, clab) in COVARIATES.items():
                base = {"part": "b", "comparison": f"{clab} within {cl}", "cluster": cl, "descriptor": dname,
                        "descriptor_label": DESCRIPTORS[dname], "involves_H12": INVOLVES_H12[dname],
                        "covariate": cov}
                d = c[c[dname].notna()].copy()
                if cov in LIGAND_COVS:
                    d = d[d.ligand_id_primary.notna()]
                if kind == "binary":
                    d["x"] = to_bool(d[cov])
                    d = d[d.x.notna()]
                    base.update(group0=f"no {clab}", group1=clab)
                    rows.append({**base, "level": "chain", **binary_test(d, dname, "x", "entry_id", "chain")})
                    ev = d.groupby("entry_id").agg(v=(dname, "median"), share=("x", "mean")).reset_index()
                    ev["x"] = (ev.share >= BIN_SHARE).astype(float)
                    ev = ev.rename(columns={"v": dname})
                    rows.append({**base, "level": "entry", **binary_test(ev, dname, "x", "entry_id", "entry")})
                else:
                    d["x"] = pd.to_numeric(d[cov], errors="coerce")
                    d = d[d.x.notna()]
                    base.update(group0="", group1="")
                    rows.append({**base, "level": "chain", **continuous_test(d, dname, "x", "entry_id", "chain")})
                    ev = d.groupby("entry_id").agg(v=(dname, "median"), x=("x", "median")).reset_index()
                    ev = ev.rename(columns={"v": dname})
                    rows.append({**base, "level": "entry", **continuous_test(ev, dname, "x", "entry_id", "entry")})
    return rows


def main():
    w = pd.read_csv(OUT / "geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    eg = entry_groups(ct)
    t = pd.DataFrame(part_a(w, eg) + part_b(w))
    t["tested"] = t.tested.fillna(False).astype(bool)
    t["q"] = np.nan
    t["bh_family"] = "(" + t.part + ") " + t.level
    for fam, idx in t[t.tested].groupby("bh_family").groups.items():
        t.loc[idx, "q"] = multipletests(t.loc[idx, "p"], method="fdr_bh")[1]
    cols = ["part", "comparison", "cluster", "descriptor", "descriptor_label", "involves_H12", "covariate", "level",
            "group0", "group1", "n0", "n1", "entries0", "entries1", "median0", "median1", "iqr0", "iqr1",
            "effect_type", "effect", "U", "gee_coef_rank", "gee_ci_low", "gee_ci_high", "statistic", "p", "q",
            "bh_family", "tested", "test", "note", "p_mannwhitney_chains", "p_spearman_chains"]
    t = t.reindex(columns=cols)
    t.to_csv(OUT / "tests.csv", index=False)

    # step 5 [DEFAULT] per (b) comparison: both levels
    b = t[t.part == "b"]
    key = ["cluster", "descriptor", "covariate"]
    wide = b[b.level == "chain"].set_index(key)[["effect_type", "n0", "n1", "entries0", "entries1", "median0",
                                                 "median1", "effect", "p", "q", "tested", "involves_H12"]].join(
        b[b.level == "entry"].set_index(key)[["n0", "n1", "median0", "median1", "effect", "p", "q", "tested"]],
        lsuffix="_chain", rsuffix="_entry").reset_index()
    wide["meets_default"] = (wide.tested_chain & wide.tested_entry & (wide.effect_chain.abs() >= R_MIN)
                             & (wide.effect_entry.abs() >= R_MIN) & (wide.q_chain < Q_MAX) & (wide.q_entry < Q_MAX))
    wide["same_sign"] = np.sign(wide.effect_chain) == np.sign(wide.effect_entry)
    wide["descriptor_label"] = wide.descriptor.map(DESCRIPTORS)
    wide.to_csv(OUT / "sensitivity_check.csv", index=False)

    a = t[t.part == "a"]
    aw = a[a.level == "chain"].set_index("descriptor")[["effect", "q", "median0", "median1", "n0", "n1"]].join(
        a[a.level == "entry"].set_index("descriptor")[["effect", "q", "median0", "median1", "n0", "n1"]],
        lsuffix="_chain", rsuffix="_entry")
    aw["separates_both_levels"] = (aw.q_chain < Q_MAX) & (aw.q_entry < Q_MAX)
    aw.to_csv(OUT / "cl7_vs_cl4_wide.csv")
    hits = wide[wide.meets_default]
    num = {"n_tests": t.groupby("bh_family").tested.agg(["size", "sum"]).astype(int).to_dict(orient="index"),
           "entry_groups": eg.value_counts().to_dict(),
           "a_separates_both_levels": aw.index[aw.separates_both_levels].tolist(),
           "a_not_separating": aw.index[~aw.separates_both_levels].tolist(),
           "b_comparisons": int(len(wide)), "b_tested_both_levels": int((wide.tested_chain & wide.tested_entry).sum()),
           "b_q_lt_005_both": int(((wide.q_chain < Q_MAX) & (wide.q_entry < Q_MAX)).sum()),
           "b_meets_default": int(len(hits)),
           "b_hits": hits[key + ["effect_chain", "effect_entry", "q_chain", "q_entry", "same_sign"]].to_dict(orient="records"),
           "verdict": ("more sensitive than RMSD: " + "; ".join(
               f"{r.descriptor} x {r.covariate} in {r.cluster}" for r in hits.itertuples())) if len(hits)
           else "no added sensitivity"}
    (OUT / "numbers_tests.json").write_text(json.dumps(num, indent=1, default=lambda v: v.item()
                                                       if hasattr(v, "item") else str(v)))
    pd.set_option("display.width", 250)
    print(aw.round(4).to_string())
    print(json.dumps({k: num[k] for k in ("n_tests", "entry_groups", "b_comparisons", "b_tested_both_levels",
                                         "b_q_lt_005_both", "b_meets_default", "verdict")}, indent=1, default=str))
    print(wide.sort_values("q_chain").head(25)[key + ["n0_chain", "n1_chain", "effect_chain", "q_chain",
                                                     "n0_entry", "n1_entry", "effect_entry", "q_entry",
                                                     "meets_default"]].round(4).to_string())


if __name__ == "__main__":
    main()
