"""Tracker Y14: hypothesis tests for the cl7 / cl4 split (test, do not assume).

Rows: the 774 chains in cl7 or cl4. Every hypothesis is tested at the chain level and, as a sensitivity check, at the
entry level (one chain per entry: the entry's first matrix chain; chains of one entry share crystal and ligand) or
the ligand level where the variable belongs to the ligand.
Binary variables: Fisher exact test, odds ratio of being in cl7 with a 95% CI (Woolf, Haldane +0.5).
Continuous variables: Mann-Whitney U, rank-biserial r (positive = larger in cl7) with a 95% bootstrap CI
(1000 resamples of entries, seed 42).
Adjusted questions (logistic regression, statsmodels, likelihood-ratio test):
  Q-a  does quality (resolution, R-free, mean B) add anything once peptide and mutation status are known?
  Q-b  does the ligand's basic amine add anything once peptide and mutation status are known?
  Q-c  Mantel-Haenszel: peptide effect stratified by construct mutation status.
All p-values (chain level) are BH-corrected as one family.
Outputs: results/ext/split_hypotheses.csv, figures/ext/Y14_split_hypotheses.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_split_hypotheses.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2, fisher_exact, mannwhitneyu
from statsmodels.stats.contingency_tables import StratifiedTable

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK, INK2, ORANGE, RES, SEED, VIOLET, bh, chain_table, mpl_setup, save_fig,
                        write_numbers)

warnings.filterwarnings("ignore")


def table():
    t = chain_table()
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status", "peptide_present", "peptide_class"]]
    ld = pd.read_csv(RES / "ligand_descriptors.csv")[["ligand_id", "heavy_atoms", "basic_amine", "mw"]]
    em = pd.read_csv(RES / "entry_metadata.csv")[["entry_id", "release_year", "ph"]]
    t = t.merge(cm, on="chain_key").merge(ld, on="ligand_id", how="left").merge(em, on="entry_id", how="left")
    t = t[t.cluster.isin((7, 4))].reset_index(drop=True)
    t["y"] = (t.cluster == 7).astype(int)
    t["p21"] = t.space_group == "P 1 21 1"
    t["h12_complete"] = t.h12_label == "complete"
    t["tier_A"] = t.tier == "A"
    return t


def fisher_row(name, x, y, level):
    x = np.asarray(x, bool)
    a = int(((x) & (y == 1)).sum()); b = int(((x) & (y == 0)).sum())
    c = int(((~x) & (y == 1)).sum()); d = int(((~x) & (y == 0)).sum())
    _, p = fisher_exact([[a, b], [c, d]])
    lor = np.log((a + .5) * (d + .5) / ((b + .5) * (c + .5)))
    se = np.sqrt(1 / (a + .5) + 1 / (b + .5) + 1 / (c + .5) + 1 / (d + .5))
    return dict(hypothesis=name, test="Fisher exact", level=level, n=len(y), effect_type="log odds ratio (cl7)",
                effect=lor, ci_low=lor - 1.96 * se, ci_high=lor + 1.96 * se, p_value=p,
                detail=f"cl7 with/without: {a}/{c}; cl4 with/without: {b}/{d}")


def mwu_row(name, v, y, groups, level):
    ok = ~np.isnan(v)
    v, y, groups = v[ok], y[ok], groups[ok]
    u, p = mannwhitneyu(v[y == 1], v[y == 0])
    r = 1 - 2 * u / ((y == 1).sum() * (y == 0).sum())
    r = -r  # positive = larger in cl7
    rng = np.random.default_rng(SEED)
    ug = np.unique(groups)
    idx = {g: np.where(groups == g)[0] for g in ug}
    bs = []
    for _ in range(1000):
        s = np.concatenate([idx[g] for g in rng.choice(ug, len(ug))])
        a, b = v[s][y[s] == 1], v[s][y[s] == 0]
        if len(a) and len(b):
            uu = mannwhitneyu(a, b)[0]
            bs.append(-(1 - 2 * uu / (len(a) * len(b))))
    return dict(hypothesis=name, test="Mann-Whitney U", level=level, n=int(ok.sum()), effect_type="rank-biserial (cl7 larger +)",
                effect=r, ci_low=np.percentile(bs, 2.5), ci_high=np.percentile(bs, 97.5), p_value=p,
                detail=f"median cl7 {np.median(v[y == 1]):.4g}, cl4 {np.median(v[y == 0]):.4g}")


def lr_test(t, base, extra, name):
    d = t.dropna(subset=base + extra)
    X0 = sm.add_constant(d[base].astype(float))
    X1 = sm.add_constant(d[base + extra].astype(float))
    m0 = sm.Logit(d.y, X0).fit(disp=0, method="bfgs", maxiter=2000)
    m1 = sm.Logit(d.y, X1).fit(disp=0, method="bfgs", maxiter=2000)
    lr = 2 * (m1.llf - m0.llf)
    p = chi2.sf(lr, len(extra))
    return dict(hypothesis=name, test="logistic LR test", level="chain", n=len(d), effect_type="McFadden R2 gain",
                effect=(m1.llf - m0.llf) / -m0.llnull if m0.llnull else np.nan, ci_low=np.nan, ci_high=np.nan,
                p_value=p, detail=f"LR={lr:.2f}, df={len(extra)}; base pseudo-R2 {m0.prsquared:.3f}, full {m1.prsquared:.3f}")


def main():
    t = table()
    y = t.y.values
    ent = t.drop_duplicates("entry_id")
    lig = t.dropna(subset=["ligand_id"]).groupby("ligand_id").agg(y=("y", "mean"), heavy=("heavy_atoms", "first"),
                                                                  amine=("basic_amine", "first"))
    lig = lig[(lig.y == 0) | (lig.y == 1)]
    rows = []
    for name, col in (("coactivator/other peptide present", "peptide_present"), ("Y537S construct", None),
                      ("L536S construct", None), ("D538G construct", None), ("wild-type LBD", None),
                      ("H12 complete (all 10 C-alpha)", "h12_complete"), ("space group P21", "p21"),
                      ("QC tier A", "tier_A")):
        if col is None:
            key = {"Y537S construct": "Y537S", "L536S construct": "L536S", "D538G construct": "D538G",
                   "wild-type LBD": "wild type"}[name]
            x, xe = t.mutation_status == key, ent.mutation_status == key
        else:
            x, xe = t[col], ent[col]
        rows.append(fisher_row(name, x, y, "chain"))
        rows.append(fisher_row(name, xe, ent.y.values, "entry"))
    rows.append(fisher_row("ligand has a basic amine", t.basic_amine.fillna(False).astype(bool), y, "chain"))
    rows.append(fisher_row("ligand has a basic amine", lig.amine.fillna(False).astype(bool), lig.y.values.astype(int), "ligand"))
    g = t.entry_id.values
    for name, col in (("ligand heavy atoms", "heavy_atoms"), ("resolution (A)", "resolution_A"), ("R-free", "r_free"),
                      ("mean B (A^2)", "mean_b_A2"), ("pocket size (P2Rank SAS points)", "sas_points"),
                      ("release year", "release_year"), ("crystallisation pH", "ph")):
        rows.append(mwu_row(name, t[col].astype(float).values, y, g, "chain"))
        rows.append(mwu_row(name, ent[col].astype(float).values, ent.y.values, ent.entry_id.values, "entry"))
    t["mut_Y537S"] = (t.mutation_status == "Y537S").astype(float)
    t["mut_L536S"] = (t.mutation_status == "L536S").astype(float)
    t["pep"] = t.peptide_present.astype(float)
    t["amine"] = t.basic_amine.astype(float)
    rows.append(lr_test(t, ["pep", "mut_Y537S", "mut_L536S"], ["resolution_A", "r_free", "mean_b_A2"],
                        "quality adds to peptide + mutation"))
    rows.append(lr_test(t, ["pep", "mut_Y537S", "mut_L536S"], ["amine"], "basic amine adds to peptide + mutation"))
    tabs = []
    for s, d in t.groupby("mutation_status"):
        tb = pd.crosstab(d.peptide_present, d.y).reindex(index=[True, False], columns=[1, 0], fill_value=0).values
        if (tb.sum(0) > 0).all() and (tb.sum(1) > 0).all():
            tabs.append(tb + 0.5)
    st = StratifiedTable(tabs)
    rows.append(dict(hypothesis="peptide, stratified by mutation status", test="Mantel-Haenszel", level="chain", n=len(t),
                     effect_type="log odds ratio (cl7)", effect=np.log(st.oddsratio_pooled),
                     ci_low=np.log(st.oddsratio_pooled_confint()[0]), ci_high=np.log(st.oddsratio_pooled_confint()[1]),
                     p_value=st.test_null_odds().pvalue, detail=f"{len(tabs)} strata with both levels"))
    df = pd.DataFrame(rows)
    df["q_value"] = np.nan
    ch = df.level == "chain"
    df.loc[ch, "q_value"] = bh(df.loc[ch, "p_value"])
    df.loc[~ch, "q_value"] = bh(df.loc[~ch, "p_value"])
    df.to_csv(RES / "split_hypotheses.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(19, 6.8), gridspec_kw={"width_ratios": [1.2, 1.1, 0.9]})
    ax = axs[0]
    b = df[(df.effect_type == "log odds ratio (cl7)")].copy()
    b["lab"] = b.hypothesis + " [" + b.level + "]"
    b = b.iloc[::-1]
    yy = np.arange(len(b))
    col = [VIOLET if l == "chain" else GREY for l in b.level]
    ax.errorbar(b.effect, yy, xerr=[b.effect - b.ci_low, b.ci_high - b.effect], fmt="none", ecolor=col, lw=1.5)
    ax.scatter(b.effect, yy, c=col, s=30, zorder=3)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_yticks(yy, b.lab, fontsize=7.5)
    ax.set(xlabel="log odds ratio of being in cl7 (95% CI)", title="Binary hypotheses (violet = chain, grey = entry/ligand)")
    ax = axs[1]
    c = df[df.effect_type.str.startswith("rank")].copy()
    c["lab"] = c.hypothesis + " [" + c.level + "]"
    c = c.iloc[::-1]
    yy = np.arange(len(c))
    col = [VIOLET if l == "chain" else GREY for l in c.level]
    ax.errorbar(c.effect, yy, xerr=[c.effect - c.ci_low, c.ci_high - c.effect], fmt="none", ecolor=col, lw=1.5)
    ax.scatter(c.effect, yy, c=col, s=30, zorder=3)
    ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_yticks(yy, c.lab, fontsize=7.5)
    ax.set(xlabel="rank-biserial r (+ = larger in cl7), 95% entry-bootstrap CI", title="Continuous hypotheses")
    ax = axs[2]
    d = df[ch].sort_values("q_value").head(16)[::-1]
    ax.barh(d.hypothesis.str.slice(0, 38), -np.log10(d.q_value.clip(lower=1e-300)), color=GREY)
    ax.axvline(-np.log10(0.05), color=INK, ls="--", lw=1)
    ax.set(xlabel="-log10 q (BH, chain-level family); q = 0 (underflow) drawn at 300", title="Strength of evidence")
    r = df.set_index(["hypothesis", "level"])
    qa = r.loc[("quality adds to peptide + mutation", "chain")]
    pe = r.loc[("coactivator/other peptide present", "entry")]
    fig.suptitle(f"Hypotheses: peptide (entry-level log OR {pe.effect:.1f}) and construct mutations explain the split; "
                 f"quality adds little once they are known (McFadden gain {qa.effect:.3f}, p={qa.p_value:.2g})",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y14_split_hypotheses", section="C. Why two clusters", task="Y14",
             shows="forest plots of binary (log OR) and continuous (rank-biserial) hypotheses at chain and entry level; BH q",
             takeaway=(f"Peptide presence, Y537S/L536S constructs and the ligand basic amine have the largest effects; "
                       f"adding resolution/R-free/mean B to peptide+mutation gains McFadden R2 {qa.effect:.3f} "
                       f"(LR p={qa.p_value:.2g})."),
             source="results/ext/split_hypotheses.csv", script="scripts/ext_split_hypotheses.py")
    write_numbers("Y14", df.set_index(df.hypothesis + "|" + df.level)[["effect", "ci_low", "ci_high", "p_value",
                                                                         "q_value", "n"]].round(6).to_dict("index"))
    print(df[["hypothesis", "level", "n", "effect", "ci_low", "ci_high", "p_value", "q_value", "detail"]].round(4).to_string())


if __name__ == "__main__":
    main()
