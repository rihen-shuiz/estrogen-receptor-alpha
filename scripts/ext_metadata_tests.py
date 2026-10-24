"""Tracker Y17: experimental metadata vs structural cluster and vs outlier status.

Inputs: results/ext/entry_metadata.csv, results/ext/chain_metadata.csv (Y02), results/ext/outlier_flags.csv (Y15),
chain_table(). Targets: (1) structural cluster, 3 levels (cl7 / cl4 / minor); (2) consensus outlier (yes/no).
Levels: chain (788 rows) and entry (one row per entry: the first matrix chain, 351 rows).
Variables (missing = 'not_reported' for categorical, kept as its own level; numeric rows with no value are left out
of that one test only): categorical - crystallisation method, precipitant class, crystal system, space group (8 most
common, else other), expression system, oligomer state, X-ray source, peptide class, construct mutation, bound
additive present, PEG/ammonium-sulfate/DMSO/glycerol flags; numeric - pH, growth temperature, cell volume, ASU
chains, polymer chains, release year, diffraction temperature, PEG molecular weight, completeness, Wilson B, R-work.
Tests: categorical - chi-square, effect Cramer's V; numeric - Kruskal-Wallis, effect epsilon^2 = H / (n - 1);
mutual information (numeric binned into quintiles) for every variable. BH per target x level.
Experimental method is constant (all X-ray diffraction) and is reported, not tested.
Outputs: results/ext/metadata_tests.csv, figures/ext/Y17_metadata_effects.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_metadata_tests.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kruskal
from sklearn.metrics import mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK, INK2, ORANGE, RED, RES, VIOLET, bh, chain_table, cramers_v, mpl_setup,
                        save_fig, write_numbers)

CAT = ["cryst_method", "additive_class", "crystal_system", "space_group8", "expression_system", "oligomer_state",
       "source", "peptide_class", "mutation_status", "bound_additive", "has_peg", "has_ammonium_sulfate", "has_dmso",
       "has_glycerol"]
NUM = ["ph", "temp_K", "cell_volume_A3", "chains_asu", "polymer_chains_asu", "release_year", "diffraction_temp_K",
       "peg_mw", "completeness_pct", "wilson_b", "r_work"]


def table():
    t = chain_table()
    em = pd.read_csv(RES / "entry_metadata.csv").drop(columns=["mutation_status", "mutations"])
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status"]]
    o = pd.read_csv(RES / "outlier_flags.csv")[["chain_key", "consensus_outlier"]]
    t = t.drop(columns=["space_group"]).merge(em, on="entry_id", how="left").merge(cm, on="chain_key").merge(o, on="chain_key")
    sg = t.space_group.fillna("not_reported")
    t["space_group8"] = sg.where(sg.isin(sg.value_counts().index[:8]), "other")
    t["bound_additive"] = t.n_bound_additive_types > 0
    t["cluster3"] = np.where(t.cluster == 7, "cl7", np.where(t.cluster == 4, "cl4", "minor"))
    for c in CAT:
        t[c] = t[c].astype(str).replace({"nan": "not_reported"})
    return t


def tests(d, target, level):
    rows = []
    y = d[target].astype(str)
    for v in CAT:
        ct = pd.crosstab(d[v], y)
        if ct.shape[0] < 2:
            continue
        chi, p, dof, _ = chi2_contingency(ct)
        rows.append(dict(variable=v, kind="categorical", target=target, level=level, test="chi-square", n=len(d),
                         n_levels=ct.shape[0], statistic=chi, effect=cramers_v(ct.values), effect_type="Cramer's V",
                         p_value=p, mutual_info=mutual_info_score(d[v], y)))
    for v in NUM:
        dd = d[[v, target]].dropna()
        groups = [g[v].values for _, g in dd.groupby(target) if len(g) >= 2]
        if len(groups) < 2 or dd[v].nunique() < 2:
            continue
        h, p = kruskal(*groups)
        b = pd.qcut(dd[v].rank(method="first"), 5, labels=False)
        rows.append(dict(variable=v, kind="numeric", target=target, level=level, test="Kruskal-Wallis", n=len(dd),
                         n_levels=np.nan, statistic=h, effect=h / (len(dd) - 1), effect_type="epsilon^2", p_value=p,
                         mutual_info=mutual_info_score(b, dd[target].astype(str))))
    return rows


def main():
    t = table()
    ent = t.drop_duplicates("entry_id")
    rows = []
    for target in ("cluster3", "consensus_outlier"):
        rows += tests(t, target, "chain")
        rows += tests(ent, target, "entry")
    df = pd.DataFrame(rows)
    df["q_value"] = np.nan
    for (tg, lv), g in df.groupby(["target", "level"]):
        df.loc[g.index, "q_value"] = bh(g.p_value)
    df["target"] = df.target.replace({"cluster3": "structural cluster", "consensus_outlier": "consensus outlier"})
    df.to_csv(RES / "metadata_tests.csv", index=False)

    plt = mpl_setup()
    fig = plt.figure(figsize=(19, 8.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.1, 1])
    for k, tg in enumerate(("structural cluster", "consensus outlier")):
        ax = fig.add_subplot(gs[k])
        e = df[(df.target == tg) & (df.level == "entry")].set_index("variable")
        c = df[(df.target == tg) & (df.level == "chain")].set_index("variable")
        order = c.sort_values("effect").index
        yy = np.arange(len(order))
        ax.barh(yy + 0.2, c.loc[order, "effect"], height=0.38, color=VIOLET, label="chain level")
        ax.barh(yy - 0.2, e.reindex(order)["effect"], height=0.38, color=GREY, label="entry level")
        for y_, v in zip(yy, order):
            q = c.loc[v, "q_value"]
            ax.text(c.loc[v, "effect"] + 0.005, y_ + 0.2, "*" if q < 0.05 else "", va="center", fontsize=9)
        ax.set_yticks(yy, [f"{v} ({c.loc[v, 'effect_type']})" for v in order], fontsize=8)
        ax.set(xlabel="effect size (Cramer's V or epsilon^2); * BH q < 0.05 (chain level)",
               title=f"Metadata vs {tg}")
        ax.legend(fontsize=8, loc="lower right")
    ax = fig.add_subplot(gs[2])
    ct = pd.crosstab(t.space_group8, t.cluster3)[["cl7", "cl4", "minor"]]
    ct = ct.loc[ct.sum(1).sort_values(ascending=False).index]
    frac = ct.div(ct.sum(1), axis=0)
    left = np.zeros(len(ct))
    for col, c in (("cl7", BLUE), ("cl4", ORANGE), ("minor", "#1baf7a")):
        ax.barh(np.arange(len(ct)), frac[col], left=left, color=c, label=col, height=0.7)
        left += frac[col].values
    ax.set_yticks(np.arange(len(ct)), [f"{s} (n={n})" for s, n in zip(ct.index, ct.sum(1))], fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="fraction of chains", title="Space group x cluster (chains)")
    ax.legend(fontsize=8, loc="lower right")
    top = df[(df.target == "structural cluster") & (df.level == "entry")].sort_values("effect", ascending=False).head(3)
    to = df[(df.target == "consensus outlier") & (df.level == "chain")]
    fig.suptitle(f"Metadata: crystal form, peptide and construct track the cluster (entry level: "
                 f"{', '.join(f'{r.variable} {r.effect:.2f}' for r in top.itertuples())}); "
                 f"{int((to.q_value < 0.05).sum())} of {len(to)} variables differ for consensus outliers",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y17_metadata_effects", section="E. Experimental metadata", task="Y17",
             shows="effect size of every mmCIF metadata variable against cluster and consensus-outlier status (chain and entry level); space group x cluster",
             takeaway=(f"Strongest cluster associations at entry level: {', '.join(r.variable for r in top.itertuples())}; "
                       f"pH and precipitant class are weaker; all entries are X-ray structures."),
             source="results/ext/metadata_tests.csv", script="scripts/ext_metadata_tests.py")
    write_numbers("Y17", dict(n_tests=len(df), n_sig=int((df.q_value < 0.05).sum()),
                              top_entry_cluster=top[["variable", "effect", "q_value"]].round(5).to_dict("records"),
                              exp_method_values=t.exp_method.unique().tolist()))
    print(df.sort_values(["target", "level", "effect"], ascending=[True, True, False]).round(4).to_string())


if __name__ == "__main__":
    main()
