"""Tracker Y16: is a structure an outlier because of quality, or because it is a different conformation?

Inputs: results/ext/outlier_flags.csv (Y15), chain_table(), results/ext/entry_metadata.csv and chain_metadata.csv (Y02),
results/ext/structural_features.csv (Y11).
Groups compared with the rest: consensus outliers (>= 3 methods), each single detector (knn, medoid, lof, iforest;
mahal is too broad, see Y15), and the 14 minor-cluster chains. Two references: (a) all other chains; (b) the other
chains of the same main cluster (so the cl7/cl4 quality difference, Y14, cannot masquerade as an outlier effect;
groups that are only minor-cluster chains have no reference b).
Variables: resolution, R-free, R-work, entry mean B, chain mean C-alpha B, completeness, missing LBD C-alpha
(246 - observed), H12 not complete, tier not A, high B-group (quality); peptide present, Y537S, L536S, space group
P21, crystallisation pH, PEG, ASU chains, release year (experimental conditions).
Tests: Mann-Whitney U (rank-biserial r, + = larger in the outlier group) and Fisher exact (log odds ratio, Haldane).
BH over every test in the file.
Per consensus outlier, a verdict: 'quality' if any of resolution, R-free, mean C-alpha B or missing C-alpha lies
above the 95th percentile of the non-flagged chains (or H12 is not modelled); otherwise 'conformation'.
Outputs: results/ext/outlier_tests.csv, results/ext/outlier_verdicts.csv, figures/ext/Y16_outlier_quality.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_outlier_quality.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK, INK2, ORANGE, RED, RES, VIOLET, bh, chain_table, mpl_setup, save_fig,
                        write_numbers)

warnings.filterwarnings("ignore")
QUAL = ["resolution_A", "r_free", "r_work", "mean_b_A2", "mean_b_chain", "completeness_pct", "missing_lbd_ca",
        "h12_not_complete", "tier_not_A", "b_group_high"]
COND = ["peptide_present", "mut_Y537S", "mut_L536S", "sg_P21", "ph", "precip_PEG", "chains_asu", "release_year"]
BINARY = {"h12_not_complete", "tier_not_A", "b_group_high", "peptide_present", "mut_Y537S", "mut_L536S", "sg_P21",
          "precip_PEG"}


def table():
    t = chain_table()
    o = pd.read_csv(RES / "outlier_flags.csv").drop(columns=["entry_id", "cluster"])
    cm = pd.read_csv(RES / "chain_metadata.csv")[["chain_key", "mutation_status", "peptide_present"]]
    em = pd.read_csv(RES / "entry_metadata.csv")[["entry_id", "r_work", "completeness_pct", "ph", "additive_class",
                                                   "chains_asu", "release_year"]]
    sf = pd.read_csv(RES / "structural_features.csv")[["chain_key", "lbd_ca_observed", "h12_shift_A"]]
    t = t.merge(o, on="chain_key").merge(cm, on="chain_key").merge(em, on="entry_id", how="left").merge(sf, on="chain_key")
    t["missing_lbd_ca"] = 246 - t.lbd_ca_observed
    t["h12_not_complete"] = t.h12_label != "complete"
    t["tier_not_A"] = t.tier != "A"
    t["b_group_high"] = t.b_group == "high"
    t["mut_Y537S"] = t.mutation_status == "Y537S"
    t["mut_L536S"] = t.mutation_status == "L536S"
    t["sg_P21"] = t.space_group == "P 1 21 1"
    t["precip_PEG"] = t.additive_class == "PEG"
    t["minor_cluster"] = ~t.cluster.isin((7, 4))
    return t


def compare(t, gname, gmask, ref, refmask):
    rows = []
    for v in QUAL + COND:
        a, b = t.loc[gmask, v], t.loc[refmask & ~gmask, v]
        if v in BINARY:
            a, b = a.astype(bool), b.astype(bool)
            tab = [[a.sum(), (~a).sum()], [b.sum(), (~b).sum()]]
            _, p = fisher_exact(tab)
            eff = np.log((tab[0][0] + .5) * (tab[1][1] + .5) / ((tab[0][1] + .5) * (tab[1][0] + .5)))
            rows.append(dict(group=gname, reference=ref, variable=v, kind="quality" if v in QUAL else "conditions",
                             test="Fisher", n_group=int(len(a)), n_ref=int(len(b)), value_group=a.mean(),
                             value_ref=b.mean(), effect=eff, effect_type="log OR", p_value=p))
        else:
            a, b = a.astype(float).dropna(), b.astype(float).dropna()
            if len(a) < 3 or len(b) < 3:
                continue
            u, p = mannwhitneyu(a, b)
            rows.append(dict(group=gname, reference=ref, variable=v, kind="quality" if v in QUAL else "conditions",
                             test="Mann-Whitney", n_group=int(len(a)), n_ref=int(len(b)), value_group=a.median(),
                             value_ref=b.median(), effect=2 * u / (len(a) * len(b)) - 1, effect_type="rank-biserial",
                             p_value=p))
    return rows


def main():
    t = table()
    groups = {"consensus (>=3)": t.consensus_outlier, "knn": t.knn_flag, "medoid": t.medoid_flag, "lof": t.lof_flag,
              "iforest": t.iforest_flag, "minor clusters": t.minor_cluster}
    rows = []
    allm = pd.Series(True, index=t.index)
    for g, m in groups.items():
        m = m.astype(bool)
        rows += compare(t, g, m, "all other chains", allm)
        if g != "minor clusters":
            mm = m & ~t.minor_cluster
            if mm.sum() >= 3:
                same = t.cluster.isin(t.loc[mm, "cluster"].unique()) & ~t.minor_cluster
                rows += compare(t, g + " (in main cl.)", mm, "same main cluster", same)
    df = pd.DataFrame(rows)
    df["q_value"] = bh(df.p_value)
    df.to_csv(RES / "outlier_tests.csv", index=False)

    clean = t.n_methods == 0
    thr = {v: t.loc[clean, v].quantile(0.95) for v in ("resolution_A", "r_free", "mean_b_chain", "missing_lbd_ca")}
    ver = []
    for r in t[t.consensus_outlier].itertuples():
        reasons = [f"{v} {getattr(r, v):.3g} > p95 {thr[v]:.3g}" for v in thr if getattr(r, v) > thr[v]]
        if r.h12_label == "not_modelled":
            reasons.append("H12 not modelled")
        ver.append(dict(chain_key=r.chain_key, entry_id=r.entry_id, cluster=r.cluster, n_methods=r.n_methods,
                        resolution_A=r.resolution_A, r_free=r.r_free, mean_b_chain=r.mean_b_chain,
                        missing_lbd_ca=r.missing_lbd_ca, h12_label=r.h12_label, h12_shift_A=r.h12_shift_A,
                        mutation_status=r.mutation_status, peptide=r.peptide_present, space_group=r.space_group,
                        verdict="quality" if reasons else "conformation", reasons="; ".join(reasons) or "none"))
    ver = pd.DataFrame(ver)
    ver.to_csv(RES / "outlier_verdicts.csv", index=False)

    plt = mpl_setup()
    fig = plt.figure(figsize=(19, 9.5))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15])
    rng = np.random.default_rng(42)
    for k, (v, lab) in enumerate((("resolution_A", "resolution (A)"), ("r_free", "R-free"),
                                  ("mean_b_chain", "chain mean C-alpha B (A^2)"), ("missing_lbd_ca", "missing LBD C-alpha"))):
        ax = fig.add_subplot(gs[0, k])
        cats = [("not flagged", t.n_methods == 0, GREY), ("1-2 methods", t.n_methods.between(1, 2), VIOLET),
                ("consensus >= 3", t.consensus_outlier, RED)]
        for i, (nm, m, col) in enumerate(cats):
            vals = t.loc[m, v].astype(float).dropna().values
            parts = ax.violinplot([vals], positions=[i], widths=0.8, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(col); pc.set_alpha(0.3); pc.set_edgecolor("none")
            ax.scatter(i + rng.uniform(-.15, .15, len(vals)), vals, s=4 if i < 2 else 14, c=col, lw=0, alpha=.7)
            ax.hlines(np.median(vals), i - .3, i + .3, color=INK, lw=1.4)
        r = df[(df.group == "consensus (>=3)") & (df.reference == "all other chains") & (df.variable == v)].iloc[0]
        ax.set_xticks([0, 1, 2], [f"{c[0]}\nn={int(c[1].sum())}" for c in cats], fontsize=8)
        ax.set(ylabel=lab, title=f"{lab}: consensus vs rest q={r.q_value:.2g}")
    ax = fig.add_subplot(gs[1, :3])
    sub = df[df.reference == "all other chains"]
    piv = sub.pivot(index="group", columns="variable", values="effect")[QUAL + COND]
    qv = sub.pivot(index="group", columns="variable", values="q_value")[QUAL + COND]
    eff_scaled = piv.copy()
    for c in piv.columns:
        if c in BINARY:
            eff_scaled[c] = np.tanh(piv[c] / 3)  # log OR squashed into [-1, 1] for one colour scale
    im = ax.imshow(eff_scaled.values, cmap="esr_div", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(piv.columns)), piv.columns, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(piv.index)), piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            q = qv.values[i, j]
            if np.isfinite(q):
                ax.text(j, i, f"{piv.values[i, j]:+.2f}" + ("*" if q < 0.05 else ""), ha="center", va="center",
                        fontsize=6.5, color=INK)
    ax.axvline(len(QUAL) - 0.5, color=INK, lw=1)
    ax.text(len(QUAL) / 2 - .5, -0.62, "structure quality", ha="center", fontsize=9, fontweight="bold")
    ax.text(len(QUAL) + len(COND) / 2 - .5, -0.62, "experimental conditions", ha="center", fontsize=9, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.02, label="effect (rank-biserial; log OR squashed by tanh(x/3)); * q < 0.05")
    ax.set_title("Outlier groups vs all other chains: effect sizes", pad=26)
    ax = fig.add_subplot(gs[1, 3])
    vc = ver.verdict.value_counts()
    ax.bar(vc.index, vc.values, color=[RED if i == "quality" else VIOLET for i in vc.index], width=0.6)
    for i, (k_, v_) in enumerate(vc.items()):
        ax.text(i, v_, str(v_), ha="center", va="bottom", fontsize=9)
    ax.set(ylabel="consensus outlier chains", title="Verdict per consensus outlier")
    nq = int((ver.verdict == "quality").sum())
    fig.suptitle(f"Outliers vs quality: {nq} of {len(ver)} consensus outliers have a quality reason "
                 f"(resolution / R-free / B / missing residues beyond the 95th percentile, or no H12); "
                 f"{len(ver) - nq} look like genuine conformations", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    fig.subplots_adjust(top=0.92)
    save_fig(fig, "Y16_outlier_quality", section="D. Outliers and quality", task="Y16",
             shows="quality metrics of consensus outliers vs the rest; effect-size heatmap for every outlier group and variable; verdict counts",
             takeaway=(f"{nq}/{len(ver)} consensus outliers are explained by quality (or a missing H12), "
                       f"{len(ver) - nq} are well-determined structures in an unusual conformation; see outlier_verdicts.csv."),
             source="results/ext/outlier_tests.csv, results/ext/outlier_verdicts.csv", script="scripts/ext_outlier_quality.py")
    write_numbers("Y16", dict(n_tests=len(df), n_sig=int((df.q_value < 0.05).sum()), n_consensus=len(ver),
                              n_quality=nq, n_conformation=len(ver) - nq, thresholds_p95=thr,
                              verdicts=ver[["chain_key", "verdict", "reasons"]].to_dict("records")))
    print(ver[["chain_key", "cluster", "verdict", "reasons"]].to_string())
    print(df[(df.group == "consensus (>=3)") & (df.q_value < 0.1)][["reference", "variable", "value_group", "value_ref", "effect", "q_value"]].to_string())


if __name__ == "__main__":
    main()
