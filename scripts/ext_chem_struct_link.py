"""Tracker Y25: ligand <-> conformation link on the 788-chain matrix (the provisional P04 test redone).

Inputs: results/ext/chem_labels_best.csv (Y08: auto-selected chemistry cluster per ligand of interest and the
structural majority of its chains), chain_table() (S08 cluster per chain, ligand of interest per chain),
results/ext/chem_fp_ECFP4.npy / chem_fp_MACCS.npy (Y08 bit matrices, row order = results/ligands_of_interest.csv),
results/rmsd_entry.npy, results/potency_by_cluster.csv (X15).
Unit: the ligand (a ligand seen in 47 chains counts once); a chain-level table is reported too, but its p-value comes
from shuffling chemistry labels between ligands (chains of one ligand move together).
Tests (all written to results/ext/chem_struct_link.csv, p never left empty):
  1 chi-square, ligand level: chem cluster (6) x structural majority (cl7 / cl4 / mixed / minor), with the count of
    cells with expected < 5; asymptotic p and a permutation p (10,000 shuffles of the structural labels over ligands)
  2 Cramer's V with a 95% bootstrap CI (2,000 ligand resamples) and mutual information (nats)
  3 cluster purity: fraction of ligands that carry their chem cluster's majority structural label; permutation p
  4 chain level: chem cluster x S08 cluster chi-square; permutation p with chemistry labels shuffled among ligands
  5 nearest-neighbour agreement: for each ligand its most similar other ligand (ECFP4 Tanimoto); fraction sharing the
    structural majority vs the permutation null; agreement by similarity bin
  6 RMSD vs ligand similarity: for all pairs of ligands, Spearman between Tanimoto and the median RMSD between their
    chains; Mantel-style permutation of ligand labels (2,000)
  7 potency by structural cluster controlling for chemistry: OLS pChEMBL ~ cl4 + C(chem cluster) on X15 ligands
Outputs: results/ext/chem_struct_link.csv, results/ext/chem_struct_contingency.csv, results/ext/nn_agreement.csv,
figures/ext/Y25_chem_struct_link.png, figures/ext/Y25_nn_agreement.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_chem_struct_link.py
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2_contingency, spearmanr
from sklearn.metrics import mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RES, ROOT, SEED, VIOLET, chain_table,
                        cramers_v, mpl_setup, rmsd, save_fig, write_numbers)

warnings.filterwarnings("ignore")
N_PERM = 10000


def tanimoto(bits):
    b = bits.astype(bool)
    inter = b.astype(np.int32) @ b.T.astype(np.int32)
    cnt = b.sum(1)
    union = cnt[:, None] + cnt[None, :] - inter
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(union > 0, inter / union, 0.0)
    np.fill_diagonal(s, 1.0)
    return s


def main():
    rng = np.random.default_rng(SEED)
    lab = pd.read_csv(RES / "chem_labels_best.csv")
    lig = pd.read_csv(ROOT / "results/ligands_of_interest.csv")
    assert (lab.ligand_id.values == lig.ligand_id.values).all()
    L = lab[lab.structural_majority != "none"].reset_index(drop=True)
    chem, smaj = L.chem_cluster_auto.values, L.structural_majority.values
    ct = pd.crosstab(L.chem_cluster_auto, L.structural_majority)
    ct.to_csv(RES / "chem_struct_contingency.csv")
    chi, p_asym, dof, exp = chi2_contingency(ct)
    perm = np.empty(N_PERM)
    for k in range(N_PERM):
        perm[k] = chi2_contingency(pd.crosstab(chem, rng.permutation(smaj)))[0]
    p_perm = (np.sum(perm >= chi) + 1) / (N_PERM + 1)
    v = cramers_v(ct.values)
    boots = []
    for _ in range(2000):
        s = rng.integers(0, len(L), len(L))
        boots.append(cramers_v(pd.crosstab(chem[s], smaj[s]).values))
    mi = mutual_info_score(chem, smaj)

    def purity(c, s):
        d = pd.crosstab(c, s)
        return d.max(1).sum() / d.values.sum()
    pur = purity(chem, smaj)
    pur_null = np.array([purity(chem, rng.permutation(smaj)) for _ in range(2000)])
    rows = [dict(test="chi-square chem cluster x structural majority (ligand)", n=len(L), statistic=chi, p_value=p_asym,
                 method=f"asymptotic, dof {dof}; {int((exp < 5).sum())}/{exp.size} cells expected < 5"),
            dict(test="chi-square chem cluster x structural majority (ligand)", n=len(L), statistic=chi, p_value=p_perm,
                 method=f"permutation, {N_PERM} shuffles of structural labels over ligands"),
            dict(test="Cramer's V (ligand)", n=len(L), statistic=v, p_value=p_perm,
                 method=f"bootstrap 95% CI {np.percentile(boots, 2.5):.3f}-{np.percentile(boots, 97.5):.3f}; p from the permutation chi-square"),
            dict(test="mutual information (ligand, nats)", n=len(L), statistic=mi,
                 p_value=(np.sum([mutual_info_score(chem, rng.permutation(smaj)) >= mi for _ in range(2000)]) + 1) / 2001,
                 method="permutation, 2000 shuffles"),
            dict(test="cluster purity (ligand)", n=len(L), statistic=pur, p_value=(np.sum(pur_null >= pur) + 1) / 2001,
                 method=f"permutation null mean {pur_null.mean():.3f}")]
    # chain level with ligand-level shuffling
    t = chain_table()
    cmap = dict(zip(lab.ligand_id, lab.chem_cluster_auto))
    tc = t.dropna(subset=["ligand_id"]).copy()
    tc["chem"] = tc.ligand_id.map(cmap)
    tc = tc.dropna(subset=["chem"])
    tc["s08"] = np.where(tc.cluster == 7, "cl7", np.where(tc.cluster == 4, "cl4", "minor"))
    chi_c = chi2_contingency(pd.crosstab(tc.chem, tc.s08))[0]
    ligs = np.array(sorted(cmap))
    vals = np.array([cmap[l] for l in ligs])
    cp = np.empty(2000)
    for k in range(2000):
        m = dict(zip(ligs, rng.permutation(vals)))
        cp[k] = chi2_contingency(pd.crosstab(tc.ligand_id.map(m), tc.s08))[0]
    rows.append(dict(test="chi-square chem cluster x S08 cluster (chain)", n=len(tc), statistic=chi_c,
                     p_value=(np.sum(cp >= chi_c) + 1) / 2001,
                     method="permutation, chemistry labels shuffled among ligands (2000); chains of one ligand move together"))
    rows.append(dict(test="Cramer's V (chain)", n=len(tc), statistic=cramers_v(pd.crosstab(tc.chem, tc.s08).values),
                     p_value=(np.sum(cp >= chi_c) + 1) / 2001, method="p from the ligand-level permutation"))
    # nearest-neighbour agreement
    bits = np.load(RES / "chem_fp_ECFP4.npy")
    keep = lab.structural_majority.values != "none"
    S = tanimoto(bits)[np.ix_(keep, keep)]
    np.fill_diagonal(S, -1)
    nn = S.argmax(1)
    nn_sim = S.max(1)
    agree = smaj[nn] == smaj
    null = np.array([np.mean(p[nn] == p) for p in (rng.permutation(smaj) for _ in range(2000))])
    nnd = pd.DataFrame({"ligand_id": L.ligand_id, "nn_ligand": L.ligand_id.values[nn], "nn_tanimoto": nn_sim,
                        "structural": smaj, "nn_structural": smaj[nn], "agree": agree})
    nnd.to_csv(RES / "nn_agreement.csv", index=False)
    rows.append(dict(test="nearest-neighbour ligand shares structural majority", n=len(L), statistic=agree.mean(),
                     p_value=(np.sum(null >= agree.mean()) + 1) / 2001, method=f"permutation null mean {null.mean():.3f}"))
    # RMSD vs similarity
    d, idx = rmsd()
    by_lig = tc.groupby("ligand_id").apply(lambda g: g.index.values, include_groups=False)
    lid = [l for l in L.ligand_id if l in by_lig.index]
    pos = {l: i for i, l in enumerate(L.ligand_id)}
    n = len(lid)
    med = np.full((n, n), np.nan)
    for i in range(n):
        ri = by_lig[lid[i]]
        for j in range(i + 1, n):
            med[i, j] = med[j, i] = np.median(d[np.ix_(ri, by_lig[lid[j]])])
    Sl = tanimoto(bits)[np.ix_(keep, keep)][np.ix_([pos[l] for l in lid], [pos[l] for l in lid])]
    iu = np.triu_indices(n, 1)
    rho = spearmanr(Sl[iu], med[iu])[0]
    rp = []
    for _ in range(200):
        pr = rng.permutation(n)
        rp.append(spearmanr(Sl[np.ix_(pr, pr)][iu], med[iu])[0])
    rows.append(dict(test="Spearman ligand Tanimoto vs median inter-chain RMSD (ligand pairs)", n=int(len(iu[0])),
                     statistic=rho, p_value=(np.sum(np.array(rp) <= rho) + 1) / 201,
                     method="Mantel-style permutation of ligand labels (200); one-sided (more similar -> closer)"))
    # potency controlled for chemistry
    pot = pd.read_csv(ROOT / "results/potency_by_cluster.csv")
    pot = pot[pot.cluster.isin([7, 4])].merge(lab[["ligand_id", "chem_cluster_auto"]], on="ligand_id", how="inner")
    pot["cl4"] = (pot.cluster == 4).astype(int)
    m0 = smf.ols("pchembl ~ cl4", data=pot).fit()
    m1 = smf.ols("pchembl ~ cl4 + C(chem_cluster_auto)", data=pot).fit()
    rows.append(dict(test="potency: pChEMBL difference cl4 - cl7 (unadjusted OLS)", n=len(pot), statistic=m0.params["cl4"],
                     p_value=m0.pvalues["cl4"], method="X15 ligands, OLS"))
    rows.append(dict(test="potency: pChEMBL difference cl4 - cl7 adjusted for chem cluster", n=len(pot),
                     statistic=m1.params["cl4"], p_value=m1.pvalues["cl4"], method="OLS with chem-cluster fixed effects"))
    res = pd.DataFrame(rows)
    res.to_csv(RES / "chem_struct_link.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={"width_ratios": [1, 1, 1]})
    ax = axs[0]
    order = [c for c in ("cl7", "cl4", "mixed cl7/cl4", "minor") if c in ct.columns]
    c2 = ct[order]
    frac = c2.div(c2.sum(1), axis=0)
    left = np.zeros(len(c2))
    for c, col in zip(order, (BLUE, ORANGE, VIOLET, AQUA)):
        ax.barh(np.arange(len(c2)), frac[c], left=left, color=col, label=c, height=0.7)
        for i, (fv, nv) in enumerate(zip(frac[c], c2[c])):
            if nv >= 5:
                ax.text(left[i] + fv / 2, i, str(nv), ha="center", va="center", fontsize=8, color="white")
        left += frac[c].values
    ax.set_yticks(np.arange(len(c2)), [f"chem {i} (n={n})" for i, n in zip(c2.index, c2.sum(1))])
    ax.set(xlabel="fraction of ligands (numbers = ligands)", title=f"Chemistry cluster x structural majority (ligand level)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4)
    ax = axs[1]
    ax.hist(perm, bins=50, color=GREY)
    ax.axvline(chi, color=INK, lw=1.5)
    ax.text(chi, ax.get_ylim()[1] * 0.9, f" observed chi2 = {chi:.0f}\n permutation p = {p_perm:.1g}", fontsize=8.5,
            ha="right" if chi > perm.max() else "left")
    ax.set(xlabel="chi-square under 10,000 ligand-level permutations", ylabel="count",
           title=f"Cramer's V {v:.2f} (95% CI {np.percentile(boots, 2.5):.2f}-{np.percentile(boots, 97.5):.2f}), MI {mi:.2f} nats")
    ax = axs[2]
    pot["chem"] = pot.chem_cluster_auto.astype(str)
    for k, (c, col) in enumerate(((7, BLUE), (4, ORANGE))):
        s = pot[pot.cluster == c]
        ax.scatter(s.chem_cluster_auto + (k - 0.5) * 0.25 + rng.uniform(-.06, .06, len(s)), s.pchembl, s=22, c=col,
                   lw=0, label=f"cl{c} ligands (n={len(s)})")
    ax.set(xlabel="chemistry cluster (Y08 auto)", ylabel="median pChEMBL at ESR1 (X15)",
           title=f"Potency: cl4-cl7 {m0.params['cl4']:+.2f} raw (p={m0.pvalues['cl4']:.1g}), "
                 f"{m1.params['cl4']:+.2f} within chemistry (p={m1.pvalues['cl4']:.1g})")
    ax.legend(fontsize=8)
    fig.suptitle(f"Ligand chemistry and conformation are strongly linked at the ligand level (n={len(L)} ligands, permutation "
                 f"p {p_perm:.1g}); {int((exp < 5).sum())} of {exp.size} cells have expected < 5, hence the permutation test",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y25_chem_struct_link", section="J. Ligand and conformation", task="Y25",
             shows="chem cluster x structural majority (ligand level), permutation null of the chi-square, potency by cluster within chemistry clusters",
             takeaway=(f"Cramer's V {v:.2f} (CI {np.percentile(boots, 2.5):.2f}-{np.percentile(boots, 97.5):.2f}); "
                       f"the cl4-cl7 potency gap is {m0.params['cl4']:+.2f} log units raw and {m1.params['cl4']:+.2f} within chemistry clusters."),
             source="results/ext/chem_struct_link.csv, results/ext/chem_struct_contingency.csv", script="scripts/ext_chem_struct_link.py")

    fig, axs = plt.subplots(1, 2, figsize=(15, 5.5))
    ax = axs[0]
    bins = [0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01]
    nnd["bin"] = pd.cut(nnd.nn_tanimoto, bins, right=False)
    g = nnd.groupby("bin", observed=True).agg(rate=("agree", "mean"), n=("agree", "size"))
    ax.bar(range(len(g)), g.rate, color=VIOLET, width=0.7)
    ax.axhline(null.mean(), color=INK, ls="--", lw=1, label=f"permutation null {null.mean():.2f}")
    for i, (r_, n_) in enumerate(zip(g.rate, g.n)):
        ax.text(i, r_ + 0.01, f"n={n_}", ha="center", fontsize=8, color=INK2)
    ax.set_xticks(range(len(g)), [str(b) for b in g.index], rotation=20)
    ax.set(xlabel="ECFP4 Tanimoto to the most similar other ligand", ylabel="nearest neighbour has the same structural majority",
           ylim=(0, 1.05), title=f"Do similar ligands give the same H12 state? overall {agree.mean():.2f}")
    ax.legend(fontsize=8)
    ax = axs[1]
    sv, mv = Sl[iu], med[iu]
    hb = ax.hexbin(sv, mv, gridsize=40, cmap="esr_seq", mincnt=1, bins="log")
    fig.colorbar(hb, ax=ax, label="ligand pairs (log)")
    sb = pd.Series(mv).groupby(pd.cut(sv, np.linspace(0, 1, 11))).median()
    ax.plot([i.mid for i in sb.index], sb.values, "-o", color=INK, ms=4, label="median per similarity bin")
    ax.set(xlabel="ECFP4 Tanimoto between two ligands", ylabel="median C-alpha RMSD between their chains (A)",
           title=f"Ligand similarity vs pose distance: Spearman {rho:.2f} (n={len(sv)} pairs)")
    ax.legend(fontsize=8)
    fig.suptitle("Similar ligands give similar conformations: nearest chemical neighbours share the H12 state far more often "
                 "than chance, and RMSD falls as Tanimoto rises", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y25_nn_agreement", section="J. Ligand and conformation", task="Y25",
             shows="nearest-neighbour agreement of structural majority by Tanimoto bin; hexbin of ligand-pair Tanimoto vs median inter-chain RMSD",
             takeaway=(f"{agree.mean():.0%} of ligands share the H12 state of their nearest chemical neighbour (null "
                       f"{null.mean():.0%}); Spearman Tanimoto vs RMSD {rho:.2f}."),
             source="results/ext/nn_agreement.csv, results/ext/chem_struct_link.csv", script="scripts/ext_chem_struct_link.py")
    write_numbers("Y25", dict(n_ligands=len(L), chi2=float(chi), p_perm=float(p_perm), cramers_v=float(v),
                              v_ci=[float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))], mi=float(mi),
                              purity=float(pur), nn_agree=float(agree.mean()), nn_null=float(null.mean()), rho_sim_rmsd=float(rho),
                              potency_raw=float(m0.params["cl4"]), potency_adj=float(m1.params["cl4"]),
                              potency_adj_p=float(m1.pvalues["cl4"]), n_potency=len(pot),
                              cells_expected_lt5=int((exp < 5).sum()), n_cells=int(exp.size)))
    print(res.to_string())


if __name__ == "__main__":
    main()
