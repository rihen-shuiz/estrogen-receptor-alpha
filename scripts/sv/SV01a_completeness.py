"""SV01a (RUN_SV1028): H12 completeness, first pass - distribution, within-entry concordance, Q/X/B block models.

Input (tracker card SV00 out): results/sv/SV00/chain_table.csv, results/sv/SV00/observed_mask.npy (+ row/col csv).
Chains flagged `remediated` are excluded (SV00: 0 chains). The S08 cluster is NOT a predictor; models run within cl7
and within cl4. H12 window 538-547 (run file section 2).

Rules fixed before any model was fitted (written here so that they cannot follow a result):
- Outcome (primary, [DEFAULT]): H12 complete = all 10 C-alpha of 538-547 observed in observed_mask.npy.
  Sensitivity outcome: fraction of 538-547 observed (fractional logit: binomial GLM, logit link).
- Predictors (SV01 step 3):
  Q  resolution_A, r_free, mean B = mean_b_entry_A2 (the entry mean B of Q05/K05, the "mean B" used in the QC tiers
     and in Y14; the chain C-alpha mean B is not used because it averages over the very residues whose absence is
     the outcome).
  X  ph, growth_temp_K, ambient_temp_K, space group P21 (space_group == 'P 1 21 1') vs other, esr1_chains_asu,
     dimer source (dimer_partner_source == 'symmetry' vs 'asu').
  B  peptide (chain-level peptide_bound), partner has peptide (dimer_partner_peptide_bound), Y537S, L536S, D538G,
     other mutations (any construct mutation other than these three, from construct_mutations), basic amine,
     carboxylic acid, ligand heavy atoms. Chains without a ligand of interest: amine 0, acid 0, heavy atoms 0.
- Missing values (pH 36 %, growth T 7 %, ambient T 3 % of chains): median imputation within the cluster plus a
  0/1 '<var> not reported' indicator in block X (created only when the cluster has a missing value). All models of a
  cluster then use the same rows, so McFadden gains are comparable.
- Continuous predictors are z-scored with the chain-level mean/SD of the cluster; odds ratios are per 1 SD (the SD is
  saved). The entry level reuses the chain-level scaling.
- Levels: chain level = every chain of the cluster, logit MLE (statsmodels Logit; fractional logit QMLE for the
  fraction outcome) with entry-clustered (sandwich) errors, i.e. GEE with an independence working correlation. Entry level = one chain per entry (the entry's first matrix chain in that
  cluster; the Y14 convention), ordinary MLE.
- Drop rule (SV01 step 3), per cluster x level x outcome: a predictor that is constant, or that separates the
  outcome (binary predictor: one level has every outcome at 1 or every outcome at 0; continuous predictor: the
  ranges for complete and incomplete H12 do not overlap) is dropped there and listed. A predictor that is an exact
  linear combination of earlier predictors (block order Q, X, B) is dropped as collinear and listed.
  If the full model still separates (statsmodels perfect-separation warning, no IRLS convergence, or any
  |log OR| > 15), the predictor with the largest |log OR| is dropped as 'joint separation' and the model refitted,
  until it is stable. The surviving predictor set is used for every sub-model of that cluster x level x outcome.
- McFadden R2 = 1 - LL/LL0 (LL from the fitted probabilities; for the fraction outcome the Bernoulli quasi-LL).
  Gain entered first = R2(block alone); entered last = R2(Q+X+B) - R2(other blocks). Blocks Q, X, B and Q+X.
- Tests: chain level, cluster-robust Wald test of the block in the full model; entry level, LR test. Odds ratios:
  Wald 95 % CI (cluster-robust at chain level). One BH family over every p-value of this step.
- [DEFAULT] verdict: "artifact" if the Q+X gain exceeds the B gain in both entry orders and in both clusters,
  "biology" if the reverse holds, otherwise "mixed". Applied to the primary outcome at chain level and at entry level;
  the step verdict is "artifact"/"biology" only if the condition holds at both levels (rule 1.4: every test per chain
  and per entry), otherwise "mixed". The fraction outcome gets the same verdict as a sensitivity result.

Outputs: results/sv/SV01a/ (completeness_profile.csv, h12_status_counts.csv, partial_ends.csv, chain_outcomes.csv,
concordance.csv, concordance_summary.csv, model_data.csv, predictor_scaling.csv, model_blocks.csv,
odds_ratios.csv, verdict_table.csv, verdict.md, numbers.json, env.txt).
Run: .venv/bin/python -I scripts/sv/SV01a_completeness.py
"""
from __future__ import annotations

import json
import platform
import sys
import warnings
from datetime import datetime, timezone, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
SV00 = ROOT / "results/sv/SV00"
OUT = ROOT / "results/sv/SV01a"
OUT.mkdir(parents=True, exist_ok=True)
H12 = (538, 547)
MAIN = (7, 4)
COEF_LIMIT = 15.0

Q = ["resolution_A", "r_free", "mean_b"]
X = ["ph", "growth_temp_K", "ambient_temp_K", "p21", "esr1_chains_asu", "dimer_symmetry"]
B = ["peptide", "partner_peptide", "mut_Y537S", "mut_L536S", "mut_D538G", "mut_other", "basic_amine",
     "carboxylic_acid", "ligand_heavy_atoms"]
CONTINUOUS = {"resolution_A", "r_free", "mean_b", "ph", "growth_temp_K", "ambient_temp_K", "esr1_chains_asu",
              "ligand_heavy_atoms"}
IMPUTE = ["ph", "growth_temp_K", "ambient_temp_K"]
LABEL = {
    "resolution_A": "resolution (A)", "r_free": "R-free", "mean_b": "mean B, entry (A^2)", "ph": "growth pH",
    "growth_temp_K": "growth temperature (K)", "ambient_temp_K": "data-collection temperature (K)",
    "p21": "space group P21", "esr1_chains_asu": "ESR1 chains in ASU", "dimer_symmetry": "dimer partner by symmetry",
    "ph_missing": "pH not reported", "growth_temp_K_missing": "growth T not reported",
    "ambient_temp_K_missing": "ambient T not reported", "peptide": "peptide bound", "partner_peptide":
    "partner has peptide", "mut_Y537S": "Y537S", "mut_L536S": "L536S", "mut_D538G": "D538G",
    "mut_other": "other mutations", "basic_amine": "basic amine", "carboxylic_acid": "carboxylic acid",
    "ligand_heavy_atoms": "ligand heavy atoms"}


def block_of(p):
    if p in Q:
        return "Q"
    if p in X or p.endswith("_missing"):
        return "X"
    return "B"


# ------------------------------------------------------------------------------------------------ inputs
def load():
    t = pd.read_csv(SV00 / "chain_table.csv", low_memory=False)
    rows = pd.read_csv(SV00 / "observed_mask_rows.csv")
    cols = pd.read_csv(SV00 / "observed_mask_cols.csv")
    mask = np.load(SV00 / "observed_mask.npy")
    assert mask.shape == (len(rows), len(cols)) and (rows.chain_key.values == t.chain_key.values).all()
    pos = cols.esr1_position.values
    keep = ~t.remediated.astype(bool).values
    t, mask = t[keep].reset_index(drop=True), mask[keep]
    w = (pos >= H12[0]) & (pos <= H12[1])
    hw = mask[:, w]
    t["n_obs_window"] = hw.sum(1)
    t["frac_obs"] = hw.mean(1)
    t["complete"] = (t.n_obs_window == w.sum()).astype(int)
    t["h12_status"] = np.where(t.n_obs_window == w.sum(), "complete",
                               np.where(t.n_obs_window == 0, "not_modelled", "partial"))
    t["cluster_group"] = np.where(t.s08_cluster == 7, "cl7", np.where(t.s08_cluster == 4, "cl4", "minor"))
    return t, mask, pos


def code_predictors(t):
    d = pd.DataFrame(index=t.index)
    d["resolution_A"] = t.resolution_A.astype(float)
    d["r_free"] = t.r_free.astype(float)
    d["mean_b"] = t.mean_b_entry_A2.astype(float)
    d["ph"] = t.ph.astype(float)
    d["growth_temp_K"] = t.growth_temp_K.astype(float)
    d["ambient_temp_K"] = t.ambient_temp_K.astype(float)
    d["p21"] = (t.space_group == "P 1 21 1").astype(float)
    d["esr1_chains_asu"] = t.esr1_chains_asu.astype(float)
    d["dimer_symmetry"] = (t.dimer_partner_source == "symmetry").astype(float)
    d["peptide"] = t.peptide_bound.astype(str).eq("True").astype(float)
    d["partner_peptide"] = t.dimer_partner_peptide_bound.astype(str).eq("True").astype(float)
    muts = t.construct_mutations.fillna("none").astype(str).str.split(";")
    d["mut_Y537S"] = muts.apply(lambda m: "Y537S" in m).astype(float)
    d["mut_L536S"] = muts.apply(lambda m: "L536S" in m).astype(float)
    d["mut_D538G"] = muts.apply(lambda m: "D538G" in m).astype(float)
    d["mut_other"] = muts.apply(lambda m: any(x not in ("Y537S", "L536S", "D538G", "none", "") for x in m)).astype(float)
    d["basic_amine"] = t.ligand_basic_amine.astype(str).eq("True").astype(float)
    d["carboxylic_acid"] = t.ligand_carboxylic_acid.astype(str).eq("True").astype(float)
    d["ligand_heavy_atoms"] = t.ligand_heavy_atoms.astype(float).fillna(0.0)
    return d


# ------------------------------------------------------------------------------------------------ step 1: distribution
def profile(t, mask, pos):
    groups = [("all", "all", np.ones(len(t), bool))]
    for g in ("cl7", "cl4", "minor"):
        groups.append(("cluster", g, (t.cluster_group == g).values))
    for m in sorted(t.mutation_status.unique()):
        groups.append(("mutation", m, (t.mutation_status == m).values))
    for g in ("cl7", "cl4"):
        for m in sorted(t.mutation_status.unique()):
            s = ((t.cluster_group == g) & (t.mutation_status == m)).values
            if s.any():
                groups.append(("cluster_mutation", f"{g}|{m}", s))
    out = []
    for gt, g, s in groups:
        sub, ent = mask[s], t.entry_id.values[s]
        ue = pd.unique(ent)
        ent_frac = np.vstack([sub[ent == e].mean(0) for e in ue])
        for j, p in enumerate(pos):
            out.append(dict(group_type=gt, group=g, residue=int(p), n_chains=int(s.sum()), n_observed=int(sub[:, j].sum()),
                            frac_observed=sub[:, j].mean(), n_entries=len(ue), frac_observed_entry=ent_frac[:, j].mean()))
    return pd.DataFrame(out)


def status_counts(t):
    out = []
    specs = [("cluster", "cluster_group", ("cl7", "cl4", "minor")),
             ("mutation", "mutation_status", sorted(t.mutation_status.unique()))]
    for level in ("chain", "entry"):
        for gt, col, vals in specs + [("all", None, ("all",))]:
            for v in vals:
                s = t if col is None else t[t[col] == v]
                if level == "entry":
                    s = s.drop_duplicates("entry_id")
                vc = s.h12_status.value_counts()
                out.append(dict(level=level, group_type=gt, group=v, n=len(s), complete=int(vc.get("complete", 0)),
                                partial=int(vc.get("partial", 0)), not_modelled=int(vc.get("not_modelled", 0)),
                                frac_complete=vc.get("complete", 0) / len(s) if len(s) else np.nan))
        for g in ("cl7", "cl4"):
            for m in sorted(t.mutation_status.unique()):
                s = t[(t.cluster_group == g) & (t.mutation_status == m)]
                if level == "entry":
                    s = s.drop_duplicates("entry_id")
                if len(s):
                    vc = s.h12_status.value_counts()
                    out.append(dict(level=level, group_type="cluster_mutation", group=f"{g}|{m}", n=len(s),
                                    complete=int(vc.get("complete", 0)), partial=int(vc.get("partial", 0)),
                                    not_modelled=int(vc.get("not_modelled", 0)),
                                    frac_complete=vc.get("complete", 0) / len(s)))
    return pd.DataFrame(out)


def partial_ends(t, mask, pos):
    w = np.where((pos >= H12[0]) & (pos <= H12[1]))[0]
    out = []
    for i in np.where(t.h12_status.values == "partial")[0]:
        ob = mask[i, w]
        miss = [int(pos[j]) for j, o in zip(w, ob) if not o]
        n_n = int(np.argmax(ob)) if ob.any() else len(ob)
        n_c = int(np.argmax(ob[::-1])) if ob.any() else len(ob)
        internal = len(miss) - n_n - n_c
        side = ("both ends" if n_n and n_c else "N-side" if n_n else "C-side" if n_c else "internal only")
        obs_all = pos[mask[i]]
        out.append(dict(chain_key=t.chain_key[i], entry_id=t.entry_id[i], s08_cluster=int(t.s08_cluster[i]),
                        cluster_group=t.cluster_group[i], mutation_status=t.mutation_status[i],
                        n_obs_window=int(ob.sum()), missing_residues=";".join(map(str, miss)), n_missing_N_end=n_n,
                        n_missing_C_end=n_c, n_missing_internal=internal, side=side,
                        last_observed_residue=int(obs_all.max()) if len(obs_all) else np.nan))
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ step 2: concordance
def concordance(t):
    out = []
    for e, g in t.groupby("entry_id", sort=False):
        if len(g) < 2:
            continue
        lab = g.h12_status.value_counts()
        cl = g.s08_cluster.value_counts()
        same_cl = g.groupby("s08_cluster").h12_status.nunique()
        out.append(dict(entry_id=e, n_chains=len(g),
                        h12_labels=";".join(f"{k}:{v}" for k, v in sorted(lab.items())),
                        n_distinct_labels=len(lab), label_differs=len(lab) > 1,
                        complete_differs=g.complete.nunique() > 1,
                        clusters=";".join(f"cl{k}:{v}" for k, v in sorted(cl.items())),
                        n_distinct_clusters=len(cl), cluster_differs=len(cl) > 1,
                        label_differs_within_a_cluster=bool((same_cl > 1).any()),
                        has_two_chains_in_one_cluster=bool((cl >= 2).any())))
    c = pd.DataFrame(out)
    s = [dict(measure="multi-chain entries", n=len(c), k=len(c), fraction=1.0)]
    for col, name in (("label_differs", "H12 label differs (complete/partial/not modelled)"),
                      ("complete_differs", "binary outcome differs (complete vs not)"),
                      ("cluster_differs", "S08 cluster differs")):
        s.append(dict(measure=name, n=len(c), k=int(c[col].sum()), fraction=c[col].mean()))
    both = c[c.label_differs & c.cluster_differs]
    s.append(dict(measure="label and cluster both differ", n=len(c), k=len(both), fraction=len(both) / len(c)))
    sub = c[c.has_two_chains_in_one_cluster]
    s.append(dict(measure="label differs between chains of the same cluster (entries with >=2 chains in one cluster)",
                  n=len(sub), k=int(sub.label_differs_within_a_cluster.sum()),
                  fraction=sub.label_differs_within_a_cluster.mean()))
    for g in ("cl7", "cl4"):
        cc = 7 if g == "cl7" else 4
        tt = t[t.s08_cluster == cc]
        ee = tt.groupby("entry_id").filter(lambda x: len(x) >= 2).groupby("entry_id").complete.nunique()
        s.append(dict(measure=f"{g}: binary outcome differs between the entry's {g} chains (entries with >=2 {g} chains)",
                      n=len(ee), k=int((ee > 1).sum()), fraction=(ee > 1).mean() if len(ee) else np.nan))
    return c, pd.DataFrame(s)


# ------------------------------------------------------------------------------------------------ step 3: models
def ll_bern(y, mu):
    mu = np.clip(mu, 1e-12, 1 - 1e-12)
    return float(np.sum(y * np.log(mu) + (1 - y) * np.log(1 - mu)))


def fit(y, Xd, groups=None):
    """Logit MLE (statsmodels Logit; for the fraction outcome the same estimator is the Papke-Wooldridge fractional
    logit QMLE). Logit's sandwich uses the score (y - mu) x and the Hessian X'diag(mu(1-mu))X, so it stays finite
    when fitted probabilities reach 0/1 (GLM's version divides by mu(1-mu)). Newton, BFGS if Newton fails."""
    if Xd.shape[1] == 0:
        mu = np.full(len(y), y.mean())
        return dict(res=None, llf=ll_bern(y, mu), sep=False, note="")
    Xc = sm.add_constant(Xd.astype(float), has_constant="add")
    kw = dict(cov_type="cluster", cov_kwds={"groups": groups}) if groups is not None else {}
    res, method_note = None, ""
    with warnings.catch_warnings(record=True) as wl:
        warnings.simplefilter("always")
        for method in ("newton", "bfgs"):
            try:
                res = sm.Logit(y, Xc).fit(method=method, maxiter=500, disp=0, **kw)
                method_note = "" if method == "newton" else "bfgs"
                break
            except Exception as ex:  # noqa: BLE001
                method_note = f"{method} error: {ex}"
    if res is None:
        return dict(res=None, llf=np.nan, sep=True, note=f"fit error: {method_note}")
    msgs = [str(w.message) for w in wl]
    perfect = any("separation" in m.lower() or "perfect" in m.lower() for m in msgs)
    conv = bool(res.mle_retvals.get("converged", True))
    big = bool(np.nanmax(np.abs(res.params.values[1:])) > COEF_LIMIT)
    sep = perfect or (not conv) or big
    note = ";".join(x for x, f in (("perfect-separation warning", perfect), ("no convergence", not conv),
                                   (f"|log OR|>{COEF_LIMIT:g}", big), (method_note, bool(method_note))) if f)
    return dict(res=res, llf=ll_bern(y, np.asarray(res.predict())), sep=sep, note=note)


def screen(d, y, preds):
    """Univariate drop rule. Returns kept list and dict of dropped -> reason."""
    kept, dropped = [], {}
    yb = (y >= 1 - 1e-12)
    for p in preds:
        x = d[p].values
        if np.nanstd(x) == 0:
            dropped[p] = "constant"
            continue
        if p in CONTINUOUS:
            a, b = x[yb], x[~yb]
            if len(a) and len(b) and (a.max() < b.min() or b.max() < a.min()):
                dropped[p] = "separation (non-overlapping ranges)"
                continue
        else:
            sep = False
            for lv in (0.0, 1.0):
                yy = y[x == lv]
                if len(yy) and ((yy >= 1 - 1e-12).all() or (yy <= 1e-12).all()):
                    sep = True
                    side = "all complete" if (yy >= 1 - 1e-12).all() else "none complete"
                    dropped[p] = f"separation ({'with' if lv == 1 else 'without'} = {int(lv)}: {len(yy)} chains, {side})"
                    break
            if sep:
                continue
        kept.append(p)
    # exact collinearity in block order
    final = []
    for p in kept:
        M = np.column_stack([np.ones(len(d))] + [d[q].values for q in final + [p]])
        if np.linalg.matrix_rank(M) < M.shape[1]:
            dropped[p] = "collinear with earlier predictors"
        else:
            final.append(p)
    return final, dropped


def r2(llf, ll0):
    return 1 - llf / ll0


def run_models(t, dcoded):
    mb_rows, or_rows, scale_rows, data_rows = [], [], [], []
    for c in MAIN:
        idx = np.where(t.s08_cluster.values == c)[0]
        tc = t.iloc[idx].reset_index(drop=True)
        dc = dcoded.iloc[idx].reset_index(drop=True).copy()
        preds = list(Q) + list(X)
        for v in IMPUTE:
            nmiss = int(dc[v].isna().sum())
            med = float(dc[v].median())
            if nmiss:
                dc[f"{v}_missing"] = dc[v].isna().astype(float)
                preds.append(f"{v}_missing")
            dc[v] = dc[v].fillna(med)
            scale_rows.append(dict(cluster=c, predictor=v, n_missing=nmiss, imputed_median=med))
        preds += list(B)
        raw = dc.copy()
        for p in preds:
            if p in CONTINUOUS:
                mu, sd = float(dc[p].mean()), float(dc[p].std(ddof=0))
                scale_rows.append(dict(cluster=c, predictor=p, mean=mu, sd=sd, unit="per 1 SD"))
                dc[p] = (dc[p] - mu) / sd if sd > 0 else 0.0
        raw.insert(0, "chain_key", tc.chain_key)
        raw.insert(1, "entry_id", tc.entry_id)
        raw.insert(2, "s08_cluster", c)
        raw.insert(3, "complete", tc.complete)
        raw.insert(4, "frac_obs", tc.frac_obs)
        raw["entry_level_row"] = ~tc.entry_id.duplicated()
        data_rows.append(raw)
        for level in ("chain", "entry"):
            s = np.ones(len(tc), bool) if level == "chain" else (~tc.entry_id.duplicated()).values
            ts, ds = tc[s].reset_index(drop=True), dc[s].reset_index(drop=True)
            groups = ts.entry_id.values if level == "chain" else None
            for outcome, ycol in (("complete", "complete"), ("fraction_observed", "frac_obs")):
                y = ts[ycol].values.astype(float)
                kept, dropped = screen(ds, y, preds)
                # joint separation loop
                while True:
                    full = fit(y, ds[kept], groups)
                    if not full["sep"] or not kept:
                        break
                    if full["res"] is None:
                        break
                    co = full["res"].params.drop("const").abs()
                    worst = co.idxmax()
                    dropped[worst] = f"joint separation ({full['note']}; |log OR| {co[worst]:.1f})"
                    kept.remove(worst)
                ll0 = fit(y, ds[[]])["llf"]
                bl = {k: [p for p in kept if block_of(p) == k] for k in "QXB"}
                combos = {"": [], "Q": bl["Q"], "X": bl["X"], "B": bl["B"], "QX": bl["Q"] + bl["X"],
                          "QB": bl["Q"] + bl["B"], "XB": bl["X"] + bl["B"], "QXB": kept}
                fits = {k: (full if k == "QXB" else fit(y, ds[v], groups)) for k, v in combos.items()}
                R = {k: r2(f["llf"], ll0) for k, f in fits.items()}
                R[""] = 0.0
                res = full["res"]
                n_ev = int((y < 1 - 1e-12).sum())
                common = dict(outcome=outcome, cluster=c, level=level, n=len(y), n_entries=ts.entry_id.nunique(),
                              n_incomplete=n_ev, mean_outcome=float(y.mean()), r2_full=R["QXB"],
                              full_model_note=full["note"], n_predictors_full=len(kept))
                others = {"Q": "XB", "X": "QB", "B": "QX", "Q+X": "B"}
                alone = {"Q": "Q", "X": "X", "B": "B", "Q+X": "QX"}
                for blk in ("Q", "X", "B", "Q+X"):
                    cols = bl["Q"] + bl["X"] if blk == "Q+X" else bl[blk]
                    # test of the block entered last
                    stat = dfree = p = np.nan
                    test = ""
                    if cols and res is not None:
                        if level == "chain":
                            names = list(res.params.index)
                            Rm = np.zeros((len(cols), len(names)))
                            for i, cc in enumerate(cols):
                                Rm[i, names.index(cc)] = 1
                            try:
                                wt = res.wald_test(Rm, scalar=True)
                                stat, p, dfree = float(wt.statistic), float(wt.pvalue), len(cols)
                                test = "cluster-robust Wald (block entered last)"
                            except Exception:  # noqa: BLE001
                                test = "Wald failed"
                        else:
                            lr = 2 * (fits["QXB"]["llf"] - fits[others[blk]]["llf"])
                            stat, dfree, p = lr, len(cols), float(chi2.sf(lr, len(cols)))
                            test = "LR (block entered last)"
                    for order in ("first", "last"):
                        g = R[alone[blk]] if order == "first" else R["QXB"] - R[others[blk]]
                        sub_note = fits[alone[blk]]["note"] if order == "first" else fits[others[blk]]["note"]
                        mb_rows.append(dict(common, block=blk, order=order, mcfadden_gain=g,
                                            r2_block_alone=R[alone[blk]], r2_without_block=R[others[blk]],
                                            n_predictors_block=len(cols), predictors=";".join(cols),
                                            test=test, statistic=stat, df=dfree, p_value=p if order == "last" else np.nan,
                                            submodel_note=sub_note))
                # odds ratios
                for pnm in preds:
                    x = ds[pnm].values
                    row = dict(outcome=outcome, cluster=c, level=level, block=block_of(pnm), predictor=pnm,
                               label=LABEL[pnm], unit="per 1 SD" if pnm in CONTINUOUS else "present vs absent",
                               n=len(y), n_incomplete=n_ev)
                    if pnm not in CONTINUOUS:
                        row.update(n_with=int((x == 1).sum()), n_incomplete_with=int(((x == 1) & (y < 1 - 1e-12)).sum()),
                                   n_without=int((x == 0).sum()),
                                   n_incomplete_without=int(((x == 0) & (y < 1 - 1e-12)).sum()))
                    if pnm in kept and res is not None:
                        b, se = float(res.params[pnm]), float(res.bse[pnm])
                        row.update(status="fitted", log_or=b, se=se, odds_ratio=np.exp(b),
                                   ci_low=np.exp(b - 1.96 * se), ci_high=np.exp(b + 1.96 * se),
                                   p_value=float(res.pvalues[pnm]))
                    else:
                        row.update(status="dropped: " + dropped.get(pnm, "not fitted"))
                    or_rows.append(row)
    mb = pd.DataFrame(mb_rows)
    orr = pd.DataFrame(or_rows)
    # one BH family over every p-value of the step
    pv = pd.concat([mb.loc[mb.p_value.notna(), "p_value"].rename("p").to_frame().assign(src="mb"),
                    orr.loc[orr.p_value.notna(), "p_value"].rename("p").to_frame().assign(src="or")])
    q = multipletests(pv.p.values, method="fdr_bh")[1]
    pv["q"] = q
    mb["q_value"] = np.nan
    orr["q_value"] = np.nan
    mb.loc[pv[pv.src == "mb"].index, "q_value"] = pv[pv.src == "mb"].q.values
    orr.loc[pv[pv.src == "or"].index, "q_value"] = pv[pv.src == "or"].q.values
    return mb, orr, pd.DataFrame(scale_rows), pd.concat(data_rows, ignore_index=True), len(pv)


def verdicts(mb):
    rows = []
    for outcome in ("complete", "fraction_observed"):
        for level in ("chain", "entry"):
            cells = []
            for c in MAIN:
                for order in ("first", "last"):
                    m = mb[(mb.outcome == outcome) & (mb.level == level) & (mb.cluster == c) & (mb.order == order)]
                    qx = float(m[m.block == "Q+X"].mcfadden_gain.iloc[0])
                    bb = float(m[m.block == "B"].mcfadden_gain.iloc[0])
                    cells.append(dict(outcome=outcome, level=level, cluster=c, order=order, gain_QX=qx, gain_B=bb,
                                      QX_exceeds_B=qx > bb, B_exceeds_QX=bb > qx))
            cdf = pd.DataFrame(cells)
            v = "artifact" if cdf.QX_exceeds_B.all() else "biology" if cdf.B_exceeds_QX.all() else "mixed"
            cdf["level_verdict"] = v
            rows.append(cdf)
    vt = pd.concat(rows, ignore_index=True)
    step = {}
    for outcome in ("complete", "fraction_observed"):
        lv = vt[vt.outcome == outcome].groupby("level").level_verdict.first()
        step[outcome] = lv.iloc[0] if lv.nunique() == 1 and lv.iloc[0] != "mixed" else "mixed"
    vt["step_verdict"] = vt.outcome.map(step)
    return vt, step


def main():
    t, mask, pos = load()
    # cross-check with the S16 label
    agree = int((t.h12_status == t.h12_label).sum())
    prof = profile(t, mask, pos)
    prof.to_csv(OUT / "completeness_profile.csv", index=False)
    sc = status_counts(t)
    sc.to_csv(OUT / "h12_status_counts.csv", index=False)
    pe = partial_ends(t, mask, pos)
    pe.to_csv(OUT / "partial_ends.csv", index=False)
    t[["chain_key", "entry_id", "s08_cluster", "cluster_group", "mutation_status", "h12_label", "h12_status",
       "n_obs_window", "frac_obs", "complete"]].rename(columns={"h12_label": "h12_label_s16"}).to_csv(
        OUT / "chain_outcomes.csv", index=False)
    conc, csum = concordance(t)
    conc.to_csv(OUT / "concordance.csv", index=False)
    csum.to_csv(OUT / "concordance_summary.csv", index=False)

    dcoded = code_predictors(t)
    mb, orr, scale, mdata, n_tests = run_models(t, dcoded)
    mb.to_csv(OUT / "model_blocks.csv", index=False)
    orr.to_csv(OUT / "odds_ratios.csv", index=False)
    scale.to_csv(OUT / "predictor_scaling.csv", index=False)
    mdata.to_csv(OUT / "model_data.csv", index=False)
    vt, step = verdicts(mb)
    vt.to_csv(OUT / "verdict_table.csv", index=False)

    # verdict.md (<= 10 lines)
    def gains(outcome, level, c):
        m = mb[(mb.outcome == outcome) & (mb.level == level) & (mb.cluster == c) & (mb.order == "first")]
        g = m.set_index("block").mcfadden_gain
        return f"Q+X {g['Q+X']:.3f} vs B {g['B']:.3f}"
    lv = vt.groupby(["outcome", "level"]).level_verdict.first()
    lines = [
        f"# SV01a verdict: **{step['complete']}** (H12 complete, window 538-547; [DEFAULT] rule, chain and entry level)",
        f"- Chain level ({lv[('complete', 'chain')]}): McFadden gain entered first, cl7 {gains('complete', 'chain', 7)}; "
        f"cl4 {gains('complete', 'chain', 4)}.",
        f"- Entry level ({lv[('complete', 'entry')]}): cl7 {gains('complete', 'entry', 7)}; cl4 {gains('complete', 'entry', 4)}.",
        "- With the same rows in every model, the 'entered last' comparison of Q+X vs B equals the 'entered first' one "
        "(both reduce to R2(Q+X) vs R2(B)); values for both orders are in model_blocks.csv.",
        f"- Sensitivity, fraction of 538-547 observed: **{step['fraction_observed']}** (chain "
        f"{lv[('fraction_observed', 'chain')]}, entry {lv[('fraction_observed', 'entry')]}).",
        "- Events vs fitted predictors (complete outcome): " + "; ".join(
            f"cl{c} {lvl} {int(m.n_incomplete.iloc[0])} incomplete / {int(m.n_predictors_full.iloc[0])} predictors "
            f"(full R2 {m.r2_full.iloc[0]:.3f})"
            for (c, lvl), m in mb[mb.outcome == "complete"].groupby(["cluster", "level"], sort=False)) + ".",
        "- Dropped at chain level (complete outcome; full list per level in odds_ratios.csv): " + "; ".join(
            f"cl{c} " + ", ".join(f"{r.predictor} [{r.block}, {r.status.replace('dropped: ', '').split(' (')[0]}]"
                                  for r in g.itertuples())
            for c, g in orr[(orr.outcome == "complete") & (orr.level == "chain") & (orr.status != "fitted")]
            .groupby("cluster", sort=False)) + ".",
        "- First pass without lattice contacts and mobility; SV01b refits with them and re-issues the verdict.",
    ]
    (OUT / "verdict.md").write_text("\n".join(lines) + "\n")

    nums = dict(n_chains=len(t), n_entries=int(t.entry_id.nunique()), mask_vs_s16_label_agree=agree,
                n_tests_bh=n_tests, step_verdict=step,
                level_verdicts={f"{k[0]}|{k[1]}": v for k, v in lv.items()})
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1))

    tz = timezone(timedelta(hours=5))
    env = [f"# SV01a environment, {datetime.now(tz).strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}",
           f"numpy {np.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}",
           f"statsmodels {statsmodels.__version__}"]
    import matplotlib
    env.append(f"matplotlib {matplotlib.__version__}")
    env.append("nothing installed, upgraded or downgraded by SV01a; no network used")
    (OUT / "env.txt").write_text("\n".join(env) + "\n")
    print(json.dumps(nums, indent=1))
    print(vt.to_string())


if __name__ == "__main__":
    main()
