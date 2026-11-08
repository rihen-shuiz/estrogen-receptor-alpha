"""SV01b (RUN_SV1028) steps 5 and 7: local mobility as a second outcome, lattice terms in block X, refit, re-issued verdict.

SV01b is a refit of the SV01a models (scripts/sv/SV01a_completeness.py: same chains, outcome, predictors, imputation,
scaling, levels, drop rule, fit, McFadden gains, tests and BH) with two lattice terms added to block X, plus the new
mobility outcome. The SV01a functions are imported, not copied; a check refits the SV01a predictor set with the code
below and compares it with results/sv/SV01a/model_blocks.csv and odds_ratios.csv (reproduce_SV01a.json).

Rules fixed before any SV01b model was fitted:
- Lattice terms (step 6, scripts/sv/common/lattice_contacts.py, results/sv/SV01b/lattice_contacts_per_chain.csv):
  lattice_asu = heavy atoms of the chain's residues 528-550 within 4.0 A of another chain's ASU copy;
  lattice_sym = the same for symmetry copies (other operations or lattice translations, the chain's own image
  included). The chain's own bound peptide copy is not a lattice contact (it is 'peptide' in block B). Both are
  continuous, z-scored within the cluster like the other continuous predictors (odds ratios per 1 SD), and appended
  to block X after the SV01a X predictors.
- Mobility outcome (step 5): per-residue C-alpha B (data/sv/backbone/backbone.npz, SV00; the existing altloc policy,
  highest occupancy then blank/A) of the observed residues 528-537 and the observed H12 residues 538-547, as
  z = (B_res - mean_core) / sd_core with the same chain's core v3 C-alpha B (core_ranges of
  esr_analysis/resources/config/alignment_core_v3.json, the S03 core; positions 338-362, 371-396, 420-439, 441-456,
  496-531, so 528-531 belong to both the core and the window, as the step names core v3 unchanged). Valid B =
  finite and >= 0, SD with ddof 0 (the existing regional-B producer, esr_analysis/_regional_bfactors.py).
  mobility_z = mean z over every observed residue of 528-547 (one value per chain); undefined if no residue of
  528-547 is observed or the core has < 2 valid B or SD 0 (those chains leave the mobility models only, listed).
  The loop and H12 means are saved but not modelled.
- Mobility models: the same blocks, imputation and scaling (computed on all chains of the cluster, as for the
  completeness models), rows = chains with mobility_z defined. Chain level: OLS with entry-clustered (sandwich) errors;
  entry level: OLS on one chain per entry (the entry's first matrix chain in the cluster among the rows with the
  outcome). Drop rule: constant or exactly collinear predictors (separation does not apply to a linear model).
  Gains: OLS R2 of the block alone (entered first) and R2(full) - R2(without the block) (entered last). Tests of the
  block entered last: cluster-robust Wald chi2 at chain level, nested F-test at entry level. Coefficients (z units
  per 1 SD or present vs absent) with Wald 95 % CI b +- 1.96 SE and normal p-values (as SV01a).
- TLS sensitivity (step 5): the mobility models refitted without chains of TLS-refined entries (SV00 tls_refined).
- Outcomes refitted: complete (primary), fraction observed (SV01a sensitivity), mobility, mobility without TLS.
  One BH family over every p-value of SV01b.
- [DEFAULT] verdict rule per outcome and level as SV01a ("artifact" if the Q+X gain exceeds the B gain in both entry
  orders and both clusters, "biology" if the reverse, else "mixed"). Re-issued step verdict: "artifact" or "biology"
  only when that condition holds for both outcomes (complete, mobility) at both levels (chain, entry); otherwise
  "mixed". Fraction observed and mobility without TLS are sensitivity verdicts. Step 8 (medoid H12 clash test) runs
  only if the re-issued step verdict is "mixed" (scripts/sv/SV01b_clash.py reads numbers.json).

Outputs (results/sv/SV01b/): mobility_per_chain.csv, mobility_per_residue.csv, model_data.csv, predictor_scaling.csv,
model_blocks.csv, odds_ratios.csv, verdict_table.csv, verdict.md, numbers.json, reproduce_SV01a.json, env.txt.
Run: .venv/bin/python -B scripts/sv/SV01b_models.py
"""
from __future__ import annotations

import json
import platform
import sys
import warnings
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.api as sm
from scipy.stats import chi2, norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV01a_completeness as A  # noqa: E402  (SV01a model code: load, code_predictors, screen, fit, block_of, ...)

OUT = ROOT / "results/sv/SV01b"
OUT.mkdir(parents=True, exist_ok=True)
CORE_CFG = ROOT / "esr_analysis/resources/config/alignment_core_v3.json"
MOB = (528, 547)
LOOP = (528, 537)
LATTICE = ["lattice_asu", "lattice_sym"]
A.LABEL.update({"lattice_asu": "lattice contacts, ASU partner (atoms of 528-550)",
                "lattice_sym": "lattice contacts, symmetry mate (atoms of 528-550)"})
A.CONTINUOUS.update(LATTICE)
X_SV01A = list(A.X)
X_SV01B = X_SV01A + LATTICE

# outcome specs: name, column, family, row filter (applied after imputation and scaling on all cluster chains)
SPECS = [
    dict(outcome="complete", col="complete", family="logit", rows=None, role="primary"),
    dict(outcome="fraction_observed", col="frac_obs", family="logit", rows=None, role="sensitivity of complete"),
    dict(outcome="mobility", col="mobility_z", family="linear", rows="mobility_defined", role="second outcome"),
    dict(outcome="mobility_no_tls", col="mobility_z", family="linear", rows="mobility_defined_no_tls",
         role="sensitivity of mobility (TLS-refined entries removed)"),
]


# ------------------------------------------------------------------------------------------------ step 5: mobility
def mobility(t, mask, pos):
    bb = np.load(ROOT / "data/sv/backbone/backbone.npz")
    assert (bb["chain_keys"] == t.chain_key.values).all() and (bb["positions"] == pos).all()
    ca = list(bb["atoms"]).index("CA")
    B = np.where(mask, bb["b_factor"][:, :, ca], np.nan)
    # cross-check with the regenerated Q04 residue mapping (B_factor_A2 of the selected C-alpha)
    mp = pd.read_csv(ROOT / "data/sv/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_ca", "B_factor_A2"])
    mp = mp[mp.chain_key.isin(t.chain_key) & mp.has_ca.astype(bool) & mp.esr1_position.between(pos[0], pos[-1])]
    ri = pd.Series(np.arange(len(t)), index=t.chain_key)
    bm = B[ri[mp.chain_key].values, (mp.esr1_position.values - pos[0]).astype(int)]
    b_check = dict(n_compared=int(len(mp)), max_abs_diff=float(np.nanmax(np.abs(bm - mp.B_factor_A2.values))),
                   n_nan_in_backbone=int(np.isnan(bm).sum()))
    cfg = json.loads(CORE_CFG.read_text())
    core = np.zeros(len(pos), bool)
    for lo, hi in cfg["core_ranges"]:
        core |= (pos >= lo) & (pos <= hi)
    win = (pos >= MOB[0]) & (pos <= MOB[1])
    loop = (pos >= LOOP[0]) & (pos <= LOOP[1])
    rows, rrows = [], []
    for i, ck in enumerate(t.chain_key):
        cb = B[i, core]
        cb = cb[np.isfinite(cb) & (cb >= 0)]
        mu, sd = (cb.mean(), cb.std(ddof=0)) if len(cb) else (np.nan, np.nan)
        ok_core = len(cb) >= 2 and sd > 0
        z = (B[i] - mu) / sd if ok_core else np.full(len(pos), np.nan)
        zv = np.where(np.isfinite(B[i]) & (B[i] >= 0), z, np.nan)
        wz, lz, hz = zv[win], zv[loop], zv[win & ~loop]
        rows.append(dict(chain_key=ck, entry_id=t.entry_id[i], s08_cluster=int(t.s08_cluster[i]),
                         cluster_group=t.cluster_group[i], h12_status=t.h12_status[i],
                         tls_refined=str(t.tls_refined[i]) == "True", n_core_ca_valid_b=len(cb), core_mean_b=mu, core_sd_b=sd,
                         n_obs_528_537=int(np.isfinite(lz).sum()), n_obs_538_547=int(np.isfinite(hz).sum()),
                         mobility_z=np.nanmean(wz) if np.isfinite(wz).any() else np.nan,
                         mobility_loop_z=np.nanmean(lz) if np.isfinite(lz).any() else np.nan,
                         mobility_h12_z=np.nanmean(hz) if np.isfinite(hz).any() else np.nan,
                         core_ok=ok_core))
        for j in np.where(win & np.isfinite(zv))[0]:
            rrows.append(dict(chain_key=ck, esr1_position=int(pos[j]), ca_b=B[i, j], b_z=zv[j]))
    m = pd.DataFrame(rows)
    m["mobility_defined"] = m.mobility_z.notna()
    return m, pd.DataFrame(rrows), b_check, cfg["core_ranges"]


def lattice_by_group(t, lat):
    """Descriptive table of the lattice terms by cluster group x H12 status (where on 528-550 the contacts sit)."""
    pr = pd.read_csv(OUT / "lattice_contacts_per_residue.csv")
    seg = {"528_537": (528, 537), "538_545": (538, 545), "546_550": (546, 550)}
    d = t[["chain_key", "cluster_group", "h12_status"]].merge(lat, on="chain_key")
    for k, (lo, hi) in seg.items():
        s = pr[pr.esr1_position.between(lo, hi)].groupby("chain_key")[["asu_atoms", "sym_atoms", "n_heavy_atoms"]].sum()
        s.columns = [f"{c}_{k}" for c in s.columns]
        d = d.merge(s, left_on="chain_key", right_index=True, how="left")
    d = d.fillna({c: 0 for c in d.columns if c.endswith(tuple(seg))})
    rows = []
    for (g, st), x in list(d.groupby(["cluster_group", "h12_status"])) + [((g, "all"), x) for g, x in
                                                                           d.groupby("cluster_group")]:
        r = dict(cluster_group=g, h12_status=st, n_chains=len(x), n_entries=x.entry_id.nunique(),
                 chains_asu_contact=int((x.asu_atoms > 0).sum()), chains_sym_contact=int((x.sym_atoms > 0).sum()),
                 mean_asu_atoms=x.asu_atoms.mean(), mean_sym_atoms=x.sym_atoms.mean(),
                 median_sym_atoms=x.sym_atoms.median(), mean_modelled_heavy_atoms_528_550=x.n_heavy_atoms.mean())
        for k in seg:
            r.update({f"mean_sym_atoms_{k}": x[f"sym_atoms_{k}"].mean(), f"mean_asu_atoms_{k}": x[f"asu_atoms_{k}"].mean(),
                      f"mean_modelled_heavy_atoms_{k}": x[f"n_heavy_atoms_{k}"].mean()})
        rows.append(r)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ fitting
def fit_ols(y, Xd, groups=None):
    Xc = sm.add_constant(Xd.astype(float), has_constant="add") if Xd.shape[1] else np.ones((len(y), 1))
    if groups is not None:
        res = sm.OLS(y, Xc).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(groups)[0]})
    else:
        res = sm.OLS(y, Xc).fit()
    return dict(res=res if Xd.shape[1] else None, res_any=res, r2=float(res.rsquared), sep=False, note="")


def screen_linear(d, preds):
    kept, dropped = [], {}
    for p in preds:
        if np.nanstd(d[p].values) == 0:
            dropped[p] = "constant"
        else:
            kept.append(p)
    final = []
    for p in kept:
        M = np.column_stack([np.ones(len(d))] + [d[q].values for q in final + [p]])
        if np.linalg.matrix_rank(M) < M.shape[1]:
            dropped[p] = "collinear with earlier predictors"
        else:
            final.append(p)
    return final, dropped


def run_models(t, dcoded, xlist, specs):
    """SV01a run_models generalised to linear outcomes and row subsets; with xlist = SV01a X and the two logit
    outcomes it is the SV01a computation."""
    A.X = list(xlist)  # block_of() and the predictor order read the module-level lists of SV01a
    mb_rows, or_rows, scale_rows, data_rows = [], [], [], []
    for c in A.MAIN:
        idx = np.where(t.s08_cluster.values == c)[0]
        tc = t.iloc[idx].reset_index(drop=True)
        dc = dcoded.iloc[idx].reset_index(drop=True).copy()
        preds = list(A.Q) + list(xlist)
        for v in A.IMPUTE:
            nmiss = int(dc[v].isna().sum())
            med = float(dc[v].median())
            if nmiss:
                dc[f"{v}_missing"] = dc[v].isna().astype(float)
                preds.append(f"{v}_missing")
            dc[v] = dc[v].fillna(med)
            scale_rows.append(dict(cluster=c, predictor=v, n_missing=nmiss, imputed_median=med))
        preds += list(A.B)
        raw = dc.copy()
        for p in preds:
            if p in A.CONTINUOUS:
                mu, sd = float(dc[p].mean()), float(dc[p].std(ddof=0))
                scale_rows.append(dict(cluster=c, predictor=p, mean=mu, sd=sd, unit="per 1 SD"))
                dc[p] = (dc[p] - mu) / sd if sd > 0 else 0.0
        raw.insert(0, "chain_key", tc.chain_key)
        raw.insert(1, "entry_id", tc.entry_id)
        raw.insert(2, "s08_cluster", c)
        for k, col in enumerate(["complete", "frac_obs", "mobility_z", "tls_refined"]):
            if col in tc:
                raw.insert(3 + k, col, tc[col].values)
        for spec in specs:
            rows = np.ones(len(tc), bool) if spec["rows"] is None else tc[spec["rows"]].values.astype(bool)
            raw[f"entry_level_row_{spec['outcome']}"] = rows & ~tc.entry_id.where(rows).duplicated().values \
                if spec["rows"] else ~tc.entry_id.duplicated().values
        data_rows.append(raw)
        for level in ("chain", "entry"):
            for spec in specs:
                outcome, ycol, fam = spec["outcome"], spec["col"], spec["family"]
                rows = np.ones(len(tc), bool) if spec["rows"] is None else tc[spec["rows"]].values.astype(bool)
                tr, dr = tc[rows].reset_index(drop=True), dc[rows].reset_index(drop=True)
                s = np.ones(len(tr), bool) if level == "chain" else (~tr.entry_id.duplicated()).values
                ts, ds = tr[s].reset_index(drop=True), dr[s].reset_index(drop=True)
                groups = ts.entry_id.values if level == "chain" else None
                y = ts[ycol].values.astype(float)
                if fam == "logit":
                    kept, dropped = A.screen(ds, y, preds)
                    while True:
                        full = A.fit(y, ds[kept], groups)
                        if not full["sep"] or not kept or full["res"] is None:
                            break
                        co = full["res"].params.drop("const").abs()
                        worst = co.idxmax()
                        dropped[worst] = f"joint separation ({full['note']}; |log OR| {co[worst]:.1f})"
                        kept.remove(worst)
                    fitter = lambda cols: A.fit(y, ds[cols], groups)  # noqa: E731
                    ll0 = A.fit(y, ds[[]])["llf"]
                    score = lambda f: A.r2(f["llf"], ll0)  # noqa: E731
                    gain_measure = "McFadden R2"
                else:
                    kept, dropped = screen_linear(ds, preds)
                    full = fit_ols(y, ds[kept], groups)
                    fitter = lambda cols: fit_ols(y, ds[cols], groups)  # noqa: E731
                    score = lambda f: f["r2"]  # noqa: E731
                    gain_measure = "R2 (OLS)"
                bl = {k: [p for p in kept if A.block_of(p) == k] for k in "QXB"}
                combos = {"": [], "Q": bl["Q"], "X": bl["X"], "B": bl["B"], "QX": bl["Q"] + bl["X"],
                          "QB": bl["Q"] + bl["B"], "XB": bl["X"] + bl["B"], "QXB": kept}
                fits = {k: (full if k == "QXB" else fitter(v)) for k, v in combos.items()}
                R = {k: score(f) for k, f in fits.items()}
                R[""] = 0.0
                res = full["res"]
                n_ev = int((y < 1 - 1e-12).sum()) if fam == "logit" else np.nan
                common = dict(outcome=outcome, role=spec["role"], family=fam, gain_measure=gain_measure, cluster=c,
                              level=level, n=len(y), n_entries=ts.entry_id.nunique(), n_incomplete=n_ev,
                              mean_outcome=float(y.mean()), sd_outcome=float(y.std(ddof=0)), r2_full=R["QXB"],
                              full_model_note=full["note"], n_predictors_full=len(kept))
                others = {"Q": "XB", "X": "QB", "B": "QX", "Q+X": "B"}
                alone = {"Q": "Q", "X": "X", "B": "B", "Q+X": "QX"}
                for blk in ("Q", "X", "B", "Q+X"):
                    cols = bl["Q"] + bl["X"] if blk == "Q+X" else bl[blk]
                    stat = dfree = p = np.nan
                    test = ""
                    if cols and res is not None:
                        if level == "chain":
                            names = list(res.params.index)
                            Rm = np.zeros((len(cols), len(names)))
                            for i, cc in enumerate(cols):
                                Rm[i, names.index(cc)] = 1
                            try:
                                wt = res.wald_test(Rm, scalar=True, use_f=False) if fam == "linear" else \
                                    res.wald_test(Rm, scalar=True)
                                stat, p, dfree = float(wt.statistic), float(wt.pvalue), len(cols)
                                test = "cluster-robust Wald (block entered last)"
                            except Exception:  # noqa: BLE001
                                test = "Wald failed"
                        elif fam == "logit":
                            lr = 2 * (fits["QXB"]["llf"] - fits[others[blk]]["llf"])
                            stat, dfree, p = lr, len(cols), float(chi2.sf(lr, len(cols)))
                            test = "LR (block entered last)"
                        else:
                            f_, p_, d_ = fits["QXB"]["res_any"].compare_f_test(fits[others[blk]]["res_any"])
                            stat, dfree, p = float(f_), int(d_), float(p_)
                            test = "F (block entered last)"
                    for order in ("first", "last"):
                        g = R[alone[blk]] if order == "first" else R["QXB"] - R[others[blk]]
                        sub_note = fits[alone[blk]]["note"] if order == "first" else fits[others[blk]]["note"]
                        mb_rows.append(dict(common, block=blk, order=order, gain=g,
                                            r2_block_alone=R[alone[blk]], r2_without_block=R[others[blk]],
                                            n_predictors_block=len(cols), predictors=";".join(cols),
                                            test=test, statistic=stat, df=dfree,
                                            p_value=p if order == "last" else np.nan, submodel_note=sub_note))
                for pnm in preds:
                    x = ds[pnm].values
                    row = dict(outcome=outcome, role=spec["role"], family=fam, cluster=c, level=level,
                               block=A.block_of(pnm), predictor=pnm, label=A.LABEL[pnm],
                               unit="per 1 SD" if pnm in A.CONTINUOUS else "present vs absent",
                               effect="odds ratio" if fam == "logit" else "coefficient (z units)",
                               n=len(y), n_incomplete=n_ev)
                    if pnm not in A.CONTINUOUS:
                        row.update(n_with=int((x == 1).sum()), n_without=int((x == 0).sum()))
                        if fam == "logit":
                            row.update(n_incomplete_with=int(((x == 1) & (y < 1 - 1e-12)).sum()),
                                       n_incomplete_without=int(((x == 0) & (y < 1 - 1e-12)).sum()))
                    if pnm in kept and res is not None:
                        b, se = float(res.params[pnm]), float(res.bse[pnm])
                        if fam == "logit":
                            row.update(status="fitted", coef=b, se=se, log_or=b, odds_ratio=np.exp(b),
                                       estimate=np.exp(b), ci_low=np.exp(b - 1.96 * se),
                                       ci_high=np.exp(b + 1.96 * se), p_value=float(res.pvalues[pnm]))
                        else:
                            row.update(status="fitted", coef=b, se=se, estimate=b, ci_low=b - 1.96 * se,
                                       ci_high=b + 1.96 * se, p_value=float(2 * norm.sf(abs(b / se))))
                    else:
                        row.update(status="dropped: " + dropped.get(pnm, "not fitted"))
                    or_rows.append(row)
    A.X = X_SV01A
    mb, orr = pd.DataFrame(mb_rows), pd.DataFrame(or_rows)
    pv = pd.concat([mb.loc[mb.p_value.notna(), "p_value"].rename("p").to_frame().assign(src="mb"),
                    orr.loc[orr.p_value.notna(), "p_value"].rename("p").to_frame().assign(src="or")])
    pv["q"] = multipletests(pv.p.values, method="fdr_bh")[1]
    mb["q_value"] = np.nan
    orr["q_value"] = np.nan
    mb.loc[pv[pv.src == "mb"].index, "q_value"] = pv[pv.src == "mb"].q.values
    orr.loc[pv[pv.src == "or"].index, "q_value"] = pv[pv.src == "or"].q.values
    return mb, orr, pd.DataFrame(scale_rows), pd.concat(data_rows, ignore_index=True), len(pv)


def verdicts(mb, outcomes):
    rows = []
    for outcome in outcomes:
        for level in ("chain", "entry"):
            cells = []
            for c in A.MAIN:
                for order in ("first", "last"):
                    m = mb[(mb.outcome == outcome) & (mb.level == level) & (mb.cluster == c) & (mb.order == order)]
                    qx = float(m[m.block == "Q+X"].gain.iloc[0])
                    bb = float(m[m.block == "B"].gain.iloc[0])
                    cells.append(dict(outcome=outcome, level=level, cluster=c, order=order, gain_QX=qx, gain_B=bb,
                                      QX_exceeds_B=qx > bb, B_exceeds_QX=bb > qx))
            cdf = pd.DataFrame(cells)
            cdf["level_verdict"] = ("artifact" if cdf.QX_exceeds_B.all() else
                                    "biology" if cdf.B_exceeds_QX.all() else "mixed")
            rows.append(cdf)
    vt = pd.concat(rows, ignore_index=True)
    lv = vt.groupby(["outcome", "level"]).level_verdict.first()

    def combine(keys):
        vals = {lv[k] for k in keys}
        return vals.pop() if len(vals) == 1 and next(iter(vals)) != "mixed" else "mixed"
    per_outcome = {o: combine([(o, "chain"), (o, "entry")]) for o in outcomes}
    step = combine([(o, lvl) for o in ("complete", "mobility") for lvl in ("chain", "entry")])
    vt["outcome_verdict"] = vt.outcome.map(per_outcome)
    vt["step_verdict"] = step
    return vt, per_outcome, step, lv


def reproduce_sv01a(t, dcoded):
    """Refit the SV01a predictor set with this code and compare with the saved SV01a tables."""
    mb, orr, *_ = run_models(t, dcoded, X_SV01A, SPECS[:2])
    a_mb = pd.read_csv(ROOT / "results/sv/SV01a/model_blocks.csv")
    a_or = pd.read_csv(ROOT / "results/sv/SV01a/odds_ratios.csv")
    k = ["outcome", "cluster", "level", "block", "order"]
    j = mb.merge(a_mb, on=k, suffixes=("", "_a"))
    ko = ["outcome", "cluster", "level", "predictor"]
    jo = orr.merge(a_or, on=ko, suffixes=("", "_a"))
    fitted = jo[jo.status == "fitted"]
    return dict(model_block_rows=len(j), model_block_rows_sv01a=len(a_mb),
                max_abs_diff_gain=float(np.nanmax(np.abs(j.gain - j.mcfadden_gain))),
                max_abs_diff_block_p=float(np.nanmax(np.abs(j.p_value - j.p_value_a))),
                odds_rows=len(jo), odds_rows_sv01a=len(a_or), status_equal=int((jo.status == jo.status_a).sum()),
                max_abs_diff_log_or=float(np.nanmax(np.abs(fitted.log_or - fitted.log_or_a))),
                max_abs_diff_se=float(np.nanmax(np.abs(fitted.se - fitted.se_a))))


def main():
    t, mask, pos = A.load()
    mob, mob_res, b_check, core_ranges = mobility(t, mask, pos)
    mob.to_csv(OUT / "mobility_per_chain.csv", index=False)
    mob_res.to_csv(OUT / "mobility_per_residue.csv", index=False)
    lat = pd.read_csv(OUT / "lattice_contacts_per_chain.csv")
    assert (lat.chain_key.values == t.chain_key.values).all()
    t["mobility_z"] = mob.mobility_z.values
    t["mobility_defined"] = mob.mobility_defined.values
    t["mobility_defined_no_tls"] = (mob.mobility_defined & ~mob.tls_refined).values
    lattice_by_group(t, lat).to_csv(OUT / "lattice_contacts_by_group.csv", index=False)
    dcoded = A.code_predictors(t)
    dcoded["lattice_asu"] = lat.asu_atoms.astype(float).values
    dcoded["lattice_sym"] = lat.sym_atoms.astype(float).values

    rep = reproduce_sv01a(t, dcoded)
    (OUT / "reproduce_SV01a.json").write_text(json.dumps(rep, indent=1))
    print("reproduce SV01a:", rep)

    mb, orr, scale, mdata, n_tests = run_models(t, dcoded, X_SV01B, SPECS)
    mb.to_csv(OUT / "model_blocks.csv", index=False)
    orr.to_csv(OUT / "odds_ratios.csv", index=False)
    scale.to_csv(OUT / "predictor_scaling.csv", index=False)
    mdata.to_csv(OUT / "model_data.csv", index=False)
    outcomes = [s["outcome"] for s in SPECS]
    vt, per_outcome, step, lv = verdicts(mb, outcomes)
    vt.to_csv(OUT / "verdict_table.csv", index=False)

    a_num = json.loads((ROOT / "results/sv/SV01a/numbers.json").read_text())
    a_mb = pd.read_csv(ROOT / "results/sv/SV01a/model_blocks.csv")

    def a_qx(c, lvl):
        return float(a_mb[(a_mb.outcome == "complete") & (a_mb.cluster == c) & (a_mb.level == lvl)
                          & (a_mb.order == "first") & (a_mb.block == "Q+X")].mcfadden_gain.iloc[0])

    def gains(outcome, level, c):
        m = mb[(mb.outcome == outcome) & (mb.level == level) & (mb.cluster == c) & (mb.order == "first")]
        g = m.set_index("block").gain
        return f"Q+X {g['Q+X']:.3f} vs B {g['B']:.3f}"

    def lat_or(outcome, level, c):
        o = orr[(orr.outcome == outcome) & (orr.level == level) & (orr.cluster == c) & orr.predictor.isin(LATTICE)]
        return ", ".join(f"{r.predictor.split('_')[1]} " + (f"{r.estimate:.3g} [{r.ci_low:.3g}, {r.ci_high:.3g}]"
                                                             if r.status == "fitted" else "dropped")
                         for r in o.itertuples())

    nm = mob[mob.s08_cluster.isin(A.MAIN)]
    lines = [
        f"# SV01b verdict: **{step}** (re-issued; outcomes H12 complete 538-547 and local mobility 528-547, "
        "lattice terms in block X; [DEFAULT] rule, chain and entry level)",
        f"- Complete ({per_outcome['complete']}; SV01a {a_num['step_verdict']['complete']}): chain "
        f"({lv[('complete', 'chain')]}) cl7 {gains('complete', 'chain', 7)}, cl4 {gains('complete', 'chain', 4)}; "
        f"entry ({lv[('complete', 'entry')]}) cl7 {gains('complete', 'entry', 7)}, cl4 {gains('complete', 'entry', 4)}"
        " (McFadden gain entered first).",
        f"- Mobility ({per_outcome['mobility']}): chain ({lv[('mobility', 'chain')]}) cl7 "
        f"{gains('mobility', 'chain', 7)}, cl4 {gains('mobility', 'chain', 4)}; entry ({lv[('mobility', 'entry')]}) "
        f"cl7 {gains('mobility', 'entry', 7)}, cl4 {gains('mobility', 'entry', 4)} (OLS R2 gain entered first).",
        f"- Sensitivity: fraction observed **{per_outcome['fraction_observed']}** (chain "
        f"{lv[('fraction_observed', 'chain')]}, entry {lv[('fraction_observed', 'entry')]}); mobility without "
        f"TLS-refined entries **{per_outcome['mobility_no_tls']}** (chain {lv[('mobility_no_tls', 'chain')]}, entry "
        f"{lv[('mobility_no_tls', 'entry')]}).",
        "- Entered last, Q+X vs B gives the same ordering as entered first (same rows in every model: both reduce to "
        "R2(Q+X) vs R2(B)); values for both orders are in model_blocks.csv and verdict_table.csv.",
        "- Complete, Q+X gain entered first without -> with the lattice terms (SV01a -> SV01b): " + "; ".join(
            f"cl{c} {lvl} {a_qx(c, lvl):.3f} -> "
            f"{float(mb[(mb.outcome == 'complete') & (mb.cluster == c) & (mb.level == lvl) & (mb.order == 'first') & (mb.block == 'Q+X')].gain.iloc[0]):.3f}"
            for lvl in ("chain", "entry") for c in A.MAIN) + ".",
        f"- Lattice terms, complete outcome, OR per 1 SD (chain level): cl7 {lat_or('complete', 'chain', 7)}; "
        f"cl4 {lat_or('complete', 'chain', 4)}.",
        f"- Lattice terms, mobility, coefficient per 1 SD (chain level): cl7 {lat_or('mobility', 'chain', 7)}; "
        f"cl4 {lat_or('mobility', 'chain', 4)}.",
        f"- Mobility defined for {int(nm.mobility_defined.sum())}/{len(nm)} cl7+cl4 chains; of them in entries without TLS refinement "
        f"{int((nm.mobility_defined & ~nm.tls_refined).sum())}.",
        "- Lattice contacts of 528-550 are counted on modelled atoms, so they are partly a consequence of the "
        "completeness outcome; step 8 tests the medoid H12 positions instead"
        + (" (run: verdict mixed)." if step == "mixed" else " (not run: verdict not mixed)."),
    ]
    (OUT / "verdict.md").write_text("\n".join(lines) + "\n")
    nums = dict(n_chains=len(t), n_entries=int(t.entry_id.nunique()), n_tests_bh=n_tests, step_verdict=step,
                outcome_verdicts=per_outcome, level_verdicts={f"{k[0]}|{k[1]}": v for k, v in lv.items()},
                sv01a_step_verdict=a_num["step_verdict"], run_step8=step == "mixed",
                b_crosscheck=b_check, core_ranges=core_ranges,
                mobility_defined={g: int(mob[(mob.cluster_group == g)].mobility_defined.sum()) for g in ("cl7", "cl4")},
                mobility_undefined_chains=mob.chain_key[~mob.mobility_defined & mob.s08_cluster.isin(A.MAIN)].tolist(),
                core_not_ok=mob.chain_key[~mob.core_ok].tolist())
    (OUT / "numbers.json").write_text(json.dumps(nums, indent=1))

    tz = timezone(timedelta(hours=5))
    import gemmi
    import matplotlib
    env = [f"# SV01b environment, {datetime.now(tz).strftime('%Y-%m-%dT%H:%M%z')}",
           f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}",
           f"numpy {np.__version__}", f"pandas {pd.__version__}", f"scipy {scipy.__version__}",
           f"statsmodels {statsmodels.__version__}", f"gemmi {gemmi.__version__}", f"matplotlib {matplotlib.__version__}",
           "nothing installed, upgraded or downgraded by SV01b; no network used"]
    (OUT / "env.txt").write_text("\n".join(env) + "\n")
    print(json.dumps({k: v for k, v in nums.items() if k not in ("mobility_undefined_chains",)}, indent=1))
    print(vt.to_string())


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        main()
