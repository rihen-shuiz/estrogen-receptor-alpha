"""SV13 (RUN_SV1028) Parts A and B: action classes by cluster, per-chain outcome table, coverage gate and the
pre-fixed class tests (C1 destabilisation, C2a agonist vs partial agonist, C2b agonists vs dpX).

Inputs (rule 1.10, from the tracker cards' out paths):
  SV12  results/sv/SV12/chain_class.csv (class per chain = class of its SV00 ligand_id_primary; primary_set,
        match_type), action_class.csv, class_coverage.csv
  SV00  results/sv/SV00/chain_table.csv (788 kept chains, predictors)
  SV01a results/sv/SV01a/chain_outcomes.csv (H12 complete 538-547, fraction observed, h12_status)
  SV01b results/sv/SV01b/mobility_per_chain.csv (mobility_z), lattice_contacts_per_chain.csv, model_data.csv (Q and
        X blocks with lattice terms, imputed within cl7 / cl4); model code scripts/sv/SV01b_models.py and
        scripts/sv/SV01a_completeness.py (screen, fit, code_predictors) are imported, not copied
  SV02  results/sv/SV02/phipsi_per_chain.csv (psi at 537)
  SV03  results/sv/SV03/geometry_per_chain.csv (angle_H12_H3_deg, dist_H12c_H11axis_A, ulen_H12 = uncensored H12
        length; H12 window 538-547)
  SV15  results/sv/SV15/affinity_per_ligand.csv (dpx_median, in_tests_primary / in_tests_sensitivity)

Rules fixed before any outcome was compared between classes (where the run file is open, declared here):
  - Sets: primary = SV12 primary_set rows (connectivity_only rows out, SV12 step 2); sensitivity = all rows
    (connectivity_only included, SV12 step 2 'used only in a sensitivity run'). Its own BH family.
  - Coverage gate [DEFAULT] (SV12): a class comparison is tested only if each side has >= 5 entries, counted on the
    chains of the comparison's stratum and set (C1: all chains, cl4; C2: cl7). Entry = PDB entry with >= 1 chain of
    that class in the stratum.
  - C2b is a rank correlation inside one class, so its gate is >= 5 agonist entries with a dpX in cl7. dpX = SV15
    dpx_median of the chain's ligand (SV12 ccd) for the ligands SV15 itself tests (in_tests_primary, or
    in_tests_sensitivity for the sensitivity set): relative same-assay dpX only; E2 (the reference, dpX 0 by
    definition) and absolute-only ligands are excluded, as in SV15 steps 4 and 6.
  - Outcomes: complete (logit), fraction observed (fractional logit), mobility_z (OLS) for C1; mobility_z,
    angle H12-H3, distance H12 centroid-H11 axis, uncensored H12 length, psi 537 as the wrapped difference from the
    cl7 circular mean (over every cl7 chain with psi 537) for C2 (OLS for C2a, Spearman for C2b).
  - Model (C1, C2a): outcome ~ class indicator (1 = SERD / antagonist / partial_agonist, 0 = the reference class
    SERM / SERM / agonist) + Q (resolution, R-free, entry mean B) + X (pH, growth T, ambient T, P21, ESR1 chains in
    ASU, dimer by symmetry, 'not reported' indicators, lattice_asu, lattice_sym) + construct mutations (Y537S, L536S,
    D538G, other) + peptide (chain-level peptide bound). Predictor values: cl7 / cl4 strata = SV01b model_data.csv
    (imputed within the cluster); all-chains stratum = the same code (SV01a code_predictors + SV01b lattice terms,
    median imputation + indicator over all chains of the set). Continuous predictors z-scored on the model rows.
    Drop rule as SV01a/SV01b: constant, separating (logit) or collinear covariates are dropped and listed; joint
    separation drops the covariate with the largest |log OR| (never the class term); if the class term itself is
    constant or separating, the class effect is 'not estimable' (listed). Chain level: entry-clustered sandwich
    errors (GEE independence); entry level: one chain per entry (the entry's first chain in chain-table order among
    the model rows), plain MLE / OLS. Effect = OR (logit) or coefficient (OLS) of the class indicator with Wald 95 % CI.
  - C2b: entry level = Spearman of the entry median outcome with the entry dpX; chain level = Spearman over chains,
    95 % CI and two-sided p from a bootstrap over entries (10000 resamples, seed 42; p = (1 + 2 min(#rho*<=0,
    #rho*>=0)) / (B + 1), capped at 1; the SV14 convention).
  - BH within each comparison (C1 SERD vs SERM, C1 antagonist vs SERM, C2a, C2b), per set.
  - [DEFAULT] verdict per comparison x stratum x outcome: 'linked' if q < 0.05 at chain and entry level, 'not linked'
    otherwise, 'not testable' if the coverage gate fails. Comparison verdict: 'not testable' if the gate fails in
    every stratum; else 'linked' if any tested row is linked, else 'not linked'.
Code check: --code-check fits the C1 models once on seeded random pseudo-labels (cl4 primary entries split at random
into two groups) only to exercise the model path for a later re-run with the teammate labels; it writes only
data/sv/SV13/model_code_check.json (statuses, no estimates) and nothing in results/.

Outputs: results/sv/SV13/class_by_cluster.csv, chain_data.csv (supporting), gate.csv (supporting), tests.csv,
numbers_tests.json.
Run: .venv/bin/python -B scripts/sv/SV13_analysis.py [--code-check]
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
import SV01b_models as Bm  # noqa: E402  (imports SV01a as Bm.A and adds the lattice terms to A.CONTINUOUS)

A = Bm.A
OUT = ROOT / "results/sv/SV13"
OUT.mkdir(parents=True, exist_ok=True)
CLASS_ORDER = ("agonist", "partial_agonist", "SERM", "SERD", "antagonist", "other", "unlabelled")
GATE_MIN = 5
B_BOOT, SEED = 10000, 42
Q_COLS = ["resolution_A", "r_free", "mean_b"]
X_COLS = ["ph", "growth_temp_K", "ambient_temp_K", "p21", "esr1_chains_asu", "dimer_symmetry",
          "ph_missing", "growth_temp_K_missing", "ambient_temp_K_missing", "lattice_asu", "lattice_sym"]
M_COLS = ["mut_Y537S", "mut_L536S", "mut_D538G", "mut_other"]
P_COLS = ["peptide"]
COVARIATES = Q_COLS + X_COLS + M_COLS + P_COLS
CONT = set(A.CONTINUOUS) | {"lattice_asu", "lattice_sym"}

C1_OUT = [("complete", "H12 complete (538-547)", "logit"), ("frac_obs", "H12 fraction observed", "logit"),
          ("mobility_z", "H12 mobility z (528-547)", "linear")]
C2_OUT = [("mobility_z", "H12 mobility z (528-547)"), ("angle_H12_H3_deg", "angle H12-H3 (deg)"),
          ("dist_H12c_H11axis_A", "distance H12 centroid to H11 axis (A)"),
          ("ulen_H12", "uncensored H12 length (residues)"), ("psi537_wrapdiff", "psi 537 minus cl7 circular mean (deg)")]
COMPARISONS = [
    dict(id="C1_SERD_vs_SERM", part="C1", test_class="SERD", ref_class="SERM", strata=("all", "cl4"), kind="model",
         outcomes=[(o, lab, fam) for o, lab, fam in C1_OUT]),
    dict(id="C1_antagonist_vs_SERM", part="C1", test_class="antagonist", ref_class="SERM", strata=("all", "cl4"),
         kind="model", outcomes=[(o, lab, fam) for o, lab, fam in C1_OUT]),
    dict(id="C2a_partial_agonist_vs_agonist", part="C2a", test_class="partial_agonist", ref_class="agonist",
         strata=("cl7",), kind="model", outcomes=[(o, lab, "linear") for o, lab in C2_OUT]),
    dict(id="C2b_agonist_rho_dpX", part="C2b", test_class="agonist", ref_class=None, strata=("cl7",), kind="rho",
         outcomes=[(o, lab, "spearman") for o, lab in C2_OUT]),
]


# ------------------------------------------------------------------------------------------------ inputs
def circ_mean_deg(a):
    r = np.deg2rad(np.asarray(a, float))
    return float(np.rad2deg(np.arctan2(np.sin(r).mean(), np.cos(r).mean())))


def wrap180(x):
    return (np.asarray(x, float) + 180.0) % 360.0 - 180.0


def load_chain_data():
    t = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", low_memory=False)
    cc = pd.read_csv(ROOT / "results/sv/SV12/chain_class.csv")
    co = pd.read_csv(ROOT / "results/sv/SV01a/chain_outcomes.csv")
    mob = pd.read_csv(ROOT / "results/sv/SV01b/mobility_per_chain.csv")
    geo = pd.read_csv(ROOT / "results/sv/SV03/geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    pp = pd.read_csv(ROOT / "results/sv/SV02/phipsi_per_chain.csv", usecols=["chain_key", "position", "psi"])
    aff = pd.read_csv(ROOT / "results/sv/SV15/affinity_per_ligand.csv")
    for x in (cc, co, mob, geo):
        assert len(x) == len(t) and (x.chain_key.values == t.chain_key.values).all()
    d = cc[["chain_key", "entry_id", "s08_cluster", "cluster", "ccd", "ligand_codes", "action_class", "source",
            "match_type", "primary_set"]].copy()
    d["primary_set"] = d.primary_set.astype(str).eq("True")
    d["complete"] = co.complete.values
    d["frac_obs"] = co.frac_obs.values
    d["h12_status"] = co.h12_status.values
    d["n_obs_538_547"] = co.n_obs_window.values
    d["mobility_z"] = mob.mobility_z.values
    for c in ("angle_H12_H3_deg", "dist_H12c_H11axis_A", "ulen_H12"):
        d[c] = geo[c].astype(float).values
    p537 = pp[pp.position == 537].set_index("chain_key").psi
    d["psi537"] = d.chain_key.map(p537).astype(float)
    cl7 = d[(d.cluster == "cl7") & d.psi537.notna()]
    mu7 = circ_mean_deg(cl7.psi537)
    d["psi537_wrapdiff"] = wrap180(d.psi537 - mu7)
    a = aff.set_index("ccd")
    d["dpx"] = d.ccd.map(a.dpx_median)
    d["affinity_basis"] = d.ccd.map(a.affinity_basis)
    d["dpx_in_tests_primary"] = d.ccd.map(a.in_tests_primary).astype(str).eq("True")
    d["dpx_in_tests_sensitivity"] = d.ccd.map(a.in_tests_sensitivity).astype(str).eq("True")
    psi_info = dict(cl7_circular_mean_psi537_deg=mu7, n_cl7_chains_with_psi537=int(len(cl7)),
                    n_cl7_chains=int((d.cluster == "cl7").sum()))
    return t, d, psi_info


def class_by_cluster(d):
    rows = []
    for st, s in (("primary", d[d.primary_set]), ("all", d)):
        for cls in CLASS_ORDER:
            x = s[s.action_class == cls]
            for cl in sorted(d.s08_cluster.unique()):
                y = x[x.s08_cluster == cl]
                rows.append(dict(set=st, action_class=cls, s08_cluster=str(int(cl)),
                                 cluster_group={7: "cl7", 4: "cl4"}.get(int(cl), "minor"), n_chains=len(y),
                                 n_entries=y.entry_id.nunique(), n_ligands=y.ccd.dropna().nunique(),
                                 n_complete_h12=int(y.complete.sum())))
            for g in ("cl7", "cl4", "minor", "all"):
                y = x if g == "all" else x[x.cluster == g]
                rows.append(dict(set=st, action_class=cls, s08_cluster=g if g != "minor" else "minor (1,2,3,5,6)",
                                 cluster_group=g, n_chains=len(y), n_entries=y.entry_id.nunique(),
                                 n_ligands=y.ccd.dropna().nunique(), n_complete_h12=int(y.complete.sum())))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ predictors
def build_predictors(t, rows_idx):
    """SV01a code_predictors + SV01b lattice terms for the chains rows_idx; median imputation + 'not reported'
    indicator (created only if a value is missing) over exactly these chains (SV01a rule)."""
    lat = pd.read_csv(ROOT / "results/sv/SV01b/lattice_contacts_per_chain.csv")
    assert (lat.chain_key.values == t.chain_key.values).all()
    dc = A.code_predictors(t)
    dc["lattice_asu"] = lat.asu_atoms.astype(float).values
    dc["lattice_sym"] = lat.sym_atoms.astype(float).values
    dc = dc.iloc[rows_idx].reset_index(drop=True)
    for v in A.IMPUTE:
        if dc[v].isna().any():
            dc[f"{v}_missing"] = dc[v].isna().astype(float)
        dc[v] = dc[v].fillna(float(dc[v].median()))
    dc.insert(0, "chain_key", t.chain_key.values[rows_idx])
    return dc


def predictor_tables(t, d):
    md = pd.read_csv(ROOT / "results/sv/SV01b/model_data.csv")
    tabs, checks = {}, {}
    for g, c in (("cl7", 7), ("cl4", 4)):
        m = md[md.s08_cluster == c].reset_index(drop=True)
        idx = np.where(t.s08_cluster.values == c)[0]
        rb = build_predictors(t, idx)
        assert (rb.chain_key.values == m.chain_key.values).all()
        cols = [x for x in COVARIATES if x in m.columns]
        diff = max(float(np.nanmax(np.abs(rb[x].values - m[x].values))) for x in cols)
        checks[g] = dict(n_chains=len(m), columns_compared=cols, max_abs_diff_rebuilt_vs_model_data=diff)
        tabs[g] = m[["chain_key"] + cols].copy()
    for st in ("primary", "all"):
        idx = np.where(d.primary_set.values)[0] if st == "primary" else np.arange(len(d))
        tabs[f"all|{st}"] = build_predictors(t, idx)
    return tabs, checks


# ------------------------------------------------------------------------------------------------ gate
def stratum_rows(d, stratum, setname):
    s = d if setname == "all" else d[d.primary_set]
    return s if stratum == "all" else s[s.cluster == stratum]


def gate_rows(d):
    rows = []
    for setname in ("primary", "all"):
        for comp in COMPARISONS:
            for stratum in comp["strata"]:
                s = stratum_rows(d, stratum, setname)
                if comp["kind"] == "model":
                    a = s[s.action_class == comp["test_class"]]
                    b = s[s.action_class == comp["ref_class"]]
                    na, nb = a.entry_id.nunique(), b.entry_id.nunique()
                    ok = na >= GATE_MIN and nb >= GATE_MIN
                    why = "" if ok else "; ".join(f"{c} {n} entries < {GATE_MIN}" for c, n in
                                                   ((comp["test_class"], na), (comp["ref_class"], nb)) if n < GATE_MIN)
                    rows.append(dict(set=setname, comparison=comp["id"], part=comp["part"], stratum=stratum,
                                     side_a=comp["test_class"], chains_a=len(a), entries_a=na,
                                     side_b=comp["ref_class"], chains_b=len(b), entries_b=nb, gate_pass=ok,
                                     gate_reason=why))
                else:
                    flag = "dpx_in_tests_primary" if setname == "primary" else "dpx_in_tests_sensitivity"
                    a = s[(s.action_class == comp["test_class"]) & s[flag]]
                    allag = s[s.action_class == comp["test_class"]]
                    na = a.entry_id.nunique()
                    ok = na >= GATE_MIN
                    lig = a.groupby("ccd").entry_id.nunique().to_dict()
                    why = "" if ok else (f"{comp['test_class']} entries with a tested SV15 dpX: {na} < {GATE_MIN} "
                                         f"({', '.join(f'{k} {v}' for k, v in sorted(lig.items()))}; "
                                         f"{allag.entry_id.nunique()} agonist entries in {stratum})")
                    rows.append(dict(set=setname, comparison=comp["id"], part=comp["part"], stratum=stratum,
                                     side_a=f"{comp['test_class']} with dpX", chains_a=len(a), entries_a=na,
                                     side_b="", chains_b=np.nan, entries_b=np.nan, gate_pass=ok, gate_reason=why,
                                     dpx_ligands=";".join(f"{k}:{v}" for k, v in sorted(lig.items()))))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ models
def fit_class_model(rows, X, y, fam, level):
    """rows: chain rows (entry_id), X: predictor frame with 'class_ind' first, y: outcome. Returns a result dict."""
    if level == "entry":
        keep = ~rows.entry_id.duplicated().values
        rows, X, y = rows[keep].reset_index(drop=True), X[keep].reset_index(drop=True), y[keep]
    X = X.copy()
    for c in X.columns:
        if c in CONT:
            mu, sd = float(X[c].mean()), float(X[c].std(ddof=0))
            X[c] = (X[c] - mu) / sd if sd > 0 else 0.0
    preds = list(X.columns)
    groups = rows.entry_id.values if level == "chain" else None
    out = dict(n_chains_model=int(len(y)), n_entries_model=int(rows.entry_id.nunique()))
    if fam == "logit":
        kept, dropped = A.screen(X, y, preds)
        if "class_ind" not in kept:
            return dict(out, status="class term not estimable: " + dropped.get("class_ind", "?"),
                        dropped=";".join(f"{k}: {v}" for k, v in dropped.items() if k != "class_ind"))
        gp = pd.factorize(groups)[0] if groups is not None else None
        while True:
            full = A.fit(y, X[kept], gp)
            if not full["sep"] or full["res"] is None:
                break
            co = full["res"].params.drop(["const", "class_ind"]).abs()
            if not len(co):
                return dict(out, status=f"class term not estimable: joint separation ({full['note']})",
                            dropped=";".join(f"{k}: {v}" for k, v in dropped.items()))
            worst = co.idxmax()
            dropped[worst] = f"joint separation ({full['note']}; |log OR| {co[worst]:.1f})"
            kept.remove(worst)
        res = full["res"]
        if res is None:
            return dict(out, status=f"fit failed ({full['note']})", dropped="")
        b, se = float(res.params["class_ind"]), float(res.bse["class_ind"])
        est = dict(effect_type="odds ratio (class vs reference)", estimate=np.exp(b), ci_low=np.exp(b - 1.96 * se),
                   ci_high=np.exp(b + 1.96 * se), coef=b, se=se, p=float(res.pvalues["class_ind"]))
    else:
        kept, dropped = Bm.screen_linear(X, preds)
        if "class_ind" not in kept:
            return dict(out, status="class term not estimable: " + dropped.get("class_ind", "?"),
                        dropped=";".join(f"{k}: {v}" for k, v in dropped.items() if k != "class_ind"))
        res = Bm.fit_ols(y, X[kept], groups)["res_any"]
        b, se = float(res.params["class_ind"]), float(res.bse["class_ind"])
        est = dict(effect_type="coefficient (class vs reference, outcome units)", estimate=b, ci_low=b - 1.96 * se,
                   ci_high=b + 1.96 * se, coef=b, se=se, p=float(2 * norm.sf(abs(b / se))))
    return dict(out, **est, status="fitted", n_covariates=len(kept) - 1,
                covariates_kept=";".join(k for k in kept if k != "class_ind"),
                dropped=";".join(f"{k}: {v}" for k, v in dropped.items()))


def run_model_comparison(d, tabs, comp, stratum, setname, outcome, fam, level, labels=None):
    s = stratum_rows(d, stratum, setname)
    lab = labels if labels is not None else s.action_class
    s = s[lab.loc[s.index].isin([comp["test_class"], comp["ref_class"]])]
    s = s[s[outcome].notna()]
    tab = tabs[stratum] if stratum in ("cl7", "cl4") else tabs[f"all|{setname}"]
    X = tab.set_index("chain_key").loc[s.chain_key, [c for c in COVARIATES if c in tab.columns]].reset_index(drop=True)
    X.insert(0, "class_ind", (lab.loc[s.index] == comp["test_class"]).astype(float).values)
    y = s[outcome].values.astype(float)
    return fit_class_model(s.reset_index(drop=True), X, y, fam, level)


# ------------------------------------------------------------------------------------------------ C2b
def boot_rho(s, x, y, rng):
    ents = s.entry_id.unique()
    by = {e: g[[x, y]].to_numpy(float) for e, g in s.groupby("entry_id")}
    vals = []
    for _ in range(B_BOOT):
        pick = rng.choice(ents, len(ents), replace=True)
        m = np.vstack([by[e] for e in pick])
        if np.unique(m[:, 0]).size < 2 or np.unique(m[:, 1]).size < 2:
            vals.append(np.nan)
            continue
        vals.append(spearmanr(m[:, 0], m[:, 1]).statistic)
    v = np.asarray(vals)
    v = v[np.isfinite(v)]
    p = min(1.0, (1 + 2 * min((v <= 0).sum(), (v >= 0).sum())) / (len(v) + 1))
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)), p, int(len(v))


def run_rho(d, comp, stratum, setname, outcome, level, rng):
    flag = "dpx_in_tests_primary" if setname == "primary" else "dpx_in_tests_sensitivity"
    s = stratum_rows(d, stratum, setname)
    s = s[(s.action_class == comp["test_class"]) & s[flag] & s[outcome].notna()]
    out = dict(n_chains_model=int(len(s)), n_entries_model=int(s.entry_id.nunique()),
               effect_type="Spearman rho (outcome vs dpX)")
    if level == "entry":
        e = s.groupby("entry_id").agg(v=(outcome, "median"), dpx=("dpx", "first"))
        r = spearmanr(e.v, e.dpx)
        return dict(out, status="fitted", estimate=float(r.statistic), p=float(r.pvalue), ci_low=np.nan,
                    ci_high=np.nan)
    rho = spearmanr(s[outcome], s.dpx).statistic
    lo, hi, p, nb = boot_rho(s, outcome, "dpx", rng)
    return dict(out, status="fitted", estimate=float(rho), ci_low=lo, ci_high=hi, p=p, n_boot_valid=nb)


# ------------------------------------------------------------------------------------------------ tests
def run_tests(d, tabs, gate):
    rng = np.random.default_rng(SEED)
    rows = []
    for setname in ("primary", "all"):
        for comp in COMPARISONS:
            for stratum in comp["strata"]:
                g = gate[(gate.set == setname) & (gate.comparison == comp["id"]) & (gate.stratum == stratum)].iloc[0]
                for outcome, olab, fam in comp["outcomes"]:
                    for level in ("chain", "entry"):
                        base = dict(set=setname, comparison=comp["id"], part=comp["part"], stratum=stratum,
                                    outcome=outcome, outcome_label=olab, family=fam, level=level,
                                    test_class=comp["test_class"], reference_class=comp["ref_class"] or "",
                                    chains_a=g.chains_a, entries_a=g.entries_a, chains_b=g.chains_b,
                                    entries_b=g.entries_b, gate_pass=bool(g.gate_pass))
                        if not g.gate_pass:
                            rows.append(dict(base, status="not testable (coverage gate)", gate_reason=g.gate_reason))
                            continue
                        if comp["kind"] == "model":
                            r = run_model_comparison(d, tabs, comp, stratum, setname, outcome, fam, level)
                        else:
                            r = run_rho(d, comp, stratum, setname, outcome, level, rng)
                        rows.append(dict(base, **r))
    t = pd.DataFrame(rows)
    for c in ("estimate", "ci_low", "ci_high", "p", "q"):
        if c not in t:
            t[c] = np.nan
    for (setname, comp), idx in t.groupby(["set", "comparison"]).groups.items():
        sub = t.loc[idx]
        ok = sub.p.notna()
        if ok.any():
            t.loc[sub.index[ok], "q"] = multipletests(sub.p[ok].values, method="fdr_bh")[1]
    # verdicts
    t["row_verdict"] = ""
    for key, sub in t.groupby(["set", "comparison", "stratum", "outcome"]):
        if not sub.gate_pass.iloc[0]:
            v = "not testable"
        else:
            qs = sub.set_index("level").q
            v = "linked" if all(qs.get(lv, np.nan) < 0.05 for lv in ("chain", "entry")) else "not linked"
        t.loc[sub.index, "row_verdict"] = v
    comp_v = {}
    for (setname, comp), sub in t.groupby(["set", "comparison"]):
        vs = set(sub.row_verdict)
        comp_v[f"{setname}|{comp}"] = ("not testable" if vs == {"not testable"} else
                                       "linked" if "linked" in vs else "not linked")
    t["comparison_verdict"] = [comp_v[f"{a}|{b}"] for a, b in zip(t.set, t.comparison)]
    return t, comp_v


def code_check(d, tabs):
    """Exercise the model path on seeded random pseudo-labels (no estimates saved)."""
    rng = np.random.default_rng(SEED)
    s = d[d.primary_set & (d.cluster == "cl4")]
    ents = s.entry_id.unique()
    pick = set(rng.choice(ents, len(ents) // 3, replace=False))
    lab = pd.Series(np.where(s.entry_id.isin(pick), "SERD", "SERM"), index=s.index)
    lab = lab.reindex(d.index).fillna("unlabelled")
    comp = COMPARISONS[0]
    out = []
    for stratum in ("cl4", "all"):
        for outcome, _, fam in comp["outcomes"]:
            for level in ("chain", "entry"):
                r = run_model_comparison(d, tabs, comp, stratum, "primary", outcome, fam, level, labels=lab)
                out.append(dict(stratum=stratum, outcome=outcome, level=level, status=r["status"],
                                n_chains=r.get("n_chains_model"), n_covariates=r.get("n_covariates"),
                                dropped=r.get("dropped", ""), p_finite=bool(np.isfinite(r.get("p", np.nan)))))
    # C2b path: seeded random pseudo-dpX per ligand for every cl7 agonist chain
    dz = d.copy()
    ligs = sorted(dz.loc[dz.action_class == "agonist", "ccd"].dropna().unique())
    pdx = dict(zip(ligs, rng.normal(0, 1, len(ligs))))
    dz["dpx"] = dz.ccd.map(pdx)
    dz["dpx_in_tests_primary"] = dz.action_class.eq("agonist") & dz.dpx.notna()
    comp2 = COMPARISONS[3]
    rng2 = np.random.default_rng(SEED)
    for outcome, _, _ in comp2["outcomes"]:
        for level in ("chain", "entry"):
            r = run_rho(dz, comp2, "cl7", "primary", outcome, level, rng2)
            out.append(dict(stratum="cl7", outcome=outcome, level=level, status=r["status"] + " (C2b path)",
                            n_chains=r.get("n_chains_model"), n_covariates=None, dropped="",
                            p_finite=bool(np.isfinite(r.get("p", np.nan)))))
    dd = ROOT / "data/sv/SV13"
    dd.mkdir(parents=True, exist_ok=True)
    (dd / "model_code_check.json").write_text(json.dumps(dict(
        note="code-path check only: seeded random pseudo-labels (1/3 of cl4 primary entries as 'SERD', rest "
             "'SERM') for the C1 models and seeded random pseudo-dpX per agonist ligand for the C2b path; no "
             "estimate is saved or reported; not a result", seed=SEED, rows=out), indent=1))
    print(pd.DataFrame(out).to_string())


def main():
    t, d, psi_info = load_chain_data()
    d.to_csv(OUT / "chain_data.csv", index=False)
    cbc = class_by_cluster(d)
    cbc.to_csv(OUT / "class_by_cluster.csv", index=False)
    # cross-check with SV12 class_coverage (primary/all x class x cl7/cl4/minor/all, source any)
    cov = pd.read_csv(ROOT / "results/sv/SV12/class_coverage.csv")
    cov = cov[cov.source == "any"]
    m = cbc[cbc.s08_cluster.isin(["cl7", "cl4", "all"]) | cbc.s08_cluster.str.startswith("minor")].copy()
    m["cl"] = m.cluster_group
    j = m.merge(cov, left_on=["set", "action_class", "cl"], right_on=["set", "action_class", "cluster"],
                suffixes=("", "_sv12"))
    cov_check = dict(rows_compared=int(len(j)), chains_equal=int((j.n_chains == j.n_chains_sv12).sum()),
                     entries_equal=int((j.n_entries == j.n_entries_sv12).sum()))
    tabs, pchecks = predictor_tables(t, d)
    gate = gate_rows(d)
    gate.to_csv(OUT / "gate.csv", index=False)
    if "--code-check" in sys.argv:
        code_check(d, tabs)
    tests, comp_v = run_tests(d, tabs, gate)
    tests.to_csv(OUT / "tests.csv", index=False)
    nums = dict(psi537=psi_info, coverage_crosscheck_vs_SV12=cov_check, predictor_rebuild_check=pchecks,
                gate=gate.to_dict(orient="records"), comparison_verdicts=comp_v,
                n_tests_rows=int(len(tests)), n_fitted=int((tests.get("status") == "fitted").sum()),
                n_not_testable=int(tests.status.str.startswith("not testable").sum()),
                c2b_agonist_chains_with_dpx_cl7=d[(d.cluster == "cl7") & (d.action_class == "agonist")
                                                  & d.dpx_in_tests_primary & d.primary_set][
                    ["chain_key", "entry_id", "ccd", "dpx"]].to_dict(orient="records"),
                gate_min_entries=GATE_MIN)
    (OUT / "numbers_tests.json").write_text(json.dumps(nums, indent=1, default=str))
    print(json.dumps({k: v for k, v in nums.items() if k != "gate"}, indent=1, default=str))
    print(gate.to_string())


if __name__ == "__main__":
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        main()
