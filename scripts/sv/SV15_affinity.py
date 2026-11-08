"""SV15: affinities across assay systems, relative to E2 in the same assay (RUN_SV1028 SV15 steps 1-6).

Inputs (rule 1.10, cards in tracker/esr1.tasks):
  SV12 card: results/sv/SV12/chembl_activities.csv (CHEMBL206 activities per matched ChEMBL parent, from the cached
             ChEMBL responses in data/sv/chembl_cache/), ligand_ids.csv (CCD -> parent, match type, E2 flag),
             action_class.csv (class per entry x CCD, primary_set), class_coverage.csv (coverage gate).
  SV10 card: results/sv/SV10/ligand_druglikeness.csv (ligand groups, >= 80 % rule; 'compared' = cl7/cl4);
             results/sv/SV10/entry_druglikeness.csv (entry group and entry ligand, same SV10 rules).
  SV14 card: results/sv/SV14/pocket_per_chain.csv (filling %, ligand volume; in_primary_analysis).
No network call: every ChEMBL field needed (assay id, type, standard type/relation, pChEMBL, validity comment,
duplicate flag) is already in chembl_activities.csv.

[DEFAULT] rules applied as written:
  clean: target CHEMBL206; standard_relation '='; pChEMBL present; no data_validity_comment; potential_duplicate 0;
         connectivity_only ligands stay out of the primary set (run file SV12 step 2: used only in a sensitivity run).
  binding = assay_type B with standard_type Ki, Kd or IC50; functional = assay_type F (listed only, never compared);
  assay types A, T, P, U dropped.
  primary: dpX = pX(ligand) - pX(E2) in the same binding assay; RBA = 100 * 10^dpX; per ligand median dpX, n assays.
  fallback: ligands without a same-assay E2 pair get the median binding pChEMBL (Ki/Kd before IC50), flag 'absolute';
            drawn hollow, never tested.
  check: ligands with >= 2 paired assays: SD of dpX vs SD of absolute pX over the same assays.
  comparisons (E2 excluded): dpX by SV10 ligand group (cl7 vs cl4), by class (SV12, coverage gate >= 5 entries per
  side), Spearman rho with filling % and ligand volume (SV14 DONE); ligand level primary, entry level secondary; BH.
Declared before looking at any outcome, where the run file is open:
  - ligand = CCD code (the SV10 / SV14 ligand unit); its ChEMBL data = all CHEMBL206 activities of its ChEMBL parent
    (parent and its forms, as fetched by SV12). E2 = parent CHEMBL135 (CCD EST).
  - assay value = median pChEMBL of the cleaned activities of one parent in one assay with one standard type;
    a pair needs the same assay_chembl_id AND the same standard_type for the ligand and E2 (pairs whose types differ
    inside one assay are not formed and are counted); one dpX per ligand x assay (step 3 'for every binding assay'):
    an assay with two paired types gives the median of its type-wise dpX; ligand dpX = median over its paired
    assays, n assays = number of paired assays.
  - fallback value = median of the cleaned binding activities' pChEMBL: Ki and Kd if the ligand has any, else IC50.
    Plot position only (never tested): fallback value minus the E2 median of the same family (E2 Ki/Kd or E2 IC50,
    all cleaned E2 binding activities), column dpx_plot.
  - check set = primary-set ligands except E2 with >= 2 distinct paired assays; SD with ddof 1 over the paired assay
    values; pooled SD = sqrt(sum((n-1) SD^2) / sum(n-1)). Descriptive only (no test).
  - tested ligands = primary set, basis 'relative', not E2. Ligand level = one value per CCD. Entry level = SV10 entry
    ligand (ligand_id) and SV10 entry group; entry value = its ligand's dpX; class of an entry = SV12 action_class.csv
    row (entry, entry ligand).
  - cluster test: SV10 'compared' ligands / entries, cl7 vs cl4; Mann-Whitney U two-sided, rank-biserial r.
  - class test: classes that pass the SV12 gate (class_coverage.csv: set primary, source any, cluster all,
    entries_ge5) AND keep >= 5 entries among the SV15 tested entries; 2 classes -> Mann-Whitney, >= 3 -> Kruskal-Wallis;
    the ligand-level test also needs >= 2 ligands per class. Otherwise a row with tested False and the reason.
  - Spearman: ligand = median over SV14 primary-analysis chains whose ligand_ids_used equals the CCD; entry = median
    over the entry's SV14 primary-analysis chains whose ligand_ids_used equals the entry ligand.
  - BH over every tested p-value of the primary set together; the sensitivity set (connectivity_only ligands
    included, every CCD a unit even when it shares its parent with another CCD) is a separate BH family.
  - rule 1.4 'per chain' is not run: the step fixes ligand (primary) and entry (secondary) levels and dpX is one value
    per ligand, so a chain level would repeat entry values.
Correction after the first run (disclosed in the SUMMARY): the first run formed one value per ligand x assay x
standard type, so an assay reporting IC50 and Ki for both compounds counted twice; the step text asks for one dpX per
binding assay, which is what this version does. The first-run unit is kept in dpx_median_type_units and its tests in
tests_first_run_type_units.csv (disclosure only, not used in any verdict).
Outputs (results/sv/SV15/): affinity_per_ligand.csv, assay_spread.csv, tests.csv, paired_assays.csv,
entry_level.csv, tests_first_run_type_units.csv, numbers.json.
Run: .venv/bin/python -B scripts/sv/SV15_affinity.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV15"
E2 = "CHEMBL135"
E2_CCD = "EST"
TARGET = "CHEMBL206"
BINDING_TYPES = ("Ki", "Kd", "IC50")
GATE_ENTRIES = 5


# --------------------------------------------------------------------------------------------- cleaning, systems
def clean(a: pd.DataFrame):
    """Sequential [DEFAULT] cleaning; returns the cleaned frame and the cascade counts."""
    steps = [("target CHEMBL206", a.target_chembl_id == TARGET),
             ("standard_relation '='", a.standard_relation == "="),
             ("pChEMBL present", a.pchembl_value.notna()),
             ("no data_validity_comment", a.data_validity_comment.isna()),
             ("not a potential duplicate", a.potential_duplicate.fillna(0).astype(int) == 0)]
    cascade = [{"step": "raw", "kept": int(len(a)), "removed": 0}]
    keep = pd.Series(True, index=a.index)
    for name, m in steps:
        before = int(keep.sum())
        keep &= m
        cascade.append({"step": name, "kept": int(keep.sum()), "removed": before - int(keep.sum())})
    c = a[keep].copy()
    c["system"] = np.select([(c.assay_type == "B") & c.standard_type.isin(BINDING_TYPES), c.assay_type == "F"],
                            ["binding", "functional"], "dropped")
    return c, cascade


def assay_values(b: pd.DataFrame) -> pd.DataFrame:
    """One value per parent x assay x standard type (median pChEMBL of its cleaned activities)."""
    g = b.groupby(["chembl_parent", "assay_chembl_id", "standard_type"])
    v = g.agg(px=("pchembl_value", "median"), n_act=("pchembl_value", "size"),
              document_chembl_id=("document_chembl_id", "first"), document_year=("document_year", "first"),
              assay_description=("assay_description", "first"),
              assay_variant_mutation=("assay_variant_mutation", "first")).reset_index()
    return v


def pairs(v: pd.DataFrame):
    """Same-assay, same-type pairs with E2, then ONE value per ligand x assay (step 3: 'for every binding assay').

    An assay that reports two standard types for both compounds (e.g. IC50 and the Ki derived from it) gives one
    dpX = median of its type-wise dpX (px and px_e2 likewise); the type-wise values are kept in type_detail.
    """
    e2 = v[v.chembl_parent == E2][["assay_chembl_id", "standard_type", "px", "n_act"]].rename(
        columns={"px": "px_e2", "n_act": "n_act_e2"})
    lig = v[v.chembl_parent != E2]
    pt = lig.merge(e2, on=["assay_chembl_id", "standard_type"], how="inner")
    pt["dpx"] = pt.px - pt.px_e2
    # same assay, different standard type -> no pair (counted)
    e2_assays = set(e2.assay_chembl_id)
    keyset = set(zip(pt.chembl_parent, pt.assay_chembl_id, pt.standard_type))
    mism = lig[lig.assay_chembl_id.isin(e2_assays)]
    mism = mism[[(r.chembl_parent, r.assay_chembl_id, r.standard_type) not in keyset for r in mism.itertuples()]]
    pt["detail"] = [f"{r.standard_type} {r.px:.2f}/{r.px_e2:.2f}" for r in pt.itertuples()]
    p = pt.groupby(["chembl_parent", "assay_chembl_id"]).agg(
        standard_types=("standard_type", lambda s: ";".join(sorted(s))), n_types=("standard_type", "size"),
        px=("px", "median"), px_e2=("px_e2", "median"), dpx=("dpx", "median"), n_act=("n_act", "sum"),
        n_act_e2=("n_act_e2", "sum"), type_detail=("detail", "; ".join),
        document_chembl_id=("document_chembl_id", "first"), document_year=("document_year", "first"),
        assay_variant_mutation=("assay_variant_mutation", "first"),
        assay_description=("assay_description", "first")).reset_index()
    p["rba"] = 100 * 10 ** p.dpx
    return p, pt, mism


def fallback_value(b: pd.DataFrame):
    """Median binding pChEMBL over activities: Ki/Kd if any, else IC50."""
    kk = b[b.standard_type.isin(["Ki", "Kd"])]
    if len(kk):
        return float(kk.pchembl_value.median()), "Ki/Kd", int(len(kk))
    ic = b[b.standard_type == "IC50"]
    if len(ic):
        return float(ic.pchembl_value.median()), "IC50", int(len(ic))
    return np.nan, "", 0


def type_summary(x: pd.DataFrame) -> str:
    if not len(x):
        return ""
    parts = []
    for t, s in x.groupby("standard_type"):
        parts.append(f"{t}: {len(s)} (median pX {s.pchembl_value.median():.2f})")
    return "; ".join(parts)


# --------------------------------------------------------------------------------------------- per-ligand table
def ligand_table(a, c, v, p, pt):
    lig = pd.read_csv(ROOT / "results/sv/SV12/ligand_ids.csv")
    s10 = pd.read_csv(ROOT / "results/sv/SV10/ligand_druglikeness.csv")[["ligand_id", "ligand_name", "group",
                                                                        "compared"]]
    ac = pd.read_csv(ROOT / "results/sv/SV12/action_class.csv")
    e10 = pd.read_csv(ROOT / "results/sv/SV10/entry_druglikeness.csv")
    p14 = pd.read_csv(ROOT / "results/sv/SV14/pocket_per_chain.csv")
    p14 = p14[p14.in_primary_analysis.astype(bool)]
    cls = ac.groupby("ccd").agg(action_class=("action_class", "first"), class_source=("source", "first"),
                                n_class_values=("action_class", "nunique"))
    assert cls.n_class_values.max() == 1, "a CCD carries several classes in action_class.csv"
    e2bind = c[(c.chembl_parent == E2) & (c.system == "binding")]
    e2ref = {"Ki/Kd": float(e2bind[e2bind.standard_type.isin(["Ki", "Kd"])].pchembl_value.median()),
             "IC50": float(e2bind[e2bind.standard_type == "IC50"].pchembl_value.median())}
    rows = []
    for r in lig.merge(s10, left_on="ccd", right_on="ligand_id", how="left").itertuples():
        par = r.chembl_parent if isinstance(r.chembl_parent, str) else None
        row = {"ccd": r.ccd, "ligand_name": r.ligand_name, "sv10_group": r.group, "sv10_compared": bool(r.compared),
               "chembl_parent": par, "match_type": r.match_type,
               "primary_set": r.match_type != "connectivity_only",
               "action_class": cls.action_class.get(r.ccd, "unlabelled"),
               "class_source": cls.class_source.get(r.ccd, "none"),
               "is_e2_reference": r.ccd == E2_CCD}
        aa = a[a.chembl_parent == par] if par else a.iloc[0:0]
        cc = c[c.chembl_parent == par] if par else c.iloc[0:0]
        bb = cc[cc.system == "binding"]
        ff = cc[cc.system == "functional"]
        pp = p[p.chembl_parent == par] if par else p.iloc[0:0]
        ppt = pt[pt.chembl_parent == par] if par else pt.iloc[0:0]
        row.update({"n_act_chembl206": int(len(aa)), "n_act_clean": int(len(cc)),
                    "n_binding_act": int(len(bb)), "n_binding_assays": int(bb.assay_chembl_id.nunique()),
                    "binding_types": type_summary(bb),
                    "n_dropped_system_act": int((cc.system == "dropped").sum()),
                    "n_paired_assays": int(len(pp)),
                    "n_paired_two_type_assays": int((pp.n_types > 1).sum()),
                    # first-run unit (assay x standard type), kept for transparency; not used downstream
                    "dpx_median_type_units": float(ppt.dpx.median()) if len(ppt) else np.nan})
        if row["is_e2_reference"]:
            row.update({"dpx_median": 0.0, "dpx_q1": 0.0, "dpx_q3": 0.0, "dpx_min": 0.0, "dpx_max": 0.0,
                        "rba_median": 100.0, "affinity_basis": "reference (E2)",
                        "n_paired_assays": int(bb.assay_chembl_id.nunique())})
        elif len(pp):
            d = pp.dpx
            row.update({"dpx_median": float(d.median()), "dpx_q1": float(d.quantile(.25)),
                        "dpx_q3": float(d.quantile(.75)), "dpx_min": float(d.min()), "dpx_max": float(d.max()),
                        "rba_median": float(100 * 10 ** d.median()), "affinity_basis": "relative"})
        else:
            basis = ("no_chembl_parent" if par is None else
                     "no_binding_data" if not len(bb) else "absolute")
            row.update({"dpx_median": np.nan, "dpx_q1": np.nan, "dpx_q3": np.nan, "dpx_min": np.nan,
                        "dpx_max": np.nan, "rba_median": np.nan, "affinity_basis": basis})
        fv, ft, fn = fallback_value(bb)
        row.update({"px_abs": fv, "px_abs_type": ft, "px_abs_n_act": fn})
        if row["affinity_basis"] == "relative" or row["is_e2_reference"]:
            row["dpx_plot"] = row["dpx_median"]
        elif row["affinity_basis"] == "absolute":
            row["e2_abs_ref"] = e2ref[ft]
            row["dpx_plot"] = fv - e2ref[ft]
        else:
            row["dpx_plot"] = np.nan
        row.update({"n_functional_act": int(len(ff)), "n_functional_assays": int(ff.assay_chembl_id.nunique()),
                    "functional_types": type_summary(ff)})
        s = p14[p14.ligand_ids_used == r.ccd]
        row.update({"sv14_chains": int(len(s)), "sv14_entries": int(s.entry_id.nunique()),
                    "filling_pct_median": float(s.filling_pct.median()) if len(s) else np.nan,
                    "ligand_volume_A3_median": float(s.ligand_volume_A3.median()) if len(s) else np.nan,
                    "n_entries_sv10_entry_ligand": int((e10.ligand_id == r.ccd).sum())})
        rows.append(row)
    t = pd.DataFrame(rows)
    t["in_tests_primary"] = t.primary_set & (t.affinity_basis == "relative") & ~t.is_e2_reference
    t["in_tests_sensitivity"] = (t.affinity_basis == "relative") & ~t.is_e2_reference
    cols = ["ccd", "ligand_name", "sv10_group", "sv10_compared", "chembl_parent", "match_type", "primary_set",
            "action_class", "class_source", "is_e2_reference", "affinity_basis", "dpx_median", "rba_median",
            "n_paired_assays", "n_paired_two_type_assays", "dpx_q1", "dpx_q3", "dpx_min", "dpx_max", "px_abs", "px_abs_type",
            "px_abs_n_act", "e2_abs_ref", "dpx_plot", "n_act_chembl206", "n_act_clean", "n_binding_act",
            "n_binding_assays", "binding_types", "n_dropped_system_act", "n_functional_act", "n_functional_assays",
            "functional_types", "dpx_median_type_units", "sv14_chains", "sv14_entries", "filling_pct_median", "ligand_volume_A3_median",
            "n_entries_sv10_entry_ligand", "in_tests_primary", "in_tests_sensitivity"]
    return t[cols], e2ref


def entry_table(t):
    e10 = pd.read_csv(ROOT / "results/sv/SV10/entry_druglikeness.csv")
    ac = pd.read_csv(ROOT / "results/sv/SV12/action_class.csv")
    p14 = pd.read_csv(ROOT / "results/sv/SV14/pocket_per_chain.csv")
    p14 = p14[p14.in_primary_analysis.astype(bool)]
    e = e10[["entry_id", "group", "compared", "ligand_id", "ligand_tie"]].rename(
        columns={"group": "sv10_entry_group", "compared": "sv10_compared", "ligand_id": "ccd"})
    e = e.merge(ac[["entry_id", "ccd", "action_class", "primary_set"]], on=["entry_id", "ccd"], how="left")
    e = e.merge(t[["ccd", "chembl_parent", "match_type", "affinity_basis", "dpx_median", "dpx_median_type_units",
                   "dpx_plot",
                   "n_paired_assays", "is_e2_reference", "in_tests_primary", "in_tests_sensitivity"]],
                on="ccd", how="left")
    vals = []
    for r in e.itertuples():
        s = p14[(p14.entry_id == r.entry_id) & (p14.ligand_ids_used == r.ccd)]
        vals.append((len(s), s.filling_pct.median() if len(s) else np.nan,
                     s.ligand_volume_A3.median() if len(s) else np.nan))
    e[["sv14_chains", "filling_pct_median", "ligand_volume_A3_median"]] = pd.DataFrame(vals, index=e.index)
    for col in ("in_tests_primary", "in_tests_sensitivity", "is_e2_reference"):
        e[col] = e[col].fillna(False).astype(bool)
    e["action_class"] = e.action_class.fillna("unlabelled")
    return e


# --------------------------------------------------------------------------------------------- tests
def mwu_row(x1, x2, g1, g2, **meta):
    base = {**meta, "groups": f"{g1} vs {g2}", "n_1": len(x1), "n_2": len(x2),
            "median_1": float(np.median(x1)) if len(x1) else np.nan,
            "median_2": float(np.median(x2)) if len(x2) else np.nan}
    if len(x1) < 1 or len(x2) < 1:
        return {**base, "tested": False, "note": "empty group"}
    u = mannwhitneyu(x1, x2, alternative="two-sided")
    rb = 2 * float(u.statistic) / (len(x1) * len(x2)) - 1
    return {**base, "test": "Mann-Whitney U (two-sided)", "statistic": float(u.statistic), "effect": rb,
            "effect_type": f"rank-biserial r (> 0: {g1} higher)",
            "median_difference": base["median_1"] - base["median_2"], "p": float(u.pvalue), "tested": True,
            "note": ""}


def spearman_row(df, x, **meta):
    df = df[["dpx", x]].dropna()
    base = {**meta, "groups": f"dpX vs {x}", "n_1": int(len(df))}
    if len(df) < 4 or df[x].nunique() < 2 or df.dpx.nunique() < 2:
        return {**base, "tested": False, "note": "fewer than 4 units or no variation"}
    rho, p = spearmanr(df.dpx, df[x])
    return {**base, "test": "Spearman", "statistic": float(rho), "effect": float(rho), "effect_type": "Spearman rho",
            "p": float(p), "tested": True, "note": ""}


def run_tests(t, e, gate12, set_name, dpx_col="dpx_median"):
    flag = "in_tests_sensitivity" if set_name.startswith("sensitivity") else "in_tests_primary"
    L = t[t[flag]].copy()
    L["dpx"] = L[dpx_col]
    E = e[e[flag]].copy()
    E["dpx"] = E[dpx_col]
    rows = []
    meta = {"set": set_name}
    # cluster
    lc = L[L.sv10_compared & L.sv10_group.isin(["cl7", "cl4"])]
    rows.append(mwu_row(lc[lc.sv10_group == "cl7"].dpx.to_numpy(), lc[lc.sv10_group == "cl4"].dpx.to_numpy(),
                        "cl7", "cl4", **meta, comparison="dpX by SV10 ligand group", level="ligand"))
    ec = E[E.sv10_compared.astype(bool) & E.sv10_entry_group.isin(["cl7", "cl4"])]
    rows.append(mwu_row(ec[ec.sv10_entry_group == "cl7"].dpx.to_numpy(),
                        ec[ec.sv10_entry_group == "cl4"].dpx.to_numpy(),
                        "cl7", "cl4", **meta, comparison="dpX by SV10 ligand group", level="entry"))
    # class
    n_ent = {c: int((E.action_class == c).sum()) for c in gate12}
    n_lig = {c: int((L.action_class == c).sum()) for c in gate12}
    gated = [c for c in gate12 if n_ent[c] >= GATE_ENTRIES]
    gate_info = {"sv12_gate_classes": gate12, "entries_in_sv15_tested_set": n_ent,
                 "ligands_in_sv15_tested_set": n_lig, "classes_passing": gated}
    if len(gated) >= 2:
        for level, D in (("ligand", L), ("entry", E)):
            vals = [D[D.action_class == c].dpx.to_numpy() for c in gated]
            if level == "ligand" and min(len(v) for v in vals) < 2:
                rows.append({**meta, "comparison": "dpX by action class", "level": level,
                             "groups": " vs ".join(gated), "tested": False,
                             "note": "fewer than 2 ligands in a class: " + ", ".join(
                                 f"{c} {len(v)}" for c, v in zip(gated, vals))})
                continue
            if len(gated) == 2:
                rows.append(mwu_row(vals[0], vals[1], gated[0], gated[1], **meta, comparison="dpX by action class",
                                    level=level))
            else:
                h = kruskal(*vals)
                rows.append({**meta, "comparison": "dpX by action class", "level": level,
                             "groups": " vs ".join(gated), "test": "Kruskal-Wallis", "statistic": float(h.statistic),
                             "p": float(h.pvalue), "tested": True, "note": ""})
    else:
        why = "; ".join(f"{c}: {n_ent[c]} entries, {n_lig[c]} ligands" for c in gate12)
        for level in ("ligand", "entry"):
            rows.append({**meta, "comparison": "dpX by action class", "level": level,
                         "groups": " vs ".join(gate12), "tested": False,
                         "note": f"coverage gate: fewer than 2 classes with >= {GATE_ENTRIES} entries in the SV15 "
                                 f"tested set ({why})"})
    # Spearman with SV14
    for x in ("filling_pct_median", "ligand_volume_A3_median"):
        name = "filling %" if x.startswith("filling") else "ligand volume"
        rows.append(spearman_row(L, x, **meta, comparison=f"dpX vs {name} (SV14)", level="ligand"))
        rows.append(spearman_row(E, x, **meta, comparison=f"dpX vs {name} (SV14)", level="entry"))
    df = pd.DataFrame(rows)
    df["tested"] = df.tested.astype(bool)
    df["q_bh"] = np.nan
    m = df.tested
    if m.any():
        df.loc[m, "q_bh"] = multipletests(df.loc[m, "p"], method="fdr_bh")[1]
    return df, gate_info


# --------------------------------------------------------------------------------------------- spread check
def spread(t, p):
    keep = t[t.in_tests_primary & (t.n_paired_assays >= 2)]
    rows = []
    for r in keep.itertuples():
        pp = p[p.chembl_parent == r.chembl_parent]
        rows.append({"ccd": r.ccd, "chembl_parent": r.chembl_parent, "sv10_group": r.sv10_group,
                     "action_class": r.action_class, "n_paired_assays": int(len(pp)),
                     "standard_types": " | ".join(pp.standard_types),
                     "sd_px_abs": float(pp.px.std(ddof=1)), "sd_dpx": float(pp.dpx.std(ddof=1)),
                     "sd_px_e2": float(pp.px_e2.std(ddof=1)), "range_px_abs": float(pp.px.max() - pp.px.min()),
                     "range_dpx": float(pp.dpx.max() - pp.dpx.min()),
                     "corr_px_lig_e2": float(np.corrcoef(pp.px, pp.px_e2)[0, 1]) if len(pp) >= 3
                     and pp.px.nunique() > 1 and pp.px_e2.nunique() > 1 else np.nan})
    s = pd.DataFrame(rows)
    s["sd_reduction"] = s.sd_px_abs - s.sd_dpx
    s["sd_ratio_dpx_over_abs"] = s.sd_dpx / s.sd_px_abs
    s["dpx_spread_smaller"] = s.sd_dpx < s.sd_px_abs
    return s


def pooled(sd, n):
    df = n - 1
    return float(np.sqrt((df * sd ** 2).sum() / df.sum()))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    a = pd.read_csv(ROOT / "results/sv/SV12/chembl_activities.csv")
    c, cascade = clean(a)
    b = c[c.system == "binding"]
    v = assay_values(b)
    p, pt, mism = pairs(v)
    t, e2ref = ligand_table(a, c, v, p, pt)
    e = entry_table(t)
    cov = pd.read_csv(ROOT / "results/sv/SV12/class_coverage.csv")
    gate12 = cov[(cov.set == "primary") & (cov.source == "any") & (cov.cluster == "all")
                 & (cov.entries_ge5.astype(str) == "True")].action_class.tolist()
    tp, gate_p = run_tests(t, e, gate12, "primary")
    ts, gate_s = run_tests(t, e, gate12, "sensitivity (connectivity_only included)")
    tests = pd.concat([tp, ts], ignore_index=True)
    order = ["set", "comparison", "level", "groups", "tested", "n_1", "n_2", "median_1", "median_2",
             "median_difference", "test", "statistic", "effect", "effect_type", "p", "q_bh", "note"]
    tests = tests[[c_ for c_ in order if c_ in tests.columns]]
    # disclosure only: the first run paired per assay x standard type (an assay with IC50 and Ki counted twice);
    # corrected to one dpX per assay as step 3 is written. Not used in any verdict.
    tf, _ = run_tests(t, e, gate12, "primary, first-run unit assay x standard type (disclosure only)",
                      dpx_col="dpx_median_type_units")
    tf[[c_ for c_ in order if c_ in tf.columns]].to_csv(OUT / "tests_first_run_type_units.csv", index=False)
    sp = spread(t, p)
    # paired assays file (supporting): ligand CCDs per parent
    p2ccd = t.dropna(subset=["chembl_parent"]).groupby("chembl_parent").ccd.apply(lambda s: ";".join(sorted(s)))
    p2prim = t.dropna(subset=["chembl_parent"]).groupby("chembl_parent").primary_set.any()
    pa = p.assign(ccd_codes=p.chembl_parent.map(p2ccd), parent_in_primary_set=p.chembl_parent.map(p2prim))
    pa = pa[["chembl_parent", "ccd_codes", "parent_in_primary_set", "assay_chembl_id", "standard_types", "n_types",
             "px", "px_e2", "dpx", "rba", "type_detail", "n_act", "n_act_e2", "document_chembl_id", "document_year",
             "assay_variant_mutation", "assay_description"]].sort_values(["chembl_parent", "assay_chembl_id"])
    t.to_csv(OUT / "affinity_per_ligand.csv", index=False)
    e.to_csv(OUT / "entry_level.csv", index=False)
    pa.to_csv(OUT / "paired_assays.csv", index=False)
    tests.to_csv(OUT / "tests.csv", index=False)
    sp.to_csv(OUT / "assay_spread.csv", index=False)

    e2v = v[v.chembl_parent == E2]
    tested = t[t.in_tests_primary]
    numbers = {
        "activities_raw": int(len(a)), "cleaning_cascade": cascade,
        "clean_by_assay_type": c.assay_type.value_counts().to_dict(),
        "clean_by_system": c.system.value_counts().to_dict(),
        "clean_B_not_Ki_Kd_IC50": int(((c.assay_type == "B") & ~c.standard_type.isin(BINDING_TYPES)).sum()),
        "clean_B_other_types": c[(c.assay_type == "B") & ~c.standard_type.isin(BINDING_TYPES)]
        .standard_type.value_counts().to_dict(),
        "binding_act_by_type": b.standard_type.value_counts().to_dict(),
        "binding_parents": int(b.chembl_parent.nunique()),
        "binding_assay_values": int(len(v)), "binding_assays": int(v.assay_chembl_id.nunique()),
        "e2_binding_act": int((b.chembl_parent == E2).sum()), "e2_binding_assays": int(e2v.assay_chembl_id.nunique()),
        "e2_binding_assay_values_by_type": e2v.standard_type.value_counts().to_dict(),
        "e2_abs_reference_medians": e2ref,
        "pairs": int(len(p)), "paired_assays": int(p.assay_chembl_id.nunique()),
        "paired_parents": int(p.chembl_parent.nunique()),
        "type_level_pairs": int(len(pt)), "type_level_pairs_by_type": pt.standard_type.value_counts().to_dict(),
        "pairs_two_types": int((p.n_types > 1).sum()),
        "two_type_assays": sorted(p[p.n_types > 1].assay_chembl_id.unique().tolist()),
        "two_type_max_abs_dpx_difference": float(pt.groupby(["chembl_parent", "assay_chembl_id"]).dpx
                                                 .agg(lambda s: s.max() - s.min()).max()),
        "pairs_variant_assays": int(p.assay_variant_mutation.notna().sum()),
        "unpaired_type_mismatch_values": int(len(mism)),
        "unpaired_type_mismatch_parents": sorted(mism.chembl_parent.unique().tolist()),
        "ligands": int(len(t)), "basis_counts": t.affinity_basis.value_counts().to_dict(),
        "basis_counts_primary_set": t[t.primary_set].affinity_basis.value_counts().to_dict(),
        "basis_counts_connectivity_only": t[~t.primary_set].affinity_basis.value_counts().to_dict(),
        "absolute_by_type": t[t.affinity_basis == "absolute"].px_abs_type.value_counts().to_dict(),
        "tested_ligands_primary": int(len(tested)),
        "tested_ligands_primary_by_group": tested.sv10_group.value_counts().to_dict(),
        "tested_ligands_primary_by_class": tested.action_class.value_counts().to_dict(),
        "tested_ligands_sensitivity": int(t.in_tests_sensitivity.sum()),
        "tested_entries_primary": int(e.in_tests_primary.sum()),
        "tested_entries_primary_by_group": e[e.in_tests_primary].sv10_entry_group.value_counts().to_dict(),
        "tested_entries_primary_by_class": e[e.in_tests_primary].action_class.value_counts().to_dict(),
        "tested_entries_sensitivity": int(e.in_tests_sensitivity.sum()),
        "e2_entries_sv10": int((e.ccd == E2_CCD).sum()),
        "dpx_tested_median": float(tested.dpx_median.median()),
        "dpx_tested_q1": float(tested.dpx_median.quantile(.25)), "dpx_tested_q3": float(tested.dpx_median.quantile(.75)),
        "dpx_tested_min": float(tested.dpx_median.min()), "dpx_tested_max": float(tested.dpx_median.max()),
        "tested_ligands_dpx_gt0": int((tested.dpx_median > 0).sum()),
        "tested_paired_assays_median": float(tested.n_paired_assays.median()),
        "tested_single_assay_ligands": int((tested.n_paired_assays == 1).sum()),
        "tested_with_sv14": int(tested.filling_pct_median.notna().sum()),
        "class_gate_primary": gate_p, "class_gate_sensitivity": gate_s,
        "spread_ligands": int(len(sp)),
        "spread_median_sd_abs": float(sp.sd_px_abs.median()) if len(sp) else np.nan,
        "spread_median_sd_dpx": float(sp.sd_dpx.median()) if len(sp) else np.nan,
        "spread_median_sd_e2": float(sp.sd_px_e2.median()) if len(sp) else np.nan,
        "spread_pooled_sd_abs": pooled(sp.sd_px_abs, sp.n_paired_assays) if len(sp) else np.nan,
        "spread_pooled_sd_dpx": pooled(sp.sd_dpx, sp.n_paired_assays) if len(sp) else np.nan,
        "spread_n_smaller": int(sp.dpx_spread_smaller.sum()) if len(sp) else 0,
        "spread_median_ratio": float(sp.sd_ratio_dpx_over_abs.median()) if len(sp) else np.nan,
        "spread_median_reduction": float(sp.sd_reduction.median()) if len(sp) else np.nan,
        "functional_ligands": int((t.n_functional_act > 0).sum()),
        "functional_act": int(t.drop_duplicates("chembl_parent").n_functional_act.sum()),
        "shared_parents": {k: v_ for k, v_ in p2ccd.items() if ";" in v_},
    }
    (OUT / "numbers.json").write_text(json.dumps(numbers, indent=1, default=float))
    print(json.dumps({k: numbers[k] for k in ("cleaning_cascade", "clean_by_system", "pairs", "basis_counts",
                                              "tested_ligands_primary", "class_gate_primary", "spread_ligands")},
                     indent=1, default=float))
    print(tests.to_string())


if __name__ == "__main__":
    main()
