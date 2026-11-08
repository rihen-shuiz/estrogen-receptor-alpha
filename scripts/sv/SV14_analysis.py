"""SV14 stage 3: references, per-chain table, correlations and filling tests (RUN_SV1028 SV14 steps 3-6).

Inputs: data/sv/SV14/grid_per_chain.csv, chains.csv, ligand_instances.csv, apo_check.csv, e2_centroid.json (SV14
stages 1-2); results/sv/SV00/chain_table.csv (cluster, construct_mutations, ligand_chemotype); SV12 card:
results/sv/SV12/action_class.csv (class per entry x CCD, primary_set), class_coverage.csv (coverage gate),
chain_class.csv (cross-check only).

[DEFAULT] as written: wild type = construct_mutations 'none' (SV00); ligand 100 % = median EST ligand volume over
wild-type EST chains; apo entries = ligand_chemotype none and the apo check passed by every kept chain of the entry;
pocket 100 % = median of entry medians over apo entries (IQR of the entry medians); if fewer than 5 apo entries pass,
pocket 100 % = median pocket of wild-type EST chains; the E2-bound pocket (median over wild-type EST chains) is
always reported. ligand_incomplete chains stay out of the primary analysis.
Declared before looking at the outcomes, where the run file is open:
  - the E2 reference chains are the wild-type chains whose analysed ligand is EST alone with a complete ligand;
  - primary analysis set = chains with a ligand of interest, ligand complete, pocket found;
  - class of a chain = SV12 action_class.csv row (entry, analysed CCD); class analyses use primary_set rows only;
    a class is analysed if it passes the SV12 coverage gate (class_coverage.csv: set primary, source any, cluster
    all, entries_ge5) and still has >= 5 entries in the SV14 primary set;
  - subsets (cl7, cl4, class) are formed on chains, then medians per ligand (analysed CCD) or per entry;
  - filling tests: cl7 vs cl4 and between gated classes; entry level = Mann-Whitney U (two-sided; Kruskal-Wallis if
    >= 3 groups) on entry medians, entry cluster group = cluster of >= 80 % of its analysed chains (SV10/SV11 rule),
    entry value = median over its chains in the group; chain level = bootstrap over entries (10000 resamples, seed
    42, entries resampled with replacement within each group, statistic = difference of chain medians, two-sided
    p = (1 + 2 min(#<=0, #>=0)) / (B + 1), capped at 1);
  - BH over every p-value of correlations.csv (rows with in_bh True) and tests.csv together; the buriedness
    sensitivity rho values are descriptive (in_bh False).
Outputs: results/sv/SV14/pocket_per_chain.csv, references.json, correlations.csv, tests.csv, filling_by_group.csv,
filling_entry_level.csv (figure input), numbers.json
Run: .venv/bin/python -B scripts/sv/SV14_analysis.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/sv/SV14"
OUT = ROOT / "results/sv/SV14"
SEED = 42
B = 10000
ENTRY_SHARE_MIN = 0.8
MIN_APO_ENTRIES = 5
GATE_ENTRIES = 5
VOXEL = 0.064


def cluster_group(c) -> str:
    return {7: "cl7", 4: "cl4"}.get(int(c), "minor")


def q(x, p):
    x = pd.Series(x).dropna()
    return float(x.quantile(p)) if len(x) else np.nan


def build_table():
    g = pd.read_csv(DATA / "grid_per_chain.csv")
    ch = pd.read_csv(DATA / "chains.csv")
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv")
    apo = pd.read_csv(DATA / "apo_check.csv")
    ac = pd.read_csv(ROOT / "results/sv/SV12/action_class.csv")
    cc = pd.read_csv(ROOT / "results/sv/SV12/chain_class.csv")
    d = ct[["chain_key", "entry_id", "s08_cluster", "construct_mutations", "ligand_chemotype", "ligand_id_primary",
            "ligand_status"]].merge(ch, on="chain_key").merge(g, on="chain_key")
    d["cluster"] = d.s08_cluster.map(cluster_group)
    d["wild_type"] = d.construct_mutations == "none"
    d["has_ligand"] = d.ligand_ids_used.notna()
    d["ligand_incomplete"] = d.ligand_incomplete.astype("boolean")
    acm = {(r.entry_id, r.ccd): (r.action_class, bool(r.primary_set)) for r in ac.itertuples()}

    def cls(r):
        if not r.has_ligand:
            return ("no ligand of interest", False)
        hits = {acm.get((r.entry_id, c), ("unlabelled", True)) for c in r.ligand_ids_used.split(";")}
        return hits.pop() if len(hits) == 1 else ("unlabelled", False)

    d[["action_class", "class_primary_set"]] = pd.DataFrame(d.apply(cls, axis=1).tolist(), index=d.index)
    d = d.merge(cc[["chain_key", "action_class"]].rename(columns={"action_class": "sv12_chain_class"}),
                on="chain_key", how="left")
    d = d.merge(apo[["chain_key", "passes_apo_check", "residues_ge6_within_4A"]], on="chain_key", how="left")
    d["pocket_volume_A3"] = d.pocket_volume_A3_b10
    d["pocket_volume_b8_A3"] = d.pocket_volume_A3_b8
    d["pocket_volume_b12_A3"] = d.pocket_volume_A3_b12
    d["border_open_points"] = d.border_open_points_b10
    d["pocket_open"] = d.border_open_points > 0
    d["pocket_touches_box"] = d.pocket_touches_box_b10.astype(bool)
    d["ligand_volume_A3"] = np.where(d.has_ligand, d.ligand_volume_A3, np.nan)
    d["ligand_inside_A3"] = np.where(d.has_ligand, d.ligand_points_inside_b10 * VOXEL, np.nan)
    d["inside_share"] = d.ligand_inside_A3 / d.ligand_volume_A3
    d["outside_share"] = 1 - d.inside_share
    d["filling_pct"] = 100 * d.ligand_inside_A3 / d.pocket_volume_A3
    reason = np.where(~d.has_ligand, "no ligand of interest",
                      np.where(d.ligand_incomplete.fillna(False), "ligand_incomplete",
                               np.where(d.pocket_volume_A3.isna(), "no pocket component", "")))
    d["exclusion_reason"] = reason
    d["in_primary_analysis"] = d.exclusion_reason == ""
    return d


def references(d):
    e2 = d[d.wild_type & (d.ligand_ids_used == "EST") & (d.ligand_incomplete == False)]  # noqa: E712
    lig_ref = float(e2.ligand_volume_A3.median())
    e2_pocket = {"median_A3": float(e2.pocket_volume_A3.median()), "q1_A3": q(e2.pocket_volume_A3, .25),
                 "q3_A3": q(e2.pocket_volume_A3, .75), "n_chains": int(len(e2)), "n_entries": int(e2.entry_id.nunique()),
                 "chains": list(e2.chain_key), "definition": "median pocket volume (>= 10 rays) over wild-type chains "
                 "whose analysed ligand is EST alone, complete"}
    cand = d[d.ligand_chemotype == "none"]
    ent = cand.groupby("entry_id").agg(n_chains=("chain_key", "size"), all_pass=("passes_apo_check", "all"),
                                      median_pocket_A3=("pocket_volume_A3", "median"),
                                      blocking=("residues_ge6_within_4A",
                                                lambda s: ";".join(sorted(set(x for x in s.dropna() if x)))))
    passed = ent[ent.all_pass]
    apo = {"candidate_entries": int(len(ent)), "candidate_chains": int(len(cand)),
           "entries_passing": int(len(passed)), "chains_passing": int(cand.passes_apo_check.sum()),
           "entries": {e: {"n_chains": int(r.n_chains), "passes": bool(r.all_pass),
                           "median_pocket_A3": float(r.median_pocket_A3), "blocking_residues": r.blocking}
                       for e, r in ent.iterrows()},
           "median_of_entry_medians_A3": float(passed.median_pocket_A3.median()) if len(passed) else np.nan,
           "q1_A3": q(passed.median_pocket_A3, .25), "q3_A3": q(passed.median_pocket_A3, .75),
           "definition": "ligand_chemotype none and no residue other than the chain's own polymer or water with "
                         ">= 6 heavy atoms within 4 A of the ligand-pocket residues (every kept chain of the entry)"}
    fallback = apo["entries_passing"] < MIN_APO_ENTRIES
    if fallback:
        pref = {"source": "wild-type EST chains (fallback: fewer than 5 apo entries pass)",
                "median_A3": e2_pocket["median_A3"], "q1_A3": e2_pocket["q1_A3"], "q3_A3": e2_pocket["q3_A3"]}
    else:
        pref = {"source": "apo baseline (median of entry medians over apo entries)",
                "median_A3": apo["median_of_entry_medians_A3"], "q1_A3": apo["q1_A3"], "q3_A3": apo["q3_A3"]}
    e2c = json.loads((DATA / "e2_centroid.json").read_text())
    return {"ligand_100pct": {"median_A3": lig_ref, "q1_A3": q(e2.ligand_volume_A3, .25),
                              "q3_A3": q(e2.ligand_volume_A3, .75), "n_chains": int(len(e2)),
                              "definition": "median EST ligand vdW volume (heavy atoms, Bondi) over wild-type EST "
                                            "chains"},
            "pocket_100pct": {**pref, "fallback_used": fallback},
            "apo_baseline": apo, "e2_bound_pocket": e2_pocket, "e2_reference_centroid": e2c,
            "grid": {"spacing_A": 0.4, "voxel_A3": VOXEL, "probe_A": 1.4, "box_margin_A": 6.0, "ray_length_A": 10.0,
                     "rays": 14, "buried_min_rays_primary": 10, "sensitivity": [8, 12]}}


def spearman_row(df, x, y, **meta):
    df = df[[x, y]].dropna()
    r = {**meta, "x": x, "y": y, "n_units": int(len(df))}
    if len(df) < 4 or df[x].nunique() < 2 or df[y].nunique() < 2:
        return {**r, "rho": np.nan, "p": np.nan, "tested": False, "note": "fewer than 4 units or no variation"}
    rho, p = spearmanr(df[x], df[y])
    return {**r, "rho": float(rho), "p": float(p), "tested": True, "note": ""}


def correlations(d, gated):
    a = d[d.in_primary_analysis]
    rows = []
    subsets = {"all": a, "cl7": a[a.cluster == "cl7"], "cl4": a[a.cluster == "cl4"]}
    for c in gated:
        subsets[f"class {c}"] = a[(a.action_class == c) & a.class_primary_set]
    for name, s in subsets.items():
        for level, unit in (("ligand", "ligand_ids_used"), ("entry", "entry_id")):
            m = s.groupby(unit)[["ligand_volume_A3", "pocket_volume_A3"]].median()
            rows.append(spearman_row(m, "ligand_volume_A3", "pocket_volume_A3", subset=name, level=level,
                                     kind="ligand vs pocket volume", in_bh=True,
                                     n_chains=int(len(s)), n_entries=int(s.entry_id.nunique())))
    for t in ("b8", "b12"):
        col = f"pocket_volume_{t}_A3"
        rows.append(spearman_row(a, col, "pocket_volume_A3", subset="all", level="chain",
                                 kind=f"buriedness sensitivity >= {t[1:]} vs >= 10 rays", in_bh=False,
                                 n_chains=int(len(a)), n_entries=int(a.entry_id.nunique())))
        m = a.groupby("entry_id")[[col, "pocket_volume_A3"]].median()
        rows.append(spearman_row(m, col, "pocket_volume_A3", subset="all", level="entry",
                                 kind=f"buriedness sensitivity >= {t[1:]} vs >= 10 rays", in_bh=False,
                                 n_chains=int(len(a)), n_entries=int(a.entry_id.nunique())))
    return pd.DataFrame(rows)


def boot_diff(s, grp_col, g1, g2, rng):
    """Chain-level difference of medians g1 - g2, entries resampled with replacement within each group."""
    parts = []
    for g in (g1, g2):
        x = s[s[grp_col] == g]
        parts.append([v.to_numpy() for _, v in x.groupby("entry_id").filling_pct])
    obs = np.median(np.concatenate(parts[0])) - np.median(np.concatenate(parts[1]))
    bs = np.empty(B)
    for b in range(B):
        m = []
        for p in parts:
            pick = rng.integers(0, len(p), len(p))
            m.append(np.median(np.concatenate([p[i] for i in pick])))
        bs[b] = m[0] - m[1]
    p = min(1.0, (1 + 2 * min((bs <= 0).sum(), (bs >= 0).sum())) / (B + 1))
    return float(obs), float(np.quantile(bs, .025)), float(np.quantile(bs, .975)), float(p)


def entry_values(s, grp_col, groups, cluster_rule):
    """Entry medians of filling % per group. cluster_rule: entry group = group of >= 80 % of its analysed chains."""
    rows = []
    for e, x in s.groupby("entry_id"):
        if cluster_rule:
            share = x[grp_col].value_counts(normalize=True)
            top = share.index[0]
            if share.iloc[0] >= ENTRY_SHARE_MIN and top in groups:
                rows.append({"entry_id": e, "group": top, "value": x[x[grp_col] == top].filling_pct.median()})
        else:
            for g, y in x.groupby(grp_col):
                if g in groups:
                    rows.append({"entry_id": e, "group": g, "value": y.filling_pct.median()})
    return pd.DataFrame(rows)


def entry_table(d, gated):
    """Entry-level filling % for the figure: clusters (>= 80 % rule) and every class (entry x class medians)."""
    a = d[d.in_primary_analysis]
    rows = []
    ev = entry_values(a, "cluster", ["cl7", "cl4", "minor"], True)
    for r in ev.itertuples():
        rows.append({"grouping": "cluster", "group": r.group, "entry_id": r.entry_id, "filling_pct": r.value,
                     "tested": r.group in ("cl7", "cl4")})
    s = a[a.class_primary_set]
    ev = entry_values(s, "action_class", sorted(set(s.action_class)), False)
    for r in ev.itertuples():
        rows.append({"grouping": "action_class", "group": r.group, "entry_id": r.entry_id, "filling_pct": r.value,
                     "tested": (r.group in gated) and len(gated) >= 2})
    return pd.DataFrame(rows)


def filling_tests(d, gated, rng):
    a = d[d.in_primary_analysis]
    rows = []
    comps = [("cluster", "cluster", ["cl7", "cl4"], a, True)]
    if len(gated) >= 2:
        comps.append(("action class", "action_class", list(gated), a[a.class_primary_set], False))
    for name, col, groups, s, crule in comps:
        ev = entry_values(s, col, groups, crule)   # cluster rule: share over all analysed chains of the entry
        s = s[s[col].isin(groups)]
        base = {"comparison": name, "groups": " vs ".join(groups), "outcome": "filling %"}
        for g in groups:
            base[f"chains_{g}"] = int((s[col] == g).sum())
            base[f"entries_{g}"] = int((ev.group == g).sum())
            base[f"chain_median_{g}"] = float(s[s[col] == g].filling_pct.median())
            base[f"entry_median_{g}"] = float(ev[ev.group == g].value.median())
        vals = [ev[ev.group == g].value.to_numpy() for g in groups]
        if len(groups) == 2:
            u = mannwhitneyu(vals[0], vals[1], alternative="two-sided")
            rb = 2 * float(u.statistic) / (len(vals[0]) * len(vals[1])) - 1
            rows.append({**base, "level": "entry", "test": "Mann-Whitney U (two-sided) on entry medians",
                         "statistic": float(u.statistic), "effect": rb,
                         "effect_type": f"rank-biserial r (> 0: {groups[0]} higher)", "ci_low": np.nan,
                         "ci_high": np.nan, "p": float(u.pvalue)})
            obs, lo, hi, p = boot_diff(s, col, groups[0], groups[1], rng)
            rows.append({**base, "level": "chain", "test": f"bootstrap over entries ({B} resamples, seed {SEED})",
                         "statistic": obs, "effect": obs,
                         "effect_type": f"difference of chain medians {groups[0]} - {groups[1]} (percentage points)",
                         "ci_low": lo, "ci_high": hi, "p": p})
        else:
            h = kruskal(*vals)
            rows.append({**base, "level": "entry", "test": "Kruskal-Wallis on entry medians",
                         "statistic": float(h.statistic), "p": float(h.pvalue)})
    return pd.DataFrame(rows)


def by_group(d):
    a = d[d.in_primary_analysis]
    rows = []
    for col in ("cluster", "action_class"):
        s = a if col == "cluster" else a[a.class_primary_set]
        for g, x in s.groupby(col):
            em = x.groupby("entry_id")[["filling_pct", "pocket_volume_A3", "ligand_volume_A3", "inside_share"]].median()
            r = {"grouping": col, "group": g, "n_chains": int(len(x)), "n_entries": int(x.entry_id.nunique()),
                 "n_ligands": int(x.ligand_ids_used.nunique())}
            for v in ("filling_pct", "pocket_volume_A3", "ligand_volume_A3", "inside_share",
                      "ligand_pct_of_E2", "pocket_pct_of_ref"):
                r[f"{v}_chain_median"] = float(x[v].median())
                r[f"{v}_chain_q1"], r[f"{v}_chain_q3"] = q(x[v], .25), q(x[v], .75)
            for v in ("filling_pct", "pocket_volume_A3", "ligand_volume_A3", "inside_share"):
                r[f"{v}_entry_median"] = float(em[v].median())
            r["pocket_touches_box_chains"] = int(x.pocket_touches_box.sum())
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    d = build_table()
    ref = references(d)
    L, P = ref["ligand_100pct"]["median_A3"], ref["pocket_100pct"]["median_A3"]
    d["ligand_pct_of_E2"] = 100 * d.ligand_volume_A3 / L
    d["pocket_pct_of_ref"] = 100 * d.pocket_volume_A3 / P
    d["pocket_pct_of_apo_baseline"] = 100 * d.pocket_volume_A3 / ref["apo_baseline"]["median_of_entry_medians_A3"]
    d["pocket_pct_of_E2_pocket"] = 100 * d.pocket_volume_A3 / ref["e2_bound_pocket"]["median_A3"]
    cov = pd.read_csv(ROOT / "results/sv/SV12/class_coverage.csv")
    gate12 = cov[(cov.set == "primary") & (cov.source == "any") & (cov.cluster == "all")
                 & (cov.entries_ge5.astype(str) == "True")].action_class.tolist()
    a = d[d.in_primary_analysis & d.class_primary_set]
    n_ent = {c: int(a[a.action_class == c].entry_id.nunique()) for c in gate12}
    gated = [c for c in gate12 if n_ent[c] >= GATE_ENTRIES]
    ref["class_gate"] = {"sv12_gate_classes": gate12, "entries_in_sv14_primary_set": n_ent, "analysed": gated}
    cols = ["chain_key", "entry_id", "s08_cluster", "cluster", "construct_mutations", "wild_type", "ligand_chemotype",
            "ligand_id_primary", "ligand_ids_used", "n_instances_q18", "n_instances_used", "instance_rule",
            "lig_heavy_modelled", "lig_heavy_ccd", "ligand_incomplete", "action_class", "class_primary_set",
            "sv12_chain_class", "pocket_volume_A3", "pocket_volume_b8_A3", "pocket_volume_b12_A3", "pocket_open",
            "border_open_points", "border_open_points_b8", "border_open_points_b12", "pocket_touches_box",
            "pocket_touches_box_b8", "pocket_touches_box_b12", "n_overlapping_components", "pocket_rule_b10",
            "ligand_volume_A3", "ligand_inside_A3", "inside_share", "outside_share", "filling_pct",
            "ligand_pct_of_E2", "pocket_pct_of_ref", "pocket_pct_of_apo_baseline", "pocket_pct_of_E2_pocket",
            "ligand_points_outside_box", "ligand_points_in_receptor_vdw", "lig_centroid_to_e2_A",
            "pocket_centroid_to_e2_A", "e2_centroid_to_nearest_point_A", "passes_apo_check",
            "residues_ge6_within_4A", "in_primary_analysis", "exclusion_reason"]
    d[cols].to_csv(OUT / "pocket_per_chain.csv", index=False)
    cor = correlations(d, gated)
    tests = filling_tests(d, gated, rng)
    pv = list(cor.loc[cor.in_bh & cor.tested, "p"]) + list(tests.p)
    qv = multipletests(pv, method="fdr_bh")[1] if pv else []
    k = int((cor.in_bh & cor.tested).sum())
    cor["q_bh"] = np.nan
    cor.loc[cor.in_bh & cor.tested, "q_bh"] = qv[:k]
    tests["q_bh"] = qv[k:]
    cor.to_csv(OUT / "correlations.csv", index=False)
    tests.to_csv(OUT / "tests.csv", index=False)
    bg = by_group(d)
    bg.to_csv(OUT / "filling_by_group.csv", index=False)
    entry_table(d, gated).to_csv(OUT / "filling_entry_level.csv", index=False)
    (OUT / "references.json").write_text(json.dumps(ref, indent=1, default=float))
    prim = d[d.in_primary_analysis]
    numbers = {
        "chains": int(len(d)), "chains_with_ligand": int(d.has_ligand.sum()),
        "ligand_incomplete_chains": int(d.ligand_incomplete.fillna(False).sum()),
        "ligand_incomplete_entries": int(d[d.ligand_incomplete.fillna(False)].entry_id.nunique()),
        "primary_chains": int(len(prim)), "primary_entries": int(prim.entry_id.nunique()),
        "primary_ligands": int(prim.ligand_ids_used.nunique()),
        "primary_by_cluster": prim.cluster.value_counts().to_dict(),
        "class_counts_primary_chains": prim[prim.class_primary_set].action_class.value_counts().to_dict(),
        "class_mismatch_vs_sv12_chain_class": int((d.has_ligand & (d.action_class != d.sv12_chain_class)).sum()),
        "class_mismatch_chains": list(d[d.has_ligand & (d.action_class != d.sv12_chain_class)].chain_key),
        "multi_instance_chains": int((d.n_instances_q18 > 1).sum()),
        "alt_instance_rule_chains": int(d.instance_rule.fillna("").str.startswith("alt id").sum()),
        "pocket_open_chains": int(d.pocket_open.sum()),
        "border_points_median": float(d.border_open_points.median()),
        "border_points_q1": q(d.border_open_points, .25), "border_points_q3": q(d.border_open_points, .75),
        "touches_box_b10_chains": int(d.pocket_touches_box.sum()),
        "touches_box_b8_chains": int(d.pocket_touches_box_b8.astype(bool).sum()),
        "touches_box_b12_chains": int(d.pocket_touches_box_b12.astype(bool).sum()),
        "touches_box_b10_primary": int(prim.pocket_touches_box.sum()),
        "touches_box_b10_by_cluster_primary": prim[prim.pocket_touches_box].cluster.value_counts().to_dict(),
        "ligand_outside_box_chains": int((d.ligand_points_outside_box > 0).sum()),
        "overlap_gt1_component_chains": int((d.n_overlapping_components > 1).sum()),
        "pocket_volume_primary_median": float(prim.pocket_volume_A3.median()),
        "pocket_volume_primary_q1": q(prim.pocket_volume_A3, .25), "pocket_volume_primary_q3": q(prim.pocket_volume_A3, .75),
        "pocket_volume_primary_min": float(prim.pocket_volume_A3.min()),
        "pocket_volume_primary_max": float(prim.pocket_volume_A3.max()),
        "ligand_volume_primary_median": float(prim.ligand_volume_A3.median()),
        "filling_primary_median": float(prim.filling_pct.median()),
        "filling_primary_q1": q(prim.filling_pct, .25), "filling_primary_q3": q(prim.filling_pct, .75),
        "inside_share_primary_median": float(prim.inside_share.median()),
        "inside_share_lt_0_9_chains": int((prim.inside_share < 0.9).sum()),
        "e2_ref_pocket_values": sorted(float(v) for v in d[d.chain_key.isin(ref["e2_bound_pocket"]["chains"])]
                                       .pocket_volume_A3),
        "e2_ref_filling_values": sorted(float(v) for v in d[d.chain_key.isin(ref["e2_bound_pocket"]["chains"])]
                                        .filling_pct),
        "e2_ref_border_points": sorted(int(v) for v in d[d.chain_key.isin(ref["e2_bound_pocket"]["chains"])]
                                       .border_open_points),
        "apo_chain_pockets": {r.chain_key: float(r.pocket_volume_A3) for r in d[~d.has_ligand].itertuples()},
    }
    (OUT / "numbers.json").write_text(json.dumps(numbers, indent=1, default=float))
    print(json.dumps({"references": {k: ref[k] for k in ("ligand_100pct", "pocket_100pct")}, "gate": ref["class_gate"]},
                     indent=1, default=float))
    print(cor.to_string())
    print(tests.to_string())


if __name__ == "__main__":
    main()
