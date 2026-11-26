"""Stage 12.4 - Effects of curated action class on H12 inside one H12 position.

DOI-block bootstrap of median differences, leave-one-DOI-out and within-DOI checks, robustness filters, and the
curated class of every examined lid outlier. Tables only; the figures are drawn by pipeline/figures/fig_*.py.
Outputs (results/12_curated_action/tables/): effects_summary.csv, lodo_within_doi.csv, robustness_sensitivity.csv,
outliers_curated_audit.csv, pocket_by_class.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402


def compute_effects(pg: pd.DataFrame, chains: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute pre-registered contrasts inside cl7 and cl4, LODO, within-DOI, and robustness filters."""
    core_pg = pg[pg["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()
    effects = []

    def add_contrast(name, group, stratum, sub_df, cls_a, cls_b, val_col="dist_H12c_H11axis_A", exploratory=False):
        da = sub_df[sub_df["action_class"] == cls_a]
        db = sub_df[sub_df["action_class"] == cls_b]
        res = cc.block_bootstrap_diff_medians(da, db, val_col, block_col="doi_block", n_boot=2000, seed=cc.SEED)
        thr = cc.ORDINARY_THR_A.get(group, np.nan)
        shared_dois = len(set(da.dropna(subset=[val_col])["doi_block"]) & set(db.dropna(subset=[val_col])["doi_block"]))
        effects.append({
            "contrast_id": name,
            "group": group,
            "stratum": stratum,
            "metric": val_col,
            "class_a": cls_a,
            "class_b": cls_b,
            "n_a": res["n_a"],
            "n_b": res["n_b"],
            "pdbs_a": int(da.dropna(subset=[val_col])["pdb_norm"].nunique()),
            "pdbs_b": int(db.dropna(subset=[val_col])["pdb_norm"].nunique()),
            "inchikeys_a": int(da.dropna(subset=[val_col])["inchikey"].nunique()),
            "inchikeys_b": int(db.dropna(subset=[val_col])["inchikey"].nunique()),
            "dois_a": res["blocks_a"],
            "dois_b": res["blocks_b"],
            "shared_dois": int(shared_dois),
            "median_a": res["med_a"],
            "median_b": res["med_b"],
            "diff_a_minus_b": res["diff"],
            "ci95_lo": res["ci_lo"],
            "ci95_hi": res["ci_hi"],
            "ci_valid": res["ci_valid"],
            "p_mwu": res["p_mwu"],
            "ordinary_threshold_A": thr,
            "exceeds_ordinary": bool(abs(res["diff"]) > thr) if np.isfinite(res["diff"]) and np.isfinite(thr) else False,
            "exploratory": exploratory,
        })

    # Pre-registered primary strata
    cl7_all = core_pg[core_pg["group"] == "cl7"]
    cl4_all = core_pg[core_pg["group"] == "cl4"]
    cl7_y537s_pep = cl7_all[(cl7_all["construct_mutations"] == "Y537S") & (cl7_all["peptide_bound"].astype(str) == "True")]
    cl4_l536s_nopep = cl4_all[(cl4_all["mutation_status"] == "L536S") & (cl4_all["peptide_bound"].astype(str) == "False")]
    cl4_l536s_exact = cl4_all[(cl4_all["construct_mutations"] == "L372S;L536S") & (cl4_all["peptide_bound"].astype(str) == "False")]

    # Primary contrasts (pre-registered)
    add_contrast("cl7_Y537S_pep: ant vs ago", "cl7", "Y537S + peptide (exact)", cl7_y537s_pep, "antagonist", "agonist")
    add_contrast("cl7_Y537S_pep: part vs ago", "cl7", "Y537S + peptide (exact)", cl7_y537s_pep, "partial_agonist", "agonist")
    add_contrast("cl7_Y537S_pep: ant vs part", "cl7", "Y537S + peptide (exact)", cl7_y537s_pep, "antagonist", "partial_agonist")
    add_contrast("cl4_L536S_nopep: ant vs part", "cl4", "L536S group, no peptide", cl4_l536s_nopep, "antagonist", "partial_agonist")
    add_contrast("cl4_L372S_L536S: ant vs part", "cl4", "L372S;L536S (exact), no pep", cl4_l536s_exact, "antagonist", "partial_agonist")

    # Crude (unstratified) within-group contrasts for comparison
    add_contrast("cl7_crude: ant vs ago", "cl7", "all cl7 (crude)", cl7_all, "antagonist", "agonist", exploratory=True)
    add_contrast("cl7_crude: part vs ago", "cl7", "all cl7 (crude)", cl7_all, "partial_agonist", "agonist", exploratory=True)
    add_contrast("cl4_crude: ant vs part", "cl4", "all cl4 (crude)", cl4_all, "antagonist", "partial_agonist", exploratory=True)

    # Secondary geometry metrics inside cl7 Y537S + peptide (exploratory)
    for col in ("angle_H12_H3_deg", "angle_H12_H11_deg", "len_H12", "mobility_h12_z"):
        add_contrast(f"cl7_Y537S_pep_{col}: ant vs ago", "cl7", "Y537S + peptide (exact)", cl7_y537s_pep, "antagonist", "agonist", val_col=col, exploratory=True)
        add_contrast(f"cl7_Y537S_pep_{col}: part vs ago", "cl7", "Y537S + peptide (exact)", cl7_y537s_pep, "partial_agonist", "agonist", val_col=col, exploratory=True)

    eff_df = pd.DataFrame(effects)
    eff_df["q_bh"] = cc.bh_fdr(eff_df["p_mwu"])
    eff_df.to_csv(cc.TABLES_DIR / "effects_summary.csv", index=False)

    # Leave-one-DOI-out (LODO) & within-DOI for main contrast (cl7 Y537S + peptide: antagonist vs agonist, and partial vs agonist)
    lodo_rows = []
    for cname, ca, cb in [("ant_vs_ago", "antagonist", "agonist"), ("part_vs_ago", "partial_agonist", "agonist")]:
        sub = cl7_y537s_pep[cl7_y537s_pep["action_class"].isin([ca, cb])].dropna(subset=["dist_H12c_H11axis_A"]).copy()
        full_diff = float(sub[sub["action_class"] == ca]["dist_H12c_H11axis_A"].median() - sub[sub["action_class"] == cb]["dist_H12c_H11axis_A"].median())
        top_dois = sub["doi_block"].value_counts().head(12).index.tolist()
        for doi in top_dois:
            s_loo = sub[sub["doi_block"] != doi]
            da, db = s_loo[s_loo["action_class"] == ca]["dist_H12c_H11axis_A"], s_loo[s_loo["action_class"] == cb]["dist_H12c_H11axis_A"]
            diff_loo = float(da.median() - db.median()) if len(da) and len(db) else np.nan
            s_in = sub[sub["doi_block"] == doi]
            in_a, in_b = s_in[s_in["action_class"] == ca]["dist_H12c_H11axis_A"], s_in[s_in["action_class"] == cb]["dist_H12c_H11axis_A"]
            lodo_rows.append({
                "contrast": cname,
                "doi_block": doi,
                "n_in_doi_a": int(len(in_a)),
                "n_in_doi_b": int(len(in_b)),
                "within_doi_diff": float(in_a.median() - in_b.median()) if len(in_a) and len(in_b) else np.nan,
                "full_diff": full_diff,
                "lodo_diff": diff_loo,
                "delta_lodo": diff_loo - full_diff if np.isfinite(diff_loo) else np.nan,
            })
    lodo_df = pd.DataFrame(lodo_rows)
    lodo_df.to_csv(cc.TABLES_DIR / "lodo_within_doi.csv", index=False)

    # Robustness sensitivity across 9 filters for main contrast: cl7 Y537S + peptide (antagonist vs agonist & partial vs agonist)
    # Build "with multi" version from chains directly
    ch_c3m = chains[(chains["in_core3_with_multi"].astype(str) == "True") & (chains["group"] == "cl7") &
                    (chains["construct_mutations"] == "Y537S") & (chains["peptide_bound"].astype(str) == "True")].copy()
    pg_multi = ch_c3m.groupby(["pdb_norm", "lig_norm"], as_index=False).agg(
        action_class=("action_class", "first"),
        doi_block=("doi_block", "first"),
        inchikey=("inchikey", "first"),
        dist_H12c_H11axis_A=("dist_H12c_H11axis_A", "median"),
    )

    top_doi_cl7 = cl7_y537s_pep["doi_block"].value_counts().index[0]
    filters = [
        ("1. Strict single (base)", cl7_y537s_pep),
        ("2. Include multi-ligand", pg_multi),
        ("3. QC Tier A only", cl7_y537s_pep[cl7_y537s_pep["qc_label"] == "A"]),
        ("4. R-free <= 0.30", cl7_y537s_pep[pd.to_numeric(cl7_y537s_pep["r_free"], errors="coerce") <= 0.30]),
        ("5. H12 complete 538-547", cl7_y537s_pep[cl7_y537s_pep["h12_complete_547"].astype(str) == "True"]),
        ("6. H12 complete 538-544", cl7_y537s_pep[cl7_y537s_pep["h12_complete_544"].astype(str) == "True"]),
        ("7. Direct lit (excl derived)", cl7_y537s_pep[cl7_y537s_pep["primary_source"].isin(["full_text", "abstract", "literature"])]),
        ("8. Excl derived_from_data", cl7_y537s_pep[cl7_y537s_pep["primary_source"] != "derived_from_data"]),
        ("9. Excl EST & WT-mismatch", cl7_y537s_pep[(cl7_y537s_pep["lig_norm"] != "EST") & (cl7_y537s_pep["wt_assay_on_mutant_note"].astype(str) == "False")]),
        (f"10. Excl top DOI ({top_doi_cl7[:18]})", cl7_y537s_pep[cl7_y537s_pep["doi_block"] != top_doi_cl7]),
        ("11. Equal InChIKey weight", cl7_y537s_pep.groupby(["action_class", "inchikey"], as_index=False).agg(
            dist_H12c_H11axis_A=("dist_H12c_H11axis_A", "median"), doi_block=("doi_block", "first"), pdb_norm=("pdb_norm", "first")
        )),
    ]
    rob_rows = []
    for fname, fdf in filters:
        for cname, ca, cb in [("ant_vs_ago", "antagonist", "agonist"), ("part_vs_ago", "partial_agonist", "agonist")]:
            da = fdf[fdf["action_class"] == ca]
            db = fdf[fdf["action_class"] == cb]
            res = cc.block_bootstrap_diff_medians(da, db, "dist_H12c_H11axis_A", block_col="doi_block", n_boot=1500, seed=cc.SEED)
            rob_rows.append({
                "filter": fname,
                "contrast": cname,
                "n_a": res["n_a"],
                "n_b": res["n_b"],
                "dois_a": res["blocks_a"],
                "dois_b": res["blocks_b"],
                "diff_A": res["diff"],
                "ci95_lo": res["ci_lo"],
                "ci95_hi": res["ci_hi"],
                "ci_valid": res["ci_valid"],
                "p_mwu": res["p_mwu"],
            })
    rob_df = pd.DataFrame(rob_rows)
    rob_df.to_csv(cc.TABLES_DIR / "robustness_sensitivity.csv", index=False)
    return eff_df, lodo_df, rob_df


def outlier_audit(chains: pd.DataFrame) -> pd.DataFrame:
    """Examined lid outliers (stage 03) with their curated action class and evidence source."""
    out = cc.read_csv_safe(cc.RESULTS_BASE / "03_lid_outliers/lid_outliers.csv")
    out = out.merge(chains[["chain_key", "action_class", "doi", "primary_source", "drug_classes"]], on="chain_key",
                    how="left", suffixes=("_old", "_cur"), validate="one_to_one")
    out.to_csv(cc.TABLES_DIR / "outliers_curated_audit.csv", index=False)
    return out


def pocket_by_class(pg: pd.DataFrame) -> pd.DataFrame:
    """Pocket volume (% of the E2 pocket) and filling by curated class inside cl7 and cl4; pairs whose ligand is fully
    modelled; Kruskal-Wallis across the classes with at least 3 pairs."""
    d = pg[pg["action_class"].isin(cc.CLASS_ORDER_CORE) & pg["pocket_complete_primary"].astype(str).eq("True")]
    rows = []
    for g in ("cl7", "cl4"):
        for col in ("pocket_pct_of_E2_pocket", "filling_pct"):
            x = {k: d[(d["group"] == g) & (d["action_class"] == k)][col].dropna() for k in cc.CLASS_ORDER_CORE}
            ok = [v for v in x.values() if len(v) >= 3]
            row = {"group": g, "metric": col, "kruskal_p": float(kruskal(*ok).pvalue) if len(ok) >= 2 else np.nan}
            for k, v in x.items():
                row[f"n_{k}"], row[f"median_{k}"] = int(len(v)), float(v.median()) if len(v) else np.nan
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(cc.TABLES_DIR / "pocket_by_class.csv", index=False)
    return out


def main():
    cc.ensure_dirs()
    chains = cc.read_csv_safe(cc.TABLES_DIR / "chains_curated.csv")
    pg = cc.read_csv_safe(cc.TABLES_DIR / "pairs_by_group.csv")
    eff_df, _, rob_df = compute_effects(pg, chains)
    outlier_audit(chains)
    pocket_by_class(pg)
    print(eff_df[["contrast_id", "n_a", "n_b", "diff_a_minus_b", "ci95_lo", "ci95_hi", "q_bh"]].round(3).to_string())


if __name__ == "__main__":
    main()
