"""Stage 12.4a - Statistical effect estimation (DOI-block bootstrap, LODO, within-DOI, robustness) and Figures F01-F12."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import mannwhitneyu, spearmanr

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


def make_priority_figures(chains, pairs, pg, eff_df, lodo_df, manifest):
    cc.setup_mpl()
    core_pairs = pairs[pairs["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()
    core_pg = pg[pg["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()

    # ========================================================================== F01: Redundancy & Accumulation Curves
    curves = cc.read_csv_safe(cc.TABLES_DIR / "redundancy_accumulation_curves.csv")
    hill = cc.read_csv_safe(cc.TABLES_DIR / "redundancy_hill_summary.csv")
    sc_steps = cc.read_csv_safe(cc.TABLES_DIR / "greedy_set_cover_steps.csv")
    sc_sum = cc.read_csv_safe(cc.TABLES_DIR / "greedy_set_cover_summary.csv")

    fig, axs = plt.subplots(2, 2, figsize=(12.2, 8.6))
    # Panel A: All 351 PDBs accumulation across 4 RMSD cuts
    ax = axs[0, 0]
    cut_colors = {1.0: "#4a3aa7", 1.5: "#2a78d6", 2.0: "#1b9e77", 2.5: "#8a8984"}
    for cut in (1.0, 1.5, 2.0, 2.5):
        d = curves[(curves["subset"] == "all_351_pdbs") & (curves["cut_A"] == cut)]
        ax.plot(d["n_pdbs"], d["mean_clusters"], color=cut_colors[cut], lw=2.0,
                label=f"{cut:.1f} Å (total={int(d['total_clusters_denom'].iloc[0])}, n90={int(d['n90_pdbs'].iloc[0])})")
        ax.fill_between(d["n_pdbs"], d["p10_clusters"], d["p90_clusters"], color=cut_colors[cut], alpha=0.12)
    ax.set_xlabel("Randomly added PDB entries (without replacement)")
    ax.set_ylabel("Unique clusters discovered")
    ax.set_title("A. All 351 PDBs: cluster discovery across RMSD cuts")
    ax.legend(loc="lower right", fontsize=8)
    cc.grid(ax)

    # Panel B: Class comparison at 1.5 A on common PDB range
    ax = axs[0, 1]
    min_pdbs = min(int(curves[(curves["subset"] == cls) & (curves["cut_A"] == 1.5)]["total_pdbs"].iloc[0])
                   for cls in cc.CLASS_ORDER_CORE)
    for cls in cc.CLASS_ORDER_CORE:
        d = curves[(curves["subset"] == cls) & (curves["cut_A"] == 1.5)]
        dc = d[d["n_pdbs"] <= min_pdbs]
        tot_c = int(d["total_clusters_denom"].iloc[0])
        n90 = int(d["n90_pdbs"].iloc[0])
        ax.plot(dc["n_pdbs"], dc["mean_clusters"], color=cc.CLASS_COLOR[cls], lw=2.2,
                label=f"{cls} (n={int(d['total_pdbs'].iloc[0])} PDB, q0={tot_c}, n90={n90})")
        ax.fill_between(dc["n_pdbs"], dc["p10_clusters"], dc["p90_clusters"], color=cc.CLASS_COLOR[cls], alpha=0.15)
    ax.axvline(min_pdbs, color=cc.MUTED, ls="--", lw=1)
    ax.set_xlabel(f"Randomly added PDB entries (common range 1–{min_pdbs})")
    ax.set_ylabel("Unique 1.5 Å clusters discovered")
    ax.set_title("B. Cluster discovery by functional class (1.5 Å cut)")
    ax.legend(loc="lower right", fontsize=8)
    cc.grid(ax)

    # Panel C: Effective clusters (Hill q1) & Shannon bits (PDB-weighted, 1.5 A) + sensitivity
    ax = axs[1, 0]
    sub_order = ["all_351_pdbs", "agonist", "partial_agonist", "antagonist", "core3_dedup_inchikey", "core3_dedup_doi", "core3_tierA"]
    sub_labels = ["All 351 PDB", "agonist", "partial_ag.", "antagonist", "1 PDB/InChIKey", "1 PDB/DOI", "QC Tier A"]
    h15 = hill[(hill["cut_A"] == 1.5) & (hill["subset"].isin(sub_order))].set_index("subset").loc[sub_order]
    y = np.arange(len(sub_order))
    bar_cols = [cc.INK2, cc.CLASS_COLOR["agonist"], cc.CLASS_COLOR["partial_agonist"], cc.CLASS_COLOR["antagonist"],
                "#6da7ec", "#6da7ec", "#6da7ec"]
    ax.barh(y, h15["hill_q1"], color=bar_cols, alpha=0.88, height=0.62)
    for i, r in enumerate(h15.itertuples()):
        ax.text(r.hill_q1 + 0.05, i, f"q1={r.hill_q1:.2f} | H={r.shannon_bits:.2f}b | q0={int(r.q0_clusters)} (n={int(r.n_pdbs)})",
                va="center", fontsize=7.8, color=cc.INK)
    ax.set_yticks(y, sub_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 3.8)
    ax.set_xlabel("Effective number of conformations (Hill q1, equal PDB weight)")
    ax.set_title("C. Effective diversity & sensitivity (1.5 Å, equal PDB weight)")
    cc.grid(ax, "x")

    # Panel D: Greedy set cover (fraction of 788 chains covered within radius r)
    ax = axs[1, 1]
    r_colors = {0.5: "#c2185b", 0.75: "#7b3294", 1.0: "#2a78d6", 1.5: "#1b9e77"}
    for r in (0.5, 0.75, 1.0, 1.5):
        g = sc_steps[sc_steps["radius_A"] == r]
        r90 = sc_sum.loc[sc_sum["radius_A"] == r, "reps_90pct"].iloc[0]
        lbl = f"r={r:.2f} Å (90% in {int(r90)} reps)" if pd.notna(r90) else f"r={r:.2f} Å"
        ax.plot(g["step"], g["cum_coverage"] * 100, color=r_colors[r], lw=2.0, marker="o", ms=3,
                markevery=max(1, len(g) // 12), label=lbl)
    ax.axhline(90, color=cc.MUTED, ls="--", lw=1)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 103)
    ax.set_xlabel("Complete-H12 representative chains selected (greedy cover)")
    ax.set_ylabel("Chains covered within RMSD radius (%)")
    ax.set_title("D. Chain coverage by representatives (vs cluster discovery in A–B)")
    ax.legend(loc="lower right", fontsize=8)
    cc.grid(ax)

    cc.wrap_title(fig, "Conformational redundancy is extreme: 2 representatives cover 90% of chains within 1.0 Å, and effective PDB diversity stays near Hill q1 ≈ 1.2–2.1 across classes.", width=105)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    cc.save_fig_bundle(
        fig, "F01_redundancy_accumulation", h15.reset_index(),
        "A–B: Unique clusters discovered by random PDB addition without replacement (200 permutations, 10–90% band; n90 uses each subset's observed cluster total as denominator). C: Equal-PDB-weight Hill q1 and Shannon entropy (bits) at 1.5 Å across classes and deduplication/QC sensitivities. D: Greedy set cover of 788 chains by complete-H12 representatives (2 reps at 1.0 Å, 48 at 0.5 Å for 90% chain coverage — note this measures chain proximity, not cluster count). Sources: results/02_lid_groups/groups.csv, rmsd.npy, esr1_dataset_final.csv.",
        manifest, "Redundancy and cluster accumulation overall and by curated class",
        "results/02_lid_groups/groups.csv; results/02_lid_groups/rmsd.npy; esr1_dataset_final.csv",
        "defense", "PDB entries (A–C), chains (D)"
    )

    # ========================================================================== F02: Chemical space (ECFP4 MDS + stability)
    chem = cc.read_csv_safe(cc.TABLES_DIR / "chem_ligands_order.csv")
    stab = cc.read_csv_safe(cc.TABLES_DIR / "chem_cluster_stability.csv")
    disc = cc.read_csv_safe(cc.TABLES_DIR / "chem_discordant_pairs.csv")

    fig, axs = plt.subplots(1, 2, figsize=(12.2, 5.4), gridspec_kw={"width_ratios": [1.45, 1.0]})
    ax = axs[0]
    for cls in cc.CLASS_ORDER_CORE:
        for grp in ("cl7", "cl4", "minor", "multiple"):
            d = chem[(chem["action_class_majority"] == cls) & (chem["h12_group_majority"] == grp)]
            if len(d) == 0:
                continue
            ax.scatter(
                d["mds_x"], d["mds_y"],
                c=cc.CLASS_COLOR[cls], marker=cc.GROUP_MARKER[grp],
                s=46, alpha=0.85, edgecolors="white", linewidths=0.5,
                label=f"{cls} | {grp} (n={len(d)})",
            )
    # Highlight conflict ligands (same CCD with multiple classes or groups)
    conf = chem[(chem["action_conflict"].astype(str) == "True") | (chem["h12_conflict"].astype(str) == "True")]
    if len(conf):
        ax.scatter(conf["mds_x"], conf["mds_y"], s=110, facecolors="none", edgecolors=cc.INK, linewidths=1.3, zorder=5,
                   label=f"multi-context CCD (n={len(conf)})")
        for _, r in conf.head(5).iterrows():
            ax.annotate(r["lig_norm"], (r["mds_x"], r["mds_y"]), xytext=(5, 5), textcoords="offset points",
                        fontsize=7.5, fontweight="bold", color=cc.INK)
    # Connect top 2 discordant close pairs
    for _, r in disc.head(2).iterrows():
        a = chem[chem["lig_norm"] == r["lig_a"]].iloc[0]
        b = chem[chem["lig_norm"] == r["lig_b"]].iloc[0]
        ax.plot([a["mds_x"], b["mds_x"]], [a["mds_y"], b["mds_y"]], color=cc.INK, lw=1.2, ls="--", zorder=4)
        mx, my = (a["mds_x"] + b["mds_x"]) / 2, (a["mds_y"] + b["mds_y"]) / 2
        ax.annotate(f"{r['lig_a']}–{r['lig_b']}\nTc={r['tanimoto_ecfp4']:.2f}", (mx, my),
                    xytext=(6, -12), textcoords="offset points", fontsize=7.5,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=cc.MUTED, lw=0.7))
    ax.set_xlabel("ECFP4 Tanimoto MDS 1")
    ax.set_ylabel("ECFP4 Tanimoto MDS 2")
    ax.set_title(f"A. Chemical space of {len(chem)} unique CCD ligands (color=class, shape=H12)")
    ax.legend(loc="best", fontsize=7.2, ncol=2)
    cc.grid(ax, "both")

    ax = axs[1]
    piv = stab.pivot(index="clustering_a", columns="clustering_b", values="ari")
    labels_short = [c.replace("cl_", "").replace("_0p", " d=") for c in piv.index]
    im = ax.imshow(piv.values, cmap="Blues", vmin=0, vmax=1)
    for i in range(len(piv)):
        for j in range(len(piv)):
            val = piv.values[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8,
                    color="white" if val > 0.55 else cc.INK)
    ax.set_xticks(range(len(labels_short)), labels_short, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(labels_short)), labels_short, fontsize=8)
    ax.set_title("B. Chemical cluster stability (ARI across cutoffs & ECFP4/MACCS)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Adjusted Rand Index")

    cc.wrap_title(fig, f"Chemical space separates agonist/partial-agonist scaffolds from antagonist scaffolds, with {len(disc)} close ligand pairs (Tc ≥ 0.55) switching class or H12 state.", width=102)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    cc.save_fig_bundle(
        fig, "F02_chem_space", chem[["lig_norm", "inchikey", "action_class_majority", "h12_group_majority", "cl_ecfp4_0p6", "mds_x", "mds_y", "n_pairs"]],
        f"A: 2D MDS of 1−Tanimoto (RDKit ECFP4, 2048 bits, radius 2) across {len(chem)} unique CCDs from 291 strict core-3 pairs (SMILES 291/291 matched). Color encodes curated functional class; marker encodes majority H12 group; black rings mark multi-context CCDs; dashed segments connect top discordant close pairs. B: Pairwise ARI across distance cutoffs (0.4–0.6) and fingerprints (ECFP4 vs MACCS). Sources: data/ligands/ligand_ids.csv, esr1_dataset_final.csv.",
        manifest, "Chemical similarity, functional class, and H12 conformation",
        "data/ligands/ligand_ids.csv; esr1_dataset_final.csv",
        "defense", "unique CCD ligands (n=240+)"
    )

    # ========================================================================== F03: Chem clusters on Sammon map + contingency
    maps = cc.read_csv_safe(cc.RESULTS_BASE / "05_maps/maps_coords.csv")
    pairs_chem = cc.read_csv_safe(cc.TABLES_DIR / "pairs_core3_with_chem.csv")
    assoc = cc.read_csv_safe(cc.TABLES_DIR / "chem_struct_class_associations.csv")
    # Map primary chemical cluster onto chains for Sammon plot
    ch_map = chains.merge(maps[["chain_key", "sammon1", "sammon2"]], on="chain_key", how="left")
    ch_map = ch_map.merge(chem[["lig_norm", "cl_ecfp4_0p6"]], on="lig_norm", how="left")
    top5_cl = chem["cl_ecfp4_0p6"].value_counts().head(5).index.tolist()

    fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.5))
    ax = axs[0]
    bg = ch_map[ch_map["cl_ecfp4_0p6"].isna() | ~ch_map["cl_ecfp4_0p6"].isin(top5_cl)]
    ax.scatter(bg["sammon1"], bg["sammon2"], c="#dcdbd6", s=16, alpha=0.6, label=f"other/unmatched (n={len(bg)})")
    pal5 = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    for idx_c, cid in enumerate(top5_cl):
        d = ch_map[ch_map["cl_ecfp4_0p6"] == cid]
        ax.scatter(d["sammon1"], d["sammon2"], c=pal5[idx_c], s=26, alpha=0.85, label=f"Chem C{int(cid)} (n={len(d)} ch)")
    ax.set_xlabel("Sammon 1 (structural RMSD map)")
    ax.set_ylabel("Sammon 2")
    ax.set_title("A. Top-5 ECFP4 clusters on structural map")
    ax.legend(loc="best", fontsize=7.5)
    cc.grid(ax, "both")

    # Panel B & C: Contingency heatmaps (top 8 chem clusters + other)
    top8 = pairs_chem["cl_ecfp4_0p6"].value_counts().head(8).index.tolist()
    pc = pairs_chem.copy()
    pc["chem_grp"] = pc["cl_ecfp4_0p6"].map(lambda x: f"C{int(x)}" if x in top8 else "Other")
    row_ord = [f"C{int(x)}" for x in top8] + ["Other"]

    for ax, col_var, col_ord, title_txt, target_key in [
        (axs[1], "h12_group_excl", ["cl7", "cl4", "minor", "multiple"], "B. Chem cluster × H12 group (pairs)", "h12_group"),
        (axs[2], "action_class", cc.CLASS_ORDER_CORE, "C. Chem cluster × functional class (pairs)", "action_class"),
    ]:
        tab = pd.crosstab(pc["chem_grp"], pc[col_var]).reindex(index=row_ord, columns=col_ord, fill_value=0)
        frac = tab.div(tab.sum(axis=1).replace(0, 1), axis=0)
        ax.imshow(frac.values, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
        for i in range(tab.shape[0]):
            for j in range(tab.shape[1]):
                cnt = int(tab.values[i, j])
                ax.text(j, i, f"{cnt}" if cnt > 0 else "·", ha="center", va="center", fontsize=8,
                        color="white" if frac.values[i, j] > 0.55 else cc.INK)
        ax.set_xticks(range(len(col_ord)), col_ord, fontsize=8.5)
        ax.set_yticks(range(len(row_ord)), [f"{r} (n={int(tab.loc[r].sum())})" for r in row_ord], fontsize=8)
        v_pair = assoc[(assoc["level"] == "pairs_core3") & (assoc["target"] == target_key)]["cramers_v"].iloc[0]
        v_lig = assoc[(assoc["level"] == "unique_ligands") & (assoc["target"] == target_key)]["cramers_v"].iloc[0]
        v_doi = assoc[(assoc["level"] == "doi_dedup_pairs") & (assoc["target"] == target_key)]["cramers_v"].iloc[0]
        ax.set_title(f"{title_txt}\nV={v_pair:.2f} (lig {v_lig:.2f}, DOI {v_doi:.2f})", fontsize=9)

    cc.wrap_title(fig, "Chemical scaffold families strongly track both H12 conformation and functional class, even after ligand and DOI deduplication.", width=102)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F03_chem_on_struct_map", assoc,
        "A: Top-5 ECFP4 chemical clusters (d=0.6) projected onto the existing Stage-05 Sammon map of 788 chains. B–C: Pair-level contingency tables (n=291 strict core-3 pairs) with Cramér's V reported at pair, unique-ligand, and DOI-deduplicated levels. Sources: results/05_maps/maps_coords.csv, chem_ligands_order.csv.",
        manifest, "Chemical clusters on structural map and contingency strength",
        "results/05_maps/maps_coords.csv; data/ligands/ligand_ids.csv; esr1_dataset_final.csv",
        "appendix", "chains (A), PDB×CCD pairs (B–C)"
    )

    # ========================================================================== F04: QC Clustering Sensitivity
    qc_df = cc.read_csv_safe(cc.TABLES_DIR / "qc_clustering_sensitivity.csv")
    qc_tr = cc.read_csv_safe(cc.TABLES_DIR / "qc_clustering_transitions_1p5.csv")

    fig, axs = plt.subplots(1, 3, figsize=(13.2, 4.4))
    ax = axs[0]
    for sub, col, mk in [("rfree_le_0p30", "#2a78d6", "o"), ("qc_all3_tierA", "#eb6834", "s")]:
        d = qc_df[qc_df["subset"] == sub]
        ax.plot(d["cut_A"], d["ari_vs_full"], color=col, marker=mk, lw=2.2, ms=7,
                label=f"{sub} (n={int(d['n_chains'].iloc[0])} ch)")
    ax.set_ylim(0.85, 1.02)
    ax.set_xticks([1.0, 1.5, 2.0, 2.5])
    ax.set_xlabel("Average-linkage RMSD cut (Å)")
    ax.set_ylabel("ARI vs full 788-chain clustering")
    ax.set_title("A. Clustering stability under QC filtering")
    ax.legend(loc="lower left", fontsize=8.5)
    cc.grid(ax)

    ax = axs[1]
    q15 = qc_df[qc_df["cut_A"] == 1.5].set_index("subset").loc[["full_788", "rfree_le_0p30", "qc_all3_tierA"]]
    x = np.arange(3)
    w = 0.35
    ax.bar(x - w / 2, q15["hill_q1_pdb"], width=w, color="#2a78d6", label="Hill q1 (PDB-wt)")
    ax.bar(x + w / 2, q15["shannon_bits_pdb"], width=w, color="#1b9e77", label="Shannon H (bits)")
    for i, r in enumerate(q15.itertuples()):
        ax.text(i - w / 2, r.hill_q1_pdb + 0.05, f"{r.hill_q1_pdb:.2f}", ha="center", fontsize=8)
        ax.text(i + w / 2, r.shannon_bits_pdb + 0.05, f"{r.shannon_bits_pdb:.2f}", ha="center", fontsize=8)
    ax.set_xticks(x, [f"Full\n(788 ch)", f"R-free≤0.30\n({int(q15.loc['rfree_le_0p30','n_chains'])} ch)", f"All 3 QC\n({int(q15.loc['qc_all3_tierA','n_chains'])} ch)"])
    ax.set_ylim(0, 2.6)
    ax.set_ylabel("Diversity metric (1.5 Å cut)")
    ax.set_title("B. Entropy & Hill q1 barely change with QC")
    ax.legend(loc="upper right", fontsize=8)
    cc.grid(ax)

    ax = axs[2]
    # Class composition of Tier A vs Full (strict core3 pairs)
    p_full = core_pairs["action_class"].value_counts(normalize=True).reindex(cc.CLASS_ORDER_CORE).fillna(0)
    p_ta = core_pairs[core_pairs["qc_label"] == "A"]["action_class"].value_counts(normalize=True).reindex(cc.CLASS_ORDER_CORE).fillna(0)
    comp_df = pd.DataFrame({"Full (n=291)": p_full * 100, "Tier A (n={})".format(int((core_pairs["qc_label"] == "A").sum())): p_ta * 100})
    bot = np.zeros(2)
    for cls in cc.CLASS_ORDER_CORE:
        vals = comp_df.loc[cls].values
        ax.bar([0, 1], vals, bottom=bot, color=cc.CLASS_COLOR[cls], width=0.5, label=cls)
        for xi, v in enumerate(vals):
            ax.text(xi, bot[xi] + v / 2, f"{v:.0f}%", ha="center", va="center", color="white", fontsize=8.5, fontweight="bold")
        bot += vals
    ax.set_xticks([0, 1], comp_df.columns)
    ax.set_ylabel("Share of core-3 pairs (%)")
    ax.set_title("C. Class mix preserved in QC Tier A")
    ax.legend(loc="upper right", fontsize=8)

    cc.wrap_title(fig, "Re-clustering after strict QC filtering preserves the two-group H12 structure (ARI = 0.998 for R-free ≤ 0.30; 0.982 for Tier A).", width=102)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F04_qc_sensitivity", qc_df,
        "Symmetrically sliced Stage-02 RMSD matrix re-clustered with average linkage at 1.0–2.5 Å across Full (788 chains, 351 PDB), known R-free ≤ 0.30 (767 chains, 21 excluded), and All-3-QC Tier A (544 chains, 244 excluded). Sources: results/02_lid_groups/rmsd.npy, data/structures/chain_table.csv.",
        manifest, "Sensitivity of H12 clustering and diversity to crystallographic QC",
        "results/02_lid_groups/rmsd.npy; data/structures/chain_table.csv",
        "backup", "chains (A–B), PDB×CCD pairs (C)"
    )

    # ========================================================================== F05: Priority 1 — Funnel, Transitions, Sources
    fig, axs = plt.subplots(1, 3, figsize=(13.8, 4.7), gridspec_kw={"width_ratios": [1.05, 1.15, 1.15]})
    # Panel A: Funnel
    ax = axs[0]
    n_geom = int(core_pg["dist_H12c_H11axis_A"].notna().sum())
    funnel_steps = [
        ("Real ligand pairs (CSV)", 368),
        ("3 core functional classes", 302),
        ("Exact primary CCD match", int(pairs[pairs["action_class"].isin(cc.CLASS_ORDER_CORE)].shape[0]) + 4),  # strict+multi check
        ("Strict single-ligand pairs", len(core_pairs)),
        ("With H12–H11 distance", int(core_pairs["dist_H12c_H11axis_A"].notna().sum())),
    ]
    # Exact primary match for core 3:
    n_prim_core3 = int(chains[chains["in_core3_with_multi"].astype(str) == "True"][["pdb_norm", "lig_norm"]].drop_duplicates().shape[0])
    funnel_steps[2] = ("Exact primary CCD match", n_prim_core3)
    y = np.arange(len(funnel_steps))
    vals = [v for _, v in funnel_steps]
    ax.barh(y, vals, color=["#8a8984", "#4a3aa7", "#2a78d6", "#1b9e77", "#0d6e52"], height=0.6)
    for i, (lbl, v) in enumerate(funnel_steps):
        ax.text(v + 5, i, f"{v} ({v/368:.0%})", va="center", fontsize=8.5, color=cc.INK)
    ax.set_yticks(y, [l for l, _ in funnel_steps], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 460)
    ax.set_xlabel("Unique PDB × CCD pairs (% of 368 real pairs)")
    ax.set_title("A. Coverage funnel (unit: PDB × CCD)")
    cc.grid(ax, "x")

    # Panel B: Old ChEMBL -> New Curated transition matrix on strict matched pairs (n=355)
    ax = axs[1]
    old_ord = ["agonist", "SERM", "SERD", "unlabelled"]
    new_ord = ["agonist", "partial_agonist", "antagonist", "unknown", "mixed"]
    tmat = pd.crosstab(pairs["old_chembl_class"], pairs["action_class"]).reindex(index=old_ord, columns=new_ord, fill_value=0)
    ax.imshow(tmat.values, cmap="Blues", aspect="auto")
    for i in range(tmat.shape[0]):
        for j in range(tmat.shape[1]):
            v = int(tmat.values[i, j])
            ax.text(j, i, str(v) if v > 0 else "·", ha="center", va="center", fontsize=8.5,
                    color="white" if v > tmat.values.max() * 0.5 else cc.INK)
    ax.set_xticks(range(len(new_ord)), [c.replace("_", "\n") for c in new_ord], fontsize=8)
    ax.set_yticks(range(len(old_ord)), [f"old: {o} ({int(tmat.loc[o].sum())})" for o in old_ord], fontsize=8)
    ax.set_title(f"B. Old ChEMBL → New curated (n={len(pairs)} strict pairs)")

    # Panel C: Primary evidence source by curated class
    ax = axs[2]
    src_ord = ["full_text", "abstract", "derived_from_data", "ligand_identity", "literature", "chembl", "pdb_metadata", "none"]
    smat = pd.crosstab(pairs["primary_source"], pairs["action_class"]).reindex(index=src_ord, columns=new_ord, fill_value=0)
    smat = smat[smat.sum(axis=1) > 0]
    bot = np.zeros(len(smat))
    y_s = np.arange(len(smat))
    for cls in new_ord:
        v = smat[cls].values
        ax.barh(y_s, v, left=bot, color=cc.CLASS_COLOR[cls], height=0.65, label=cls)
        bot += v
    for i, tot in enumerate(bot):
        ax.text(tot + 2, i, str(int(tot)), va="center", fontsize=8, color=cc.INK)
    ax.set_yticks(y_s, smat.index, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, bot.max() * 1.18)
    ax.set_xlabel("Strict single-ligand pairs (n=355)")
    ax.set_title("C. Evidence source by curated class")
    ax.legend(loc="lower right", fontsize=7.5)
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Manual curation increases functional coverage nearly 6-fold (50 → 291 strict core-3 pairs) and introduces partial_agonist (n=78).", width=105)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F05_annotation_funnel_transitions", tmat.reset_index(),
        "A: Pair-level attrition from 368 real ligand pairs in esr1_dataset_final.csv to 291 strict single-ligand core-3 pairs. B: Transition matrix from historical ChEMBL mechanism labels (data/ligands/chain_class.csv) to curated action_class on the 355 strict matched pairs. C: Primary evidence source breakdown across the 355 strict pairs (unknown=62 pairs, mostly primary_source=none/abstract without functional assay).",
        manifest, "Annotation gain, old-to-new label transitions, and evidence sources",
        "esr1_dataset_final.csv; data/structures/chain_table.csv; data/ligands/chain_class.csv",
        "defense", "PDB × CCD pairs"
    )

    # ========================================================================== F06: Priority 2 — Functional Class x H12 Group
    fig, axs = plt.subplots(1, 3, figsize=(13.6, 4.6), gridspec_kw={"width_ratios": [1.2, 1.1, 1.0]})
    # Panel A: P(H12 group | class) — mutually exclusive states on core3 + unknown
    ax = axs[0]
    cls_show = ["agonist", "partial_agonist", "antagonist", "unknown"]
    grp_states = ["cl7", "cl4", "minor", "multiple"]
    ct_gc = pd.crosstab(pairs["action_class"], pairs["h12_group_excl"]).reindex(index=cls_show, columns=grp_states, fill_value=0)
    pct_gc = ct_gc.div(ct_gc.sum(axis=1), axis=0) * 100
    y = np.arange(len(cls_show))
    left = np.zeros(len(cls_show))
    for g in grp_states:
        v = pct_gc[g].values
        ax.barh(y, v, left=left, color=cc.GROUP_COLOR[g], height=0.58, label=cc.GROUP_LABEL[g])
        for i, (pv, cnt) in enumerate(zip(v, ct_gc[g].values)):
            if pv >= 7:
                ax.text(left[i] + pv / 2, i, f"{cnt}\n({pv:.0f}%)", ha="center", va="center",
                        color="white" if g != "minor" else cc.INK, fontsize=8, fontweight="bold")
        left += v
    ax.set_yticks(y, [f"{c}\n(n={int(ct_gc.loc[c].sum())})" for c in cls_show], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of PDB × CCD pairs (%)")
    ax.set_title("A. P(H12 group | class) — mutually exclusive states")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7.8)

    # Panel B: P(class | H12 group) on core3 pairs
    ax = axs[1]
    ct_cg = pd.crosstab(core_pairs["h12_group_excl"], core_pairs["action_class"]).reindex(index=["cl7", "cl4", "minor", "multiple"], columns=cc.CLASS_ORDER_CORE, fill_value=0)
    pct_cg = ct_cg.div(ct_cg.sum(axis=1).replace(0, 1), axis=0) * 100
    y2 = np.arange(4)
    left = np.zeros(4)
    for c in cc.CLASS_ORDER_CORE:
        v = pct_cg[c].values
        ax.barh(y2, v, left=left, color=cc.CLASS_COLOR[c], height=0.58, label=c)
        for i, (pv, cnt) in enumerate(zip(v, ct_cg[c].values)):
            if pv >= 8:
                ax.text(left[i] + pv / 2, i, f"{cnt}\n({pv:.0f}%)", ha="center", va="center",
                        color="white", fontsize=8, fontweight="bold")
        left += v
    ax.set_yticks(y2, [f"{g}\n(n={int(ct_cg.loc[g].sum())})" for g in ct_cg.index], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Share of core-3 pairs in group (%)")
    ax.set_title("B. P(class | H12 group) on core-3 (n=291)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7.8)

    # Panel C: Historical L7 comparison (Old ChEMBL vs New Curated on shared 50-pair subset vs full 291)
    ax = axs[2]
    old_lab = pairs[pairs["old_chembl_class"].isin(["agonist", "SERM", "SERD"])]
    comp_rows = [
        ("Old ChEMBL\nagonist", (old_lab["old_chembl_class"] == "agonist").sum(), ((old_lab["old_chembl_class"] == "agonist") & (old_lab["h12_group_excl"] == "cl7")).sum()),
        ("New curated\nagonist", (core_pairs["action_class"] == "agonist").sum(), ((core_pairs["action_class"] == "agonist") & (core_pairs["h12_group_excl"] == "cl7")).sum()),
        ("Old ChEMBL\nSERM+SERD", (old_lab["old_chembl_class"].isin(["SERM", "SERD"])).sum(), ((old_lab["old_chembl_class"].isin(["SERM", "SERD"])) & (old_lab["h12_group_excl"] == "cl4")).sum()),
        ("New curated\nantagonist", (core_pairs["action_class"] == "antagonist").sum(), ((core_pairs["action_class"] == "antagonist") & (core_pairs["h12_group_excl"] == "cl4")).sum()),
    ]
    y3 = np.arange(len(comp_rows))
    shares = [k / n * 100 if n else 0 for _, n, k in comp_rows]
    ax.barh(y3, shares, color=[cc.CL7, cc.CL7, cc.CL4, cc.CL4], alpha=[0.45, 0.95, 0.45, 0.95][0], height=0.55)
    for idx_b, bar in enumerate(ax.patches):
        bar.set_alpha(0.45 if idx_b in (0, 2) else 0.95)
    for i, (lbl, n, k) in enumerate(comp_rows):
        ax.text(shares[i] - 3 if shares[i] > 25 else shares[i] + 2, i, f"{k}/{n} ({shares[i]:.0f}%)",
                va="center", ha="right" if shares[i] > 25 else "left",
                color="white" if shares[i] > 25 else cc.INK, fontsize=8.2, fontweight="bold")
    ax.set_yticks(y3, [r[0] for r in comp_rows], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 105)
    ax.set_xlabel("Share in canonical H12 group (%)")
    ax.set_title("C. Historical L7 vs curated mapping")
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Partial agonists split across cl7 (83%) and cl4 (17%), while 22 antagonist pairs (19%) sit in cl7—breaking the near-1:1 ChEMBL L7 lock.", width=106)
    fig.tight_layout(rect=[0, 0.12, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F06_class_vs_h12_group", ct_gc.reset_index(),
        "A: Mutually exclusive H12-group distribution by curated class (strict PDB×CCD pairs; n=353 excluding 2 mixed). B: Functional class composition within each H12 group (n=291 core-3 pairs). C: Comparison with historical claim L7 (ChEMBL n=50 strict pairs vs curated n=291). Sources: pairs_curated.csv.",
        manifest, "Distribution of curated functional classes across H12 conformations",
        "esr1_dataset_final.csv; data/structures/chain_table.csv; data/ligands/chain_class.csv",
        "defense", "PDB × CCD pairs"
    )

    # ========================================================================== F07: Priority 3 — H12 within-group spread & forest plot
    fig, axs = plt.subplots(1, 3, figsize=(13.8, 4.8), gridspec_kw={"width_ratios": [1.05, 1.05, 1.25]})
    rng = np.random.default_rng(cc.SEED)
    for ax, grp, thr in [(axs[0], "cl7", 0.91), (axs[1], "cl4", 1.08)]:
        sub = core_pg[(core_pg["group"] == grp) & core_pg["dist_H12c_H11axis_A"].notna()]
        med_all = float(sub["dist_H12c_H11axis_A"].median())
        ax.axhspan(med_all - thr, med_all + thr, color=cc.GROUP_COLOR[grp], alpha=0.10, label=f"ordinary ±{thr:.2f} Å")
        ax.axhline(med_all, color=cc.GROUP_COLOR[grp], ls="--", lw=1)
        xt_labels = []
        for i, cls in enumerate(cc.CLASS_ORDER_CORE):
            d = sub[sub["action_class"] == cls]["dist_H12c_H11axis_A"].to_numpy()
            if len(d) == 0:
                xt_labels.append(cls.replace("_", "\n") + "\n(n=0)")
                continue
            xj = i + rng.uniform(-0.18, 0.18, size=len(d))
            ax.scatter(xj, d, color=cc.CLASS_COLOR[cls], marker=cc.CLASS_MARKER[cls], s=26, alpha=0.72, edgecolors="white", lw=0.4)
            q25, q50, q75 = np.percentile(d, [25, 50, 75])
            ax.plot([i - 0.25, i + 0.25], [q50, q50], color=cc.INK, lw=2.3, zorder=5)
            ax.plot([i, i], [q25, q75], color=cc.INK, lw=2.0, zorder=5)
            xt_labels.append(f"{cls.replace('_', ' ')}\nn={len(d)}, med={q50:.2f}")
        ax.set_xticks(range(3), xt_labels, fontsize=8.0)
        ax.set_ylabel("dist_H12c_H11axis_A (Å, pair×group median)")
        ax.set_title(f"{'A' if grp=='cl7' else 'B'}. Inside {grp} (ordinary ±{thr:.2f} Å)")
        cc.grid(ax)

    # Panel C: Forest plot of pre-registered contrasts (with 95% DOI-block bootstrap CI)
    ax = axs[2]
    fsub = eff_df[eff_df["metric"] == "dist_H12c_H11axis_A"].head(8).copy()
    y = np.arange(len(fsub))
    ax.axvspan(-0.91, 0.91, color=cc.CL7, alpha=0.08, label="cl7 ordinary ±0.91 Å")
    ax.axvline(0, color=cc.INK2, lw=0.9)
    for i, r in enumerate(fsub.itertuples()):
        col = cc.CL7 if r.group == "cl7" else cc.CL4
        ax.scatter(r.diff_a_minus_b, i, color=col, s=50, zorder=4)
        if r.ci_valid:
            ax.plot([r.ci95_lo, r.ci95_hi], [i, i], color=col, lw=2.0, zorder=3)
        lbl_txt = f"{r.diff_a_minus_b:+.2f} Å ({r.n_a}v{r.n_b}, {r.shared_dois} shDOI)"
        ax.text(max(r.ci95_hi if r.ci_valid else r.diff_a_minus_b, r.diff_a_minus_b) + 0.05, i, lbl_txt,
                va="center", fontsize=7.2, color=cc.INK)
    ax.set_yticks(y, fsub["contrast_id"], fontsize=7.8)
    ax.invert_yaxis()
    ax.set_xlim(-1.2, 1.85)
    ax.set_xlabel("Median difference Δ dist_H12c_H11axis (Å, 95% DOI-block CI)")
    ax.set_title("C. Comparable-stratum vs crude effects")
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Inside cl7 and cl4, functional class shifts H12 by ≤ 0.25 Å — well inside ordinary crystal-copy variation (0.91 / 1.08 Å).", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F07_h12_within_group_spread", fsub,
        "A–B: Pair×group medians of dist_H12c_H11axis_A (Å) with median/IQR bars and shaded ordinary variation band (±0.91 Å cl7, ±1.08 Å cl4). C: Forest plot of median differences with 95% DOI-block bootstrap CIs (2000 resamples; shown when ≥4 DOI blocks per group) and shared-DOI counts. Sources: pairs_by_group.csv, effects_summary.csv.",
        manifest, "Within-conformation H12 shift by functional class in comparable strata",
        "data/lid/geometry_per_chain.csv; esr1_dataset_final.csv",
        "defense", "PDB × CCD × H12-group"
    )

    # ========================================================================== F08: Priority 4 — Confounding, Support, LODO, Within-DOI
    fig, axs = plt.subplots(1, 3, figsize=(13.8, 4.8), gridspec_kw={"width_ratios": [1.25, 1.0, 1.05]})
    # Panel A: Support matrix (H12 group x mutation_status x peptide_bound) x action_class
    ax = axs[0]
    core_pg["stratum_lbl"] = core_pg["group"] + " | " + core_pg["mutation_status"].str.replace("Cys->Ser/Ala only", "3xCys") + " | " + np.where(core_pg["peptide_bound"].astype(str) == "True", "pep+", "pep-")
    stab_s = pd.crosstab(core_pg["stratum_lbl"], core_pg["action_class"]).reindex(columns=cc.CLASS_ORDER_CORE, fill_value=0)
    stab_s["tot"] = stab_s.sum(axis=1)
    stab_s = stab_s[stab_s["tot"] >= 3].sort_values("tot", ascending=False).drop(columns="tot")
    ax.imshow(stab_s.values, cmap="YlOrRd", aspect="auto")
    for i in range(stab_s.shape[0]):
        for j in range(stab_s.shape[1]):
            v = int(stab_s.values[i, j])
            ax.text(j, i, str(v) if v > 0 else "·", ha="center", va="center", fontsize=8.5,
                    color="white" if v > stab_s.values.max() * 0.55 else cc.INK)
    ax.set_xticks(range(3), [c.replace("_", "\n") for c in cc.CLASS_ORDER_CORE], fontsize=8.5)
    ax.set_yticks(range(len(stab_s)), stab_s.index, fontsize=7.8)
    ax.set_title("A. Support: H12 × construct × peptide (n≥3)")

    # Panel B: LODO for cl7 Y537S+pep (ant vs ago & part vs ago)
    ax = axs[1]
    for cname, col, mk in [("ant_vs_ago", cc.CLASS_COLOR["antagonist"], "D"), ("part_vs_ago", cc.CLASS_COLOR["partial_agonist"], "s")]:
        d = lodo_df[lodo_df["contrast"] == cname].reset_index(drop=True)
        y = np.arange(len(d))
        ax.scatter(d["lodo_diff"], y + (0.15 if cname == "part_vs_ago" else -0.15), color=col, marker=mk, s=42, label=f"LODO {cname}")
        ax.axvline(d["full_diff"].iloc[0], color=col, ls="--", lw=1.1, alpha=0.7)
    ax.axvspan(-0.91, 0.91, color=cc.CL7, alpha=0.07)
    ax.axvline(0, color=cc.MUTED, lw=0.8)
    d0 = lodo_df[lodo_df["contrast"] == "ant_vs_ago"].reset_index(drop=True)
    ax.set_yticks(np.arange(len(d0)), [s[-22:] for s in d0["doi_block"]], fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlabel("Leave-one-DOI-out median Δ (Å)")
    ax.set_title("B. Leave-one-DOI-out (cl7 Y537S + pep)")
    ax.legend(loc="lower right", fontsize=7.5)
    cc.grid(ax, "x")

    # Panel C: Within-DOI direct comparisons where both classes exist in the same DOI
    ax = axs[2]
    wdoi = lodo_df[lodo_df["within_doi_diff"].notna()].copy()
    if len(wdoi):
        y = np.arange(len(wdoi))
        cols = [cc.CLASS_COLOR["antagonist"] if c == "ant_vs_ago" else cc.CLASS_COLOR["partial_agonist"] for c in wdoi["contrast"]]
        ax.barh(y, wdoi["within_doi_diff"], color=cols, height=0.55, alpha=0.85)
        ax.axvline(0, color=cc.INK, lw=0.9)
        ax.axvspan(-0.91, 0.91, color=cc.CL7, alpha=0.08)
        for i, r in enumerate(wdoi.itertuples()):
            ax.text(r.within_doi_diff + (0.03 if r.within_doi_diff >= 0 else -0.03), i,
                    f"{r.within_doi_diff:+.2f} Å ({r.n_in_doi_a}v{r.n_in_doi_b})",
                    va="center", ha="left" if r.within_doi_diff >= 0 else "right", fontsize=7.5)
        ax.set_yticks(y, [f"{r.contrast[:4]} | {r.doi_block[-18:]}" for r in wdoi.itertuples()], fontsize=7.3)
        ax.invert_yaxis()
        ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Within-DOI median Δ dist_H12c_H11axis (Å)")
    ax.set_title("C. Direct within-DOI pairs (both classes)")
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Only one stratum (cl7 Y537S + peptide) has real 3-class overlap, and both LODO and within-DOI shifts stay inside ±0.91 Å.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F08_confounding_support_lodo", lodo_df,
        "A: Pair×group support across (H12 group × mutation_status × peptide_bound) strata with n≥3. B: Leave-one-DOI-out median shifts in cl7 Y537S+peptide. C: Within-DOI shifts restricted to DOIs containing both compared classes. Shaded band: ±0.91 Å cl7 ordinary variation. Sources: pairs_by_group.csv, lodo_within_doi.csv.",
        manifest, "Confounding by construct/peptide, LODO stability, and within-DOI contrasts",
        "data/structures/chain_table.csv; data/lid/geometry_per_chain.csv; esr1_dataset_final.csv",
        "defense", "PDB × CCD × H12-group"
    )

    # ========================================================================== F09: Priority 5 — H12 Geometry, Order & Packing
    fig, axs = plt.subplots(2, 2, figsize=(12.2, 8.2))
    # A: H12-H3 vs H12-H11 angles colored by class, shaped by group
    ax = axs[0, 0]
    d9 = core_pg.dropna(subset=["angle_H12_H3_deg", "angle_H12_H11_deg"])
    for cls in cc.CLASS_ORDER_CORE:
        for grp in ("cl7", "cl4"):
            s = d9[(d9["action_class"] == cls) & (d9["group"] == grp)]
            if len(s):
                ax.scatter(s["angle_H12_H3_deg"], s["angle_H12_H11_deg"],
                           c=cc.CLASS_COLOR[cls], marker=cc.GROUP_MARKER[grp], s=32, alpha=0.8,
                           label=f"{cls} | {grp} ({len(s)})")
    ax.set_xlabel("H12–H3 angle (°)")
    ax.set_ylabel("H12–H11 angle (°)")
    ax.set_title("A. H12 orientation angles by class & group")
    ax.legend(loc="best", fontsize=7.5, ncol=2)
    cc.grid(ax, "both")

    # B: H12 completeness in 538-547 vs 538-544 windows with 95% Wilson CI
    ax = axs[0, 1]
    comp_rows = []
    for grp in ("cl7", "cl4"):
        for cls in cc.CLASS_ORDER_CORE:
            s = core_pg[(core_pg["group"] == grp) & (core_pg["action_class"] == cls)]
            if len(s) < 3:
                continue
            for win_col, win_lbl in [("h12_complete_547", "538–547"), ("h12_complete_544", "538–544")]:
                k = int((s[win_col].astype(str) == "True").sum())
                n = len(s)
                lo, hi = cc.wilson_ci(k, n)
                comp_rows.append(dict(grp=grp, cls=cls, win=win_lbl, k=k, n=n, p=k / n, lo=lo, hi=hi))
    cdf = pd.DataFrame(comp_rows)
    cats = cdf[["grp", "cls"]].drop_duplicates().reset_index(drop=True)
    x = np.arange(len(cats))
    for wi, (wlbl, col, off) in enumerate([("538–547", "#4a3aa7", -0.18), ("538–544", "#1b9e77", 0.18)]):
        sub_w = cdf[cdf["win"] == wlbl].set_index(["grp", "cls"])
        ps = [sub_w.loc[(r.grp, r.cls), "p"] * 100 for r in cats.itertuples()]
        los = [max(0.0, (sub_w.loc[(r.grp, r.cls), "p"] - sub_w.loc[(r.grp, r.cls), "lo"]) * 100) for r in cats.itertuples()]
        his = [max(0.0, (sub_w.loc[(r.grp, r.cls), "hi"] - sub_w.loc[(r.grp, r.cls), "p"]) * 100) for r in cats.itertuples()]
        ax.bar(x + off, ps, width=0.34, color=col, alpha=0.85, label=f"Window {wlbl}")
        ax.errorbar(x + off, ps, yerr=[los, his], fmt="none", ecolor=cc.INK, elinewidth=1.1, capsize=3)
    ax.set_xticks(x, [f"{r.grp}\n{r.cls[:7]}" for r in cats.itertuples()], fontsize=8)
    ax.set_ylim(0, 112)
    ax.set_ylabel("Complete H12 pairs (%, 95% Wilson CI)")
    ax.set_title("B. H12 completeness: 538–547 vs 538–544 window")
    ax.legend(loc="lower left", fontsize=8)
    cc.grid(ax)

    # C: mobility_h12_z by class & group
    ax = axs[1, 0]
    dm = core_pg.dropna(subset=["mobility_h12_z"])
    pos_x = 0
    xt, xtl = [], []
    for grp in ("cl7", "cl4"):
        for cls in cc.CLASS_ORDER_CORE:
            v = dm[(dm["group"] == grp) & (dm["action_class"] == cls)]["mobility_h12_z"].to_numpy()
            if len(v) < 3:
                continue
            xj = pos_x + rng.uniform(-0.18, 0.18, size=len(v))
            ax.scatter(xj, v, color=cc.CLASS_COLOR[cls], marker=cc.GROUP_MARKER[grp], s=24, alpha=0.72)
            med = float(np.median(v))
            ax.plot([pos_x - 0.25, pos_x + 0.25], [med, med], color=cc.INK, lw=2.2)
            xt.append(pos_x)
            xtl.append(f"{grp}\n{cls[:6]}\n(n={len(v)})")
            pos_x += 1
        pos_x += 0.5
    ax.axhline(0, color=cc.MUTED, ls="--", lw=0.9)
    ax.set_xticks(xt, xtl, fontsize=7.8)
    ax.set_ylabel("mobility_h12_z (normalized B-factor z-score)")
    ax.set_title("C. Crystallographic H12 mobility z-score")
    cc.grid(ax)

    # D: mobility_h12_z vs sym_atoms & resolution (model artifact check)
    ax = axs[1, 1]
    d_sym = core_pg.dropna(subset=["mobility_h12_z", "sym_atoms", "resolution_A"])
    for grp in ("cl7", "cl4"):
        s = d_sym[d_sym["group"] == grp]
        rho_s, p_s = spearmanr(s["sym_atoms"], s["mobility_h12_z"])
        rho_r, p_r = spearmanr(s["resolution_A"], s["mobility_h12_z"])
        ax.scatter(s["sym_atoms"], s["mobility_h12_z"], c=cc.GROUP_COLOR[grp], s=26, alpha=0.7,
                   label=f"{grp}: ρ(sym)={rho_s:.2f}, ρ(res)={rho_r:.2f}")
    ax.set_xlabel("Symmetry-mate contact atoms at H12 (528–550)")
    ax.set_ylabel("mobility_h12_z")
    ax.set_title("D. Mobility tracks crystal packing & refinement")
    ax.legend(loc="upper right", fontsize=8)
    cc.grid(ax, "both")

    cc.wrap_title(fig, "H12 orientation, completeness, and B-factor mobility track conformation and crystal packing rather than functional class.", width=104)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    cc.save_fig_bundle(
        fig, "F09_h12_geometry_order_packing", cdf,
        "A: H12–H3 vs H12–H11 angles across core-3 pair×group entries. B: Share of complete H12 models in 538–547 vs consensus 538–544 windows with 95% Wilson CIs. C: Crystallographic H12 B-factor z-score (mobility_h12_z). D: Association of mobility_h12_z with symmetry-mate contacts and resolution. Sources: geometry_per_chain.csv, mobility_per_chain.csv, lattice_contacts_per_chain.csv.",
        manifest, "H12 geometry, completeness windows, mobility, and crystal packing",
        "data/lid/geometry_per_chain.csv; data/lid/mobility_per_chain.csv; data/lid/lattice_contacts_per_chain.csv",
        "appendix", "PDB × CCD × H12-group"
    )

    # ========================================================================== F10: Priority 6 — Same Ligand & Same Crystal
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 5.2))
    # Panel A: Repeated InChIKeys across different constructs or peptide states
    ax = axs[0]
    ik_counts = core_pg.groupby("inchikey").agg(
        n_rows=("pdb_norm", "size"),
        n_mut=("mutation_status", "nunique"),
        n_pep=("peptide_bound", "nunique"),
        n_grp=("group", "nunique"),
        rng_d=("dist_H12c_H11axis_A", lambda s: s.max() - s.min()),
        ccd=("lig_norm", "first"),
    ).reset_index()
    rep_iks = ik_counts[(ik_counts["n_rows"] >= 2) & ((ik_counts["n_mut"] > 1) | (ik_counts["n_pep"] > 1) | (ik_counts["n_grp"] > 1))].sort_values("rng_d", ascending=False).head(12)
    rep_sub = core_pg[core_pg["inchikey"].isin(rep_iks["inchikey"])].copy()
    for i, r in enumerate(rep_iks.itertuples()):
        s = rep_sub[rep_sub["inchikey"] == r.inchikey].dropna(subset=["dist_H12c_H11axis_A"])
        if len(s) >= 2:
            ax.plot([s["dist_H12c_H11axis_A"].min(), s["dist_H12c_H11axis_A"].max()], [i, i], color=cc.MUTED, lw=1.6, zorder=1)
        for _, pt in s.iterrows():
            mk = "o" if str(pt["peptide_bound"]) == "True" else "s"
            ax.scatter(pt["dist_H12c_H11axis_A"], i, c=cc.CLASS_COLOR[pt["action_class"]], marker=mk, s=55, edgecolors=cc.INK, lw=0.6, zorder=3)
    ax.set_yticks(range(len(rep_iks)), [f"{r.ccd} (n={r.n_rows})" for r in rep_iks.itertuples()], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("dist_H12c_H11axis_A (Å)  [circle=pep+, square=pep−]")
    ax.set_title("A. Same ligand (InChIKey) across constructs / peptide")
    cc.grid(ax, "x")

    # Panel B: Same crystal (PDB) multi-chain divergence in dist_H12c_H11axis_A
    ax = axs[1]
    ch_c = chains[chains["in_strict_single"].astype(str) == "True"].dropna(subset=["dist_H12c_H11axis_A"]).copy()
    crys = ch_c.groupby("pdb_norm").agg(
        n_ch=("chain_key", "size"),
        min_d=("dist_H12c_H11axis_A", "min"),
        max_d=("dist_H12c_H11axis_A", "max"),
        spread=("dist_H12c_H11axis_A", lambda s: s.max() - s.min()),
        lig=("lig_norm", "first"),
        cls=("action_class", "first"),
    ).reset_index()
    top_crys = crys[crys["n_ch"] >= 2].sort_values("spread", ascending=False).head(12)
    for i, r in enumerate(top_crys.itertuples()):
        s = ch_c[ch_c["pdb_norm"] == r.pdb_norm]
        ax.plot([r.min_d, r.max_d], [i, i], color=cc.INK2, lw=1.6, zorder=1)
        for _, pt in s.iterrows():
            ax.scatter(pt["dist_H12c_H11axis_A"], i, c=cc.GROUP_COLOR[pt["group"]], s=48, edgecolors="white", lw=0.5, zorder=3)
        ax.text(r.max_d + 0.2, i, f"Δ={r.spread:.2f} Å ({r.cls[:4]})", va="center", fontsize=7.6, color=cc.INK2)
    ax.set_yticks(range(len(top_crys)), [f"{r.pdb_norm} ({r.lig})" for r in top_crys.itertuples()], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(10.5, 27.5)
    ax.set_xlabel("dist_H12c_H11axis_A (Å) across chains of one crystal")
    ax.set_title("B. Top within-crystal chain divergences (same PDB)")
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Identical ligands can switch H12 state with construct/peptide (e.g., GW5, EST), and copies in one crystal diverge up to 8.4 Å (9B25).", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F10_same_ligand_same_crystal", top_crys,
        "A: Top 12 repeated InChIKeys spanning multiple constructs, peptide states, or H12 groups (circles=peptide bound, squares=no peptide; colored by curated action_class). B: Top 12 crystals with largest within-PDB spread of dist_H12c_H11axis_A across chains (colored by H12 group). Sources: pairs_by_group.csv, chains_curated.csv.",
        manifest, "Same-ligand context switching and within-crystal chain divergence",
        "data/lid/geometry_per_chain.csv; esr1_dataset_final.csv",
        "backup", "PDB×CCD×group (A), chains within PDB (B)"
    )

    # ========================================================================== F11: Priority 7 — Outliers Audit + 3D C-alpha H11/H12 Superposition
    out_df = cc.read_csv_safe(cc.RESULTS_BASE / "03_lid_outliers/lid_outliers.csv")
    out_cur = out_df.merge(
        chains[["chain_key", "action_class", "doi", "primary_source", "drug_classes"]],
        on="chain_key", how="left", suffixes=("_old", "_cur")
    )
    out_cur.to_csv(cc.TABLES_DIR / "outliers_curated_audit.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(13.2, 5.4), gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax = axs[0]
    oplot = out_cur.dropna(subset=["lid_to_h11_A"]).sort_values("lid_to_h11_A").reset_index(drop=True)
    ax.axvspan(12.48 - 0.91, 12.48 + 0.91, color=cc.CL7, alpha=0.12, label="cl7 ordinary (12.48±0.91 Å)")
    ax.axvspan(22.26 - 1.08, 22.26 + 1.08, color=cc.CL4, alpha=0.12, label="cl4 ordinary (22.26±1.08 Å)")
    for i, r in enumerate(oplot.itertuples()):
        cls = r.action_class_cur if pd.notna(r.action_class_cur) else "unlabelled"
        col = cc.CLASS_COLOR.get(cls, cc.MUTED)
        # Plot crystal mates as hollow grey circles
        mates = chains[(chains["pdb_norm"] == r.entry_id) & (chains["chain_key"] != r.chain_key)]["dist_H12c_H11axis_A"].dropna()
        if len(mates):
            ax.scatter(mates, [i] * len(mates), facecolors="white", edgecolors=cc.MUTED, s=34, lw=1.1, zorder=2)
        ax.scatter(r.lid_to_h11_A, i, color=col, marker=cc.GROUP_MARKER.get(r.group, "o"), s=62, edgecolors=cc.INK, lw=0.6, zorder=4)
    ax.set_yticks(range(len(oplot)),
                  [f"{r.entry_id} {r.chain_key.split('|')[1]} ({r.ligand}, {str(r.action_class_cur)[:4]})" for r in oplot.itertuples()],
                  fontsize=7.6)
    ax.set_xlabel("dist_H12c_H11axis_A (Å)  [hollow grey = other chains in same crystal]")
    ax.set_title("A. Lid outliers (5KCU B, 5KD9 B, cl4 extremes, minor, 9B25)")
    ax.legend(loc="lower right", fontsize=7.6)
    cc.grid(ax, "x")

    # Panel B: Reproducible C-alpha trace superposition of H11 (497-526) & H12 (528-547) for key outliers vs reference 7NFB B (cl7) and 1ERR A (cl4)
    ax = axs[1]
    z_ca = np.load(cc.RESULTS_BASE / "02_lid_groups/aligned_ca.npz", allow_pickle=True)
    ca_keys = [str(k) for k in z_ca["chain_keys"]]
    ca_pos = z_ca["positions"].astype(int)
    ca_xyz = z_ca["coords"].astype(float)
    m_h11 = (ca_pos >= 497) & (ca_pos <= 526)
    m_h12 = (ca_pos >= 528) & (ca_pos <= 547)
    show_chains = [
        ("7NFB|B|1", "7NFB B (cl7 ref, agonist)", cc.CL7, "-", 2.0),
        ("1ERR|A|1", "1ERR A (cl4 ref, RAL ant)", cc.CL4, "-", 2.0),
        ("5KCU|B|1", "5KCU B (cl7 outlier +1.53 Å)", "#1b9e77", "--", 2.2),
        ("5KD9|B|1", "5KD9 B (cl7 outlier +1.65 Å)", "#7b3294", "--", 2.2),
        ("9B25|A|1", "9B25 A (minor intermediate)", "#c2185b", "-", 2.6),
        ("9B25|B|1", "9B25 B (minor intermediate)", "#5d4037", ":", 2.4),
    ]
    # Project onto 2D plane (PC1-PC2 of H11+H12 coordinates across shown chains)
    pts_all = []
    for ck, *_ in show_chains:
        if ck in ca_keys:
            xyz = ca_xyz[ca_keys.index(ck)][m_h11 | m_h12]
            pts_all.append(xyz[np.isfinite(xyz[:, 0])])
    P = np.vstack(pts_all)
    mu = P.mean(0)
    _, _, vt = np.linalg.svd(P - mu, full_matrices=False)
    basis = vt[:2].T  # (3, 2)
    for ck, lbl, col, ls, lw in show_chains:
        if ck not in ca_keys:
            continue
        xyz = ca_xyz[ca_keys.index(ck)]
        h11_2d = (xyz[m_h11] - mu) @ basis
        h12_2d = (xyz[m_h12] - mu) @ basis
        ax.plot(h11_2d[:, 0], h11_2d[:, 1], color="#cfd0cc", lw=1.2, alpha=0.7)
        ax.plot(h12_2d[:, 0], h12_2d[:, 1], color=col, ls=ls, lw=lw, label=lbl)
        ok = np.isfinite(h12_2d[:, 0])
        if ok.any():
            end_pt = h12_2d[ok][-1]
            ax.scatter(end_pt[0], end_pt[1], color=col, s=30, zorder=5)
    ax.set_xlabel("H11–H12 principal axis 1 (Å, core-v3 aligned)")
    ax.set_ylabel("H11–H12 principal axis 2 (Å)")
    ax.set_title("B. Aligned Cα superposition: H11 (grey) & H12 (colored)")
    ax.legend(loc="best", fontsize=7.4)
    cc.grid(ax, "both")

    cc.wrap_title(fig, "Outlier audit: 5KCU B and 5KD9 B are single-chain cl7 shifts (+1.5–1.7 Å), while 9B25 holds intermediate and groove-like H12 copies in one crystal.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F11_outliers_audit", oplot[["chain_key", "entry_id", "group", "reason_flagged", "lid_to_h11_A", "ligand", "action_class_cur", "doi", "primary_source"]],
        "A: Examined lid outliers from Stage 03 annotated with curated action_class (filled) and their crystal mates (hollow grey circles). B: Core-v3 aligned Cα trace of H11 (497–526, grey) and H11-H12/H12 (528–547, colored) for references and key outliers. Sources: results/03_lid_outliers/lid_outliers.csv, results/02_lid_groups/aligned_ca.npz, esr1_dataset_final.csv.",
        manifest, "Audit of lid outliers (5KCU B, 5KD9 B, cl4 extremes, minor, 9B25)",
        "results/03_lid_outliers/lid_outliers.csv; results/02_lid_groups/aligned_ca.npz; esr1_dataset_final.csv",
        "defense", "chains"
    )

    # ========================================================================== F12: Priority 8 — Structure <-> Efficacy
    eff_sub = pairs.dropna(subset=["luc_efficacy_pct", "proliferation_pct"]).copy()
    eff_pg = pg.dropna(subset=["luc_efficacy_pct", "proliferation_pct"]).copy()
    top_dois_eff = eff_sub["doi"].value_counts().index.tolist()

    fig, axs = plt.subplots(2, 3, figsize=(13.5, 8.2))
    # Top row + bottom-left: faceted by DOI (up to 4 panels) + bottom-middle/right: geometry vs continuous efficacy
    for idx_d in range(min(4, len(top_dois_eff))):
        ax = axs[idx_d // 3, idx_d % 3]
        doi_val = top_dois_eff[idx_d]
        d = eff_sub[eff_sub["doi"] == doi_val]
        # Rule boundaries from DATA_DICTIONARY.md (for explanation only)
        ax.axvline(80, color=cc.MUTED, ls=":", lw=1)
        ax.axvline(20, color=cc.MUTED, ls=":", lw=1)
        ax.axhline(0, color=cc.MUTED, ls=":", lw=1)
        for cls in cc.CLASS_ORDER_ALL:
            s = d[d["action_class"] == cls]
            if len(s):
                ax.scatter(s["luc_efficacy_pct"], s["proliferation_pct"], c=cc.CLASS_COLOR[cls],
                           marker=cc.CLASS_MARKER[cls], s=40, alpha=0.85, edgecolors="white", lw=0.4, label=f"{cls} ({len(s)})")
        ax.set_xlabel("luc_efficacy_pct (% of E2)")
        ax.set_ylabel("proliferation_pct (% of E2)")
        n_der = int((d["primary_source"] == "derived_from_data").sum())
        ax.set_title(f"DOI: {doi_val[:26]} (n={len(d)}, derived={n_der})", fontsize=8.8)
        ax.legend(loc="best", fontsize=7.2)
        cc.grid(ax, "both")

    # Bottom-middle: dist_H12c_H11axis_A vs luc_efficacy_pct inside cl7 (highlighting WT-assay on mutant)
    ax = axs[1, 1]
    e7 = eff_pg[(eff_pg["group"] == "cl7") & eff_pg["dist_H12c_H11axis_A"].notna()]
    for is_wtm, mk, ec, lbl_m in [(False, "o", "white", "matched context"), (True, "X", cc.INK, "WT assay / mutant struct")]:
        s = e7[e7["wt_assay_on_mutant_note"].astype(str).eq(str(is_wtm))]
        if len(s):
            ax.scatter(s["luc_efficacy_pct"], s["dist_H12c_H11axis_A"],
                       c=[cc.CLASS_COLOR[c] for c in s["action_class"]], marker=mk, s=46, edgecolors=ec, lw=0.6, alpha=0.85, label=f"{lbl_m} ({len(s)})")
    if len(e7) >= 5:
        rho, pv = spearmanr(e7["luc_efficacy_pct"], e7["dist_H12c_H11axis_A"])
        ax.text(0.03, 0.94, f"cl7 Spearman ρ={rho:.2f}, p={pv:.2g}\n(n={len(e7)} pairs)", transform=ax.transAxes, va="top", fontsize=8.2)
    ax.axhspan(12.48 - 0.91, 12.48 + 0.91, color=cc.CL7, alpha=0.08)
    ax.set_xlabel("luc_efficacy_pct (% of E2)")
    ax.set_ylabel("dist_H12c_H11axis_A (Å, cl7)")
    ax.set_title("Geometry vs Luc efficacy inside cl7")
    ax.legend(loc="lower right", fontsize=7.5)
    cc.grid(ax, "both")

    # Bottom-right: dist_H12c_H11axis_A vs proliferation_pct inside cl7
    ax = axs[1, 2]
    for is_wtm, mk, ec, lbl_m in [(False, "o", "white", "matched context"), (True, "X", cc.INK, "WT assay / mutant struct")]:
        s = e7[e7["wt_assay_on_mutant_note"].astype(str).eq(str(is_wtm))]
        if len(s):
            ax.scatter(s["proliferation_pct"], s["dist_H12c_H11axis_A"],
                       c=[cc.CLASS_COLOR[c] for c in s["action_class"]], marker=mk, s=46, edgecolors=ec, lw=0.6, alpha=0.85, label=f"{lbl_m} ({len(s)})")
    if len(e7) >= 5:
        rho, pv = spearmanr(e7["proliferation_pct"], e7["dist_H12c_H11axis_A"])
        ax.text(0.03, 0.94, f"cl7 Spearman ρ={rho:.2f}, p={pv:.2g}\n(n={len(e7)} pairs)", transform=ax.transAxes, va="top", fontsize=8.2)
    ax.axhspan(12.48 - 0.91, 12.48 + 0.91, color=cc.CL7, alpha=0.08)
    ax.set_xlabel("proliferation_pct (% of E2)")
    ax.set_ylabel("dist_H12c_H11axis_A (Å, cl7)")
    ax.set_title("Geometry vs Proliferation inside cl7")
    cc.grid(ax, "both")

    cc.wrap_title(fig, f"Continuous efficacy data cover {len(eff_sub)} strict pairs (mostly 2 DOIs, derived_from_data); within cl7, H12 distance stays flat across efficacy.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    cc.save_fig_bundle(
        fig, "F12_structure_vs_efficacy", eff_pg[["pdb_norm", "lig_norm", "doi", "group", "action_class", "primary_source", "luc_efficacy_pct", "proliferation_pct", "dist_H12c_H11axis_A", "wt_assay_on_mutant_note"]],
        f"Top & bottom-left: luc_efficacy_pct vs proliferation_pct (% of E2, preserving negative values) faceted by DOI with dotted DATA_DICTIONARY.md rule thresholds (20%, 80%, 0%). Bottom-middle/right: H12-to-H11 distance vs continuous efficacy inside cl7, marking WT-assay/mutant-structure discrepancies (X). Sources: esr1_dataset_final.csv, geometry_per_chain.csv.",
        manifest, "Continuous functional efficacy (Luc & proliferation) vs H12 geometry",
        "esr1_dataset_final.csv; data/lid/geometry_per_chain.csv",
        "backup", "PDB × CCD pairs"
    )


def main():
    cc.ensure_dirs()
    chains = cc.read_csv_safe(cc.TABLES_DIR / "chains_curated.csv")
    pairs = cc.read_csv_safe(cc.TABLES_DIR / "pairs_curated.csv")
    pg = cc.read_csv_safe(cc.TABLES_DIR / "pairs_by_group.csv")
    eff_df, lodo_df, rob_df = compute_effects(pg, chains)
    manifest: list[dict] = []
    make_priority_figures(chains, pairs, pg, eff_df, lodo_df, manifest)
    pd.DataFrame(manifest).to_csv(cc.TABLES_DIR / "manifest_part1.csv", index=False)
    print(f"=== Stage 12.4a done: {len(manifest)} figures (F01-F12) ===")


if __name__ == "__main__":
    main()
