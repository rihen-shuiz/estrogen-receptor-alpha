"""Stage 12.4b - Appendix Figures F13-F20 and unified figure manifest."""
from __future__ import annotations

import subprocess
import sys
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402


def main():
    cc.ensure_dirs()
    cc.setup_mpl()
    chains = cc.read_csv_safe(cc.TABLES_DIR / "chains_curated.csv")
    pairs = cc.read_csv_safe(cc.TABLES_DIR / "pairs_curated.csv")
    pg = cc.read_csv_safe(cc.TABLES_DIR / "pairs_by_group.csv")
    rob_df = cc.read_csv_safe(cc.TABLES_DIR / "robustness_sensitivity.csv")

    core_pairs = pairs[pairs["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()
    core_pg = pg[pg["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()
    manifest: list[dict] = []
    rng = np.random.default_rng(cc.SEED)

    # ========================================================================== F13: Pocket & Filling
    pf = core_pg[(core_pg["pocket_complete_primary"].astype(str) == "True") &
                 core_pg["ligand_volume_A3"].notna() & core_pg["pocket_volume_A3"].notna()].copy()
    fig, axs = plt.subplots(1, 3, figsize=(13.6, 4.6))
    # Panel A: Ligand volume vs Pocket volume with iso-filling lines
    ax = axs[0]
    xv = np.linspace(150, 520, 100)
    for pct in (20, 30, 40, 50):
        ax.plot(xv, xv / (pct / 100.0), color=cc.GRID, lw=1.0, ls="--")
        ax.text(510, 510 / (pct / 100.0), f"{pct}%", fontsize=7.5, color=cc.MUTED, va="center")
    for cls in cc.CLASS_ORDER_CORE:
        for grp in ("cl7", "cl4"):
            s = pf[(pf["action_class"] == cls) & (pf["group"] == grp)]
            if len(s):
                ax.scatter(s["ligand_volume_A3"], s["pocket_volume_A3"], c=cc.CLASS_COLOR[cls],
                           marker=cc.GROUP_MARKER[grp], s=30, alpha=0.8, label=f"{cls}|{grp} ({len(s)})")
    ax.set_xlim(150, 560)
    ax.set_ylim(400, 1800)
    ax.set_xlabel("Ligand volume (Å³)")
    ax.set_ylabel("Pocket volume (Å³)")
    ax.set_title("A. Ligand vs pocket volume (iso-filling %)")
    ax.legend(loc="upper left", fontsize=7.0, ncol=2)
    cc.grid(ax, "both")

    # Panel B: Residual pocket volume (linear fit on ligand volume) by class
    ax = axs[1]
    coef = np.polyfit(pf["ligand_volume_A3"], pf["pocket_volume_A3"], 1)
    pf["pocket_resid_A3"] = pf["pocket_volume_A3"] - np.polyval(coef, pf["ligand_volume_A3"])
    for i, cls in enumerate(cc.CLASS_ORDER_CORE):
        v = pf[pf["action_class"] == cls]["pocket_resid_A3"].to_numpy()
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, c=cc.CLASS_COLOR[cls], s=26, alpha=0.75)
        ax.plot([i - 0.25, i + 0.25], [np.median(v), np.median(v)], color=cc.INK, lw=2.2)
    ax.axhline(0, color=cc.MUTED, ls="--", lw=1)
    ax.set_xticks(range(3), [c.replace("_", "\n") for c in cc.CLASS_ORDER_CORE])
    ax.set_ylabel("Pocket volume residual (Å³, size-adjusted)")
    ax.set_title("B. Pocket volume residual after ligand size")
    cc.grid(ax)

    # Panel C: filling_pct & pocket_pct_of_E2_pocket by class (equal ligand weight check)
    ax = axs[2]
    pf_lig = pf.groupby(["action_class", "inchikey"], as_index=False).agg(
        filling_pct=("filling_pct", "median"), pocket_pct_e2=("pocket_pct_of_E2_pocket", "median")
    )
    for i, cls in enumerate(cc.CLASS_ORDER_CORE):
        v = pf_lig[pf_lig["action_class"] == cls]["filling_pct"].to_numpy()
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, c=cc.CLASS_COLOR[cls], s=28, alpha=0.8)
        ax.plot([i - 0.25, i + 0.25], [np.median(v), np.median(v)], color=cc.INK, lw=2.2)
        ax.text(i, np.median(v) + 1.5, f"{np.median(v):.1f}%\n(n={len(v)})", ha="center", fontsize=7.8)
    ax.set_xticks(range(3), [c.replace("_", "\n") for c in cc.CLASS_ORDER_CORE])
    ax.set_ylabel("Pocket filling (%) — unique ligands")
    ax.set_title("C. Filling % (complete ligands, equal InChIKey)")
    cc.grid(ax)

    cc.wrap_title(fig, "Pocket volume scales tightly with ligand volume; filling % stays near ~33–36% across all three functional classes.", width=104)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F13_pocket_filling", pf[["pdb_norm", "lig_norm", "group", "action_class", "ligand_volume_A3", "pocket_volume_A3", "filling_pct", "pocket_pct_of_E2_pocket", "pocket_resid_A3"]],
        f"Complete-ligand primary pocket set (n={len(pf)} pair×group rows; n={len(pf_lig)} unique InChIKeys in C). Sources: data/ligands/pocket_filling_per_chain.csv, esr1_dataset_final.csv.",
        manifest, "Pocket volume, ligand volume, and pocket filling across classes",
        "data/ligands/pocket_filling_per_chain.csv; esr1_dataset_final.csv",
        "appendix", "PDB×CCD×group (A–B), unique InChIKey (C)"
    )

    # ========================================================================== F14: Chemistry descriptors & Murcko scaffolds
    chem = cc.read_csv_safe(cc.TABLES_DIR / "chem_ligands_order.csv")
    ld = cc.read_csv_safe(cc.DATA / "ligands/ligand_descriptors.csv").rename(columns={"ligand_id": "lig_norm"})
    cl = chem.merge(ld[["lig_norm", "clogp", "tpsa", "hbd", "hba", "rotb", "aromatic_rings"]], on="lig_norm", how="left")

    fig, axs = plt.subplots(2, 2, figsize=(12.4, 8.0))
    for ax, col, ylab, ttl in [
        (axs[0, 0], "heavy_atoms", "Heavy atoms (count)", "A. Ligand size (heavy atoms)"),
        (axs[0, 1], "clogp", "cLogP", "B. Lipophilicity (cLogP)"),
        (axs[1, 0], "tpsa", "TPSA (Å²)", "C. Polar surface area (TPSA)"),
    ]:
        for i, cls in enumerate(cc.CLASS_ORDER_CORE):
            v = cl[cl["action_class_majority"] == cls][col].dropna().to_numpy()
            ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, c=cc.CLASS_COLOR[cls], s=26, alpha=0.75)
            ax.plot([i - 0.25, i + 0.25], [np.median(v), np.median(v)], color=cc.INK, lw=2.2)
            ax.text(i, np.max(v) + (np.max(v) - np.min(v)) * 0.04, f"med={np.median(v):.1f}", ha="center", fontsize=8)
        ax.set_xticks(range(3), [c.replace("_", "\n") for c in cc.CLASS_ORDER_CORE])
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        cc.grid(ax)

    # Panel D: Top 8 Murcko scaffolds by class
    ax = axs[1, 1]
    top_scaf = core_pairs["murcko"].replace("", np.nan).dropna().value_counts().head(8).index.tolist()
    sc_tab = pd.crosstab(core_pairs["murcko"], core_pairs["action_class"]).reindex(index=top_scaf, columns=cc.CLASS_ORDER_CORE, fill_value=0)
    sc_dois = core_pairs.groupby("murcko")["doi_block"].nunique().reindex(top_scaf).fillna(0).astype(int)
    y = np.arange(len(top_scaf))
    left = np.zeros(len(top_scaf))
    for cls in cc.CLASS_ORDER_CORE:
        v = sc_tab[cls].values
        ax.barh(y, v, left=left, color=cc.CLASS_COLOR[cls], height=0.6, label=cls)
        left += v
    for i, (tot, nd) in enumerate(zip(left, sc_dois.values)):
        ax.text(tot + 0.6, i, f"n={int(tot)} ({nd} DOI)", va="center", fontsize=7.6)
    ax.set_yticks(y, [f"S{i+1}: {s[:20]}…" if len(s) > 20 else f"S{i+1}: {s}" for i, s in enumerate(top_scaf)], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, left.max() * 1.32)
    ax.set_xlabel("Core-3 strict pairs")
    ax.set_title("D. Top-8 Murcko scaffolds by class & DOI count")
    ax.legend(loc="lower right", fontsize=7.5)
    cc.grid(ax, "x")

    cc.wrap_title(fig, "Antagonists are substantially larger and more lipophilic than agonists/partial agonists; top scaffolds are concentrated in specific DOI series.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    cc.save_fig_bundle(
        fig, "F14_chemistry_descriptors_scaffolds", cl[["lig_norm", "action_class_majority", "heavy_atoms", "clogp", "tpsa", "basic_amine", "murcko"]],
        f"A–C: Molecular descriptors aggregated at unique CCD ligand level (n={len(cl)}). D: Top-8 Murcko scaffolds across 291 strict core-3 pairs with DOI counts. Sources: data/ligands/ligand_descriptors.csv, esr1_dataset_final.csv.",
        manifest, "Chemical descriptors and Murcko scaffolds by curated class",
        "data/ligands/ligand_descriptors.csv; esr1_dataset_final.csv",
        "appendix", "unique CCD ligands (A–C), PDB×CCD pairs (D)"
    )

    # ========================================================================== F15: Pharmacology (drug_classes x action_class)
    cur_real = cc.read_csv_safe(cc.ROOT / "esr1_dataset_final.csv")
    cur_real = cur_real[cur_real["is_real_ligand"].astype(str) == "True"].copy()
    filled = cur_real[cur_real["drug_classes"].notna()].copy()
    # Expand multi-label drug_classes (separated by '|')
    dc_rows = []
    for _, r in filled.iterrows():
        for dc in str(r["drug_classes"]).split("|"):
            dc = dc.strip()
            if dc:
                dc_rows.append({"pdb_id": r["pdb_id"], "ligand_id": r["ligand_id"], "drug_class": dc, "action_class": r["action_class"]})
    dc_df = pd.DataFrame(dc_rows)
    dc_ord = ["SERD", "SERM", "SERCA", "CERAN", "DMERI"]
    ac_ord = ["antagonist", "agonist", "partial_agonist", "mixed", "unknown"]
    mat_dc = pd.crosstab(dc_df["drug_class"], dc_df["action_class"]).reindex(index=dc_ord, columns=ac_ord, fill_value=0)

    fig, axs = plt.subplots(1, 2, figsize=(11.6, 4.4), gridspec_kw={"width_ratios": [0.7, 1.4]})
    ax = axs[0]
    cov = [len(filled), len(cur_real) - len(filled)]
    ax.bar(["Annotated\ndrug_classes", "Empty\ndrug_classes"], cov, color=["#4a3aa7", "#d9d8d2"], width=0.52)
    for i, v in enumerate(cov):
        ax.text(i, v + 6, f"{v}\n({v/len(cur_real):.1%})", ha="center", fontsize=9, fontweight="bold")
    ax.set_ylim(0, 360)
    ax.set_ylabel("Real ligand pairs (total n=368)")
    ax.set_title("A. Pharmacological axis coverage")
    cc.grid(ax)

    ax = axs[1]
    ax.imshow(mat_dc.values, cmap="Purples", aspect="auto")
    for i in range(mat_dc.shape[0]):
        for j in range(mat_dc.shape[1]):
            v = int(mat_dc.values[i, j])
            ax.text(j, i, str(v) if v > 0 else "·", ha="center", va="center", fontsize=9.5,
                    color="white" if v > mat_dc.values.max() * 0.5 else cc.INK, fontweight="bold" if v > 0 else "normal")
    ax.set_xticks(range(len(ac_ord)), [c.replace("_", "\n") for c in ac_ord], fontsize=8.5)
    ax.set_yticks(range(len(dc_ord)), [f"{d} (n={int(mat_dc.loc[d].sum())})" for d in dc_ord], fontsize=8.8)
    ax.set_title(f"B. drug_classes × action_class (54 annotated pairs; note 9 unknown SERDs)")

    cc.wrap_title(fig, "Pharmacological drug_classes are populated for only 54/368 real pairs (14.7%), including 9 SERDs with unknown functional action.", width=102)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    cc.save_fig_bundle(
        fig, "F15_pharmacology_upset", mat_dc.reset_index(),
        "A: Coverage of drug_classes across 368 real ligand pairs in esr1_dataset_final.csv. B: Intersection counts of pipe-separated drug_classes vs action_class (n=54 annotated pairs, including 9 SERDs with action_class='unknown' where degradation was tested without separate antagonism assay). Sources: esr1_dataset_final.csv.",
        manifest, "Pharmacological drug_classes coverage and intersection with action_class",
        "esr1_dataset_final.csv",
        "appendix", "PDB × CCD real ligand pairs (n=368)"
    )

    # ========================================================================== F16: Non-H12 Footprint (Stage 07 geometry + Sammon 4-panel)
    gf = cc.read_csv_safe(cc.RESULTS_BASE / "07_geometry_free/geometry_free_per_chain.csv")
    maps = cc.read_csv_safe(cc.RESULTS_BASE / "05_maps/maps_coords.csv")
    ch_m = chains.merge(maps[["chain_key", "sammon1", "sammon2"]], on="chain_key", how="left")
    ch_m = ch_m.merge(gf.drop(columns=["entry_id", "group"], errors="ignore"), on="chain_key", how="left")

    # Compute inter-helix distance effect (antagonist - agonist) inside cl7 Y537S+pep vs overall cl4-cl7
    dist_cols = [c for c in gf.columns if c.startswith("d_H") and "H12" not in c]
    nonh12_rows = []
    sub7 = ch_m[(ch_m["in_core3_strict"].astype(str) == "True") & (ch_m["group"] == "cl7") &
                (ch_m["construct_mutations"] == "Y537S") & (ch_m["peptide_bound"].astype(str) == "True")]
    sub7_p = sub7.groupby(["pdb_norm", "lig_norm"], as_index=False).agg(
        action_class=("action_class", "first"), **{c: (c, "median") for c in dist_cols}
    )
    for col in dist_cols:
        a = sub7_p[sub7_p["action_class"] == "antagonist"][col].dropna()
        b = sub7_p[sub7_p["action_class"] == "agonist"][col].dropna()
        if len(a) >= 3 and len(b) >= 3:
            nonh12_rows.append({"pair": col.replace("d_", ""), "diff_A": float(a.median() - b.median()),
                                "p": float(mannwhitneyu(a, b).pvalue)})
    nh_df = pd.DataFrame(nonh12_rows)
    if len(nh_df):
        nh_df["q"] = cc.bh_fdr(nh_df["p"])
        nh_df = nh_df.sort_values("diff_A")

    fig, axs = plt.subplots(2, 2, figsize=(12.2, 8.2))
    # Panel A: Non-H12 inter-helix shifts (ant vs ago inside cl7 Y537S+pep)
    ax = axs[0, 0]
    if len(nh_df):
        top_nh = pd.concat([nh_df.head(6), nh_df.tail(6)]).drop_duplicates()
        y = np.arange(len(top_nh))
        ax.barh(y, top_nh["diff_A"], color=["#c2185b" if v > 0 else "#1b9e77" for v in top_nh["diff_A"]], height=0.6)
        ax.axvline(0, color=cc.INK, lw=0.9)
        ax.set_yticks(y, top_nh["pair"], fontsize=8)
        ax.set_xlabel("Median shift ant − ago in cl7 Y537S+pep (Å)")
        ax.set_title("A. Non-H12 inter-helix shifts (all ≤ 0.3 Å)")
        cc.grid(ax, "x")

    # Panels B, C, D: Same Sammon coordinates colored by curated class, H12 group, and construct
    for ax, col_by, ttl in [
        (axs[0, 1], "action_class", "B. Sammon map colored by curated class"),
        (axs[1, 0], "group", "C. Sammon map colored by H12 group"),
        (axs[1, 1], "mutation_status", "D. Sammon map colored by construct"),
    ]:
        if col_by == "action_class":
            for k in cc.CLASS_ORDER_ALL:
                s = ch_m[ch_m["action_class"] == k]
                ax.scatter(s["sammon1"], s["sammon2"], c=cc.CLASS_COLOR[k], s=18, alpha=0.8, label=f"{k} ({len(s)})")
        elif col_by == "group":
            for k in ("cl7", "cl4", "minor"):
                s = ch_m[ch_m["group"] == k]
                ax.scatter(s["sammon1"], s["sammon2"], c=cc.GROUP_COLOR[k], s=18, alpha=0.8, label=f"{k} ({len(s)})")
        else:
            mut_cols = {"Y537S": "#2a78d6", "L536S": "#eb6834", "wild type": "#1b9e77", "Cys->Ser/Ala only": "#7b3294", "other": "#8a8984"}
            for k, col in mut_cols.items():
                s = ch_m[ch_m["mutation_status"] == k]
                if len(s):
                    ax.scatter(s["sammon1"], s["sammon2"], c=col, s=18, alpha=0.8, label=f"{k} ({len(s)})")
        ax.set_xlabel("Sammon 1")
        ax.set_ylabel("Sammon 2")
        ax.set_title(ttl)
        ax.legend(loc="best", fontsize=7.3)
        cc.grid(ax, "both")

    cc.wrap_title(fig, "Outside H12, inter-helix distances shift by < 0.3 Å between classes; the Sammon map is dominated by the H12/construct axis.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    cc.save_fig_bundle(
        fig, "F16_non_h12_footprint", nh_df,
        "A: Inter-helix COM distance differences (excluding H12, Stage 07) between antagonist and agonist inside cl7 Y537S+peptide. B–D: Same Stage-05 Sammon coordinates (788 chains) colored by curated action_class, H12 group, and construct mutation_status. Sources: results/07_geometry_free/geometry_free_per_chain.csv, results/05_maps/maps_coords.csv.",
        manifest, "Structural footprint of functional class outside H12",
        "results/07_geometry_free/geometry_free_per_chain.csv; results/05_maps/maps_coords.csv; esr1_dataset_final.csv",
        "appendix", "PDB×CCD pairs (A), chains (B–D)"
    )

    # ========================================================================== F17: cl4 Heterogeneity (1.0 A subgroups)
    c4 = core_pg[(core_pg["group"] == "cl4") & core_pg["cl4_sub_cluster"].notna()].copy()
    fig, axs = plt.subplots(1, 3, figsize=(13.4, 4.5))
    ax = axs[0]
    ct4 = pd.crosstab(c4["cl4_sub_cluster"].astype(int), c4["action_class"]).reindex(columns=cc.CLASS_ORDER_CORE, fill_value=0)
    y = np.arange(len(ct4))
    left = np.zeros(len(ct4))
    for cls in cc.CLASS_ORDER_CORE:
        v = ct4[cls].values
        ax.barh(y, v, left=left, color=cc.CLASS_COLOR[cls], height=0.58, label=cls)
        left += v
    ax.set_yticks(y, [f"cl4 sub-{k} (n={int(ct4.loc[k].sum())})" for k in ct4.index], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("Core-3 pairs in cl4")
    ax.set_title("A. cl4 1.0 Å subgroups × curated class")
    ax.legend(loc="lower right", fontsize=8)
    cc.grid(ax, "x")

    for ax, col, ylab, ttl in [
        (axs[1], "h12_shift_A", "h12_shift_A (Å)", "B. H12 shift across cl4 subgroups"),
        (axs[2], "k362_e542_A", "K362–E542 Cα distance (Å)", "C. K362–E542 salt-bridge distance"),
    ]:
        subs = sorted(c4["cl4_sub_cluster"].dropna().astype(int).unique())[:4]
        for i, sc_id in enumerate(subs):
            s = c4[c4["cl4_sub_cluster"].astype(int) == sc_id].dropna(subset=[col])
            for cls in cc.CLASS_ORDER_CORE:
                v = s[s["action_class"] == cls][col].to_numpy()
                if len(v):
                    ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, c=cc.CLASS_COLOR[cls], s=28, alpha=0.8)
            if len(s):
                ax.plot([i - 0.24, i + 0.24], [s[col].median(), s[col].median()], color=cc.INK, lw=2.2)
        ax.set_xticks(range(len(subs)), [f"sub-{k}" for k in subs])
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        cc.grid(ax)

    cc.wrap_title(fig, "Inside cl4, 1.0 Å subgroups reflect L536S construct geometry rather than separating antagonists from partial agonists.", width=104)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    cc.save_fig_bundle(
        fig, "F17_cl4_heterogeneity", c4[["pdb_norm", "lig_norm", "action_class", "mutation_status", "cl4_sub_cluster", "h12_shift_A", "k362_e542_A"]],
        f"cl4 1.0 Å subgroups (data/lid/cl4_subgroups.csv) across n={len(c4)} core-3 cl4 pairs. Sources: data/lid/cl4_subgroups.csv, esr1_dataset_final.csv.",
        manifest, "Heterogeneity inside cl4 (1.0 Å subgroups) vs functional class",
        "data/lid/cl4_subgroups.csv; esr1_dataset_final.csv",
        "appendix", "PDB × CCD pairs in cl4"
    )

    # ========================================================================== F18: Affinity Coverage & dpx_median
    aff = cc.read_csv_safe(cc.DATA / "ligands/affinity_per_ligand.csv").rename(columns={"ccd": "lig_norm"})
    chem = cc.read_csv_safe(cc.TABLES_DIR / "chem_ligands_order.csv")
    af_m = chem.merge(aff[["lig_norm", "affinity_basis", "dpx_median", "px_abs", "px_abs_type", "filling_pct_median"]], on="lig_norm", how="left")
    af_m["dpx_median"] = pd.to_numeric(af_m["dpx_median"], errors="coerce")
    af_m["filling_pct_median"] = pd.to_numeric(af_m["filling_pct_median"], errors="coerce")
    # Also attach median dist_H12c_H11axis_A per ligand
    dist_lig = core_pg.groupby("lig_norm")["dist_H12c_H11axis_A"].median().reset_index()
    af_m = af_m.merge(dist_lig, on="lig_norm", how="left")

    fig, axs = plt.subplots(1, 2, figsize=(11.8, 4.6))
    ax = axs[0]
    bc = af_m["affinity_basis"].fillna("no_chembl_affinity").value_counts()
    ax.barh(range(len(bc)), bc.values, color=["#8a8984", "#2a78d6", "#1b9e77", "#eb6834"][:len(bc)], height=0.55)
    for i, v in enumerate(bc.values):
        ax.text(v + 2, i, f"{v} ({v/len(af_m):.0%})", va="center", fontsize=8.5)
    ax.set_yticks(range(len(bc)), bc.index, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, bc.max() * 1.25)
    ax.set_xlabel("Unique core-3 CCD ligands")
    ax.set_title("A. ChEMBL affinity coverage (relative vs absolute)")
    cc.grid(ax, "x")

    ax = axs[1]
    rel = af_m[af_m["affinity_basis"] == "relative"].dropna(subset=["dpx_median", "dist_H12c_H11axis_A"])
    for cls in cc.CLASS_ORDER_CORE:
        s = rel[rel["action_class_majority"] == cls]
        if len(s):
            ax.scatter(s["dpx_median"], s["dist_H12c_H11axis_A"], c=cc.CLASS_COLOR[cls],
                       marker=cc.CLASS_MARKER[cls], s=60, edgecolors=cc.INK, lw=0.6, label=f"{cls} ({len(s)})")
            for _, r in s.iterrows():
                ax.annotate(r["lig_norm"], (r["dpx_median"], r["dist_H12c_H11axis_A"]),
                            xytext=(4, 4), textcoords="offset points", fontsize=7.5)
    ax.axvline(0, color=cc.MUTED, ls="--", lw=1)
    ax.set_xlabel("dpx_median (pX − pX_E2, same-assay relative affinity)")
    ax.set_ylabel("dist_H12c_H11axis_A (Å)")
    ax.set_title(f"B. Relative affinity vs H12 distance (descriptive, n={len(rel)})")
    ax.legend(loc="best", fontsize=8)
    cc.grid(ax, "both")

    cc.wrap_title(fig, f"Same-assay relative affinity (dpx_median) covers only {len(rel)} core-3 ligands, keeping affinity–structure analysis descriptive.", width=102)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    cc.save_fig_bundle(
        fig, "F18_affinity_coverage", af_m[["lig_norm", "action_class_majority", "h12_group_majority", "affinity_basis", "dpx_median", "px_abs", "dist_H12c_H11axis_A"]],
        f"A: Coverage of ChEMBL binding affinity across {len(af_m)} unique core-3 ligands. B: Same-assay E2-normalized relative affinity (dpx_median) vs H12-to-H11 distance (n={len(rel)}). Sources: data/ligands/affinity_per_ligand.csv, esr1_dataset_final.csv.",
        manifest, "Affinity coverage and relative binding affinity vs H12 geometry",
        "data/ligands/affinity_per_ligand.csv; esr1_dataset_final.csv",
        "appendix", "unique CCD ligands"
    )

    # ========================================================================== F19: ESR2 Context (Stage 11 + sv1028-run displacement)
    bm = cc.read_csv_safe(cc.RESULTS_BASE / "11_esr2/joint_block_medians.csv")
    fig, axs = plt.subplots(1, 2, figsize=(12.2, 4.6))
    ax = axs[0]
    y = np.arange(len(bm))
    cols_bm = ["#2a78d6", "#eb6834", "#6da7ec", "#f29c76", "#4a3aa7"]
    ax.barh(y, bm["median_rmsd_A"], color=cols_bm[:len(bm)], height=0.58)
    for i, r in enumerate(bm.itertuples()):
        ax.text(r.median_rmsd_A + 0.08, i, f"{r.median_rmsd_A:.2f} Å", va="center", fontsize=8.8, fontweight="bold")
    ax.set_yticks(y, [f"{r.group_a} ↔ {r.group_b}" for r in bm.itertuples()], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 4.3)
    ax.set_xlabel("Median joint Cα RMSD (Å, Stage 11)")
    ax.set_title("A. ESR1 (788 ch) vs ESR2 (86 ch) group distances")
    cc.grid(ax, "x")

    ax = axs[1]
    try:
        disp_txt = subprocess.check_output(["git", "show", "origin/sv1028-run:results/displacement_matrix.csv"], cwd=cc.ROOT, text=True)
        disp = pd.read_csv(StringIO(disp_txt), index_col=0)
        pos_idx = disp.index.astype(int)
        med_disp = disp.median(axis=1)
        q25 = disp.quantile(0.25, axis=1)
        q75 = disp.quantile(0.75, axis=1)
        ax.plot(pos_idx, med_disp, color="#4a3aa7", lw=1.8, label=f"Median ({disp.shape[1]} ESR1–ESR2 pairs)")
        ax.fill_between(pos_idx, q25, q75, color="#4a3aa7", alpha=0.18, label="IQR")
        ax.axvspan(528, 550, color=cc.CL4, alpha=0.12, label="H12 region (528–550)")
        ax.set_xlabel("ESR1 UniProt residue position")
        ax.set_ylabel("Equivalent Cα displacement (Å)")
        ax.set_title("B. Shared-ligand ESR1 vs ESR2 Cα displacement")
        ax.legend(loc="upper left", fontsize=7.8)
        cc.grid(ax)
    except Exception as e:
        ax.text(0.5, 0.5, f"displacement_matrix unavailable:\n{e}", ha="center", va="center", transform=ax.transAxes)

    cc.wrap_title(fig, "ESR2 (86 chains, 39 PDBs) mirrors the same two-state H12 split: ESR2 main is 1.33 Å from ESR1 cl7 vs 3.56 Å from cl4.", width=104)
    fig.tight_layout(rect=[0, 0, 1, 0.91])
    cc.save_fig_bundle(
        fig, "F19_esr2_context", bm,
        "A: Stage-11 joint Cα RMSD block medians between ESR1 (788 chains, 351 PDB) and ESR2 (86 chains, 39 PDB). B: Residue-wise Cα displacement across shared-ligand ESR1–ESR2 pairs from sv1028-run (results/displacement_matrix.csv). Sources: results/11_esr2/joint_block_medians.csv, origin/sv1028-run:results/displacement_matrix.csv.",
        manifest, "ESR2 structural context and shared-ligand displacement",
        "results/11_esr2/joint_block_medians.csv; origin/sv1028-run:results/displacement_matrix.csv",
        "backup", "chains / shared-ligand pairs"
    )

    # ========================================================================== F20: Robustness Sensitivity Forest Plot
    fig, axs = plt.subplots(1, 2, figsize=(13.2, 5.4), sharey=True)
    for ax, cname, ttl in [
        (axs[0], "ant_vs_ago", "A. cl7 Y537S+pep: antagonist − agonist"),
        (axs[1], "part_vs_ago", "B. cl7 Y537S+pep: partial_agonist − agonist"),
    ]:
        d = rob_df[rob_df["contrast"] == cname].reset_index(drop=True)
        y = np.arange(len(d))
        ax.axvspan(-0.91, 0.91, color=cc.CL7, alpha=0.09, label="cl7 ordinary ±0.91 Å")
        ax.axvline(0, color=cc.INK2, lw=0.9)
        col = cc.CLASS_COLOR["antagonist"] if cname == "ant_vs_ago" else cc.CLASS_COLOR["partial_agonist"]
        for i, r in enumerate(d.itertuples()):
            if pd.notna(r.diff_A):
                ax.scatter(r.diff_A, i, color=col, s=46, zorder=4)
                if str(r.ci_valid) == "True" and pd.notna(r.ci95_lo):
                    ax.plot([r.ci95_lo, r.ci95_hi], [i, i], color=col, lw=1.9, zorder=3)
                txt = f"{r.diff_A:+.2f} Å ({int(r.n_a)}v{int(r.n_b)}, {int(r.dois_a)}/{int(r.dois_b)} DOI)"
            else:
                txt = f"n/a ({int(r.n_a)}v{int(r.n_b)})"
            ax.text(0.95, i, txt, va="center", ha="left", fontsize=7.5, color=cc.INK)
        ax.set_yticks(y, d["filter"], fontsize=8.2)
        ax.invert_yaxis()
        ax.set_xlim(-1.05, 1.95)
        ax.set_xlabel("Median Δ dist_H12c_H11axis_A (Å, 95% DOI-block CI)")
        ax.set_title(ttl)
        cc.grid(ax, "x")
    axs[0].legend(loc="lower left", fontsize=8)

    cc.wrap_title(fig, "Robustness check: across all 11 quality, completeness, evidence, and deduplication filters, within-cl7 shifts stay well below ±0.91 Å.", width=106)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    cc.save_fig_bundle(
        fig, "F20_robustness_sensitivity", rob_df,
        "Median shifts (with 95% DOI-block bootstrap CIs where ≥4 DOI blocks exist per class) for the main comparable stratum (cl7 Y537S + peptide) across 11 sensitivity filters. Sources: robustness_sensitivity.csv.",
        manifest, "Robustness of within-group H12 effect across 11 analytical filters",
        "results/12_curated_action/tables/robustness_sensitivity.csv",
        "appendix", "PDB × CCD pairs (or unique InChIKeys in filter 11)"
    )

    # Combine manifest part1 + part2
    m1 = cc.read_csv_safe(cc.TABLES_DIR / "manifest_part1.csv")
    m2 = pd.DataFrame(manifest)
    full_m = pd.concat([m1, m2], ignore_index=True)
    full_m.to_csv(cc.TABLES_DIR / "figure_manifest.csv", index=False)
    print(f"=== Stage 12.4b done: {len(m2)} appendix figures (F13-F20); total manifest = {len(full_m)} figures ===")


if __name__ == "__main__":
    main()
