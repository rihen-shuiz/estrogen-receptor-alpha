"""Stage 12.1 - Exact dataset join, coverage audit, construct/peptide reconciliation, and master tables.

Usage:
  python pipeline/curated_action/01_join_and_audit.py --dataset esr1_dataset_final.csv
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402


def get_git_sha() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=cc.ROOT, text=True).strip()
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default=str(cc.ROOT / "esr1_dataset_final.csv"))
    args = ap.parse_args()
    cc.ensure_dirs()

    dataset_path = Path(args.dataset).resolve()
    cur = cc.read_csv_safe(dataset_path)
    cur["pdb_norm"] = cur["pdb_id"].astype(str).str.strip().str.upper()
    cur["lig_norm"] = cur["ligand_id"].astype(str).str.strip().str.upper()
    cur["is_real_bool"] = cur["is_real_ligand"].astype(str).str.strip().eq("True")
    cur["is_mutant_bool"] = cur["is_mutant"].astype(str).str.strip().eq("True")
    cur["has_peptide_bool"] = cur["has_peptide"].astype(str).str.strip().eq("True")
    for num_col in ("luc_efficacy_pct", "proliferation_pct"):
        cur[num_col] = pd.to_numeric(cur[num_col], errors="coerce")

    assert cur[["pdb_norm", "lig_norm"]].duplicated().sum() == 0, "Curated CSV must be unique on (pdb_id, ligand_id)"

    cur_real = cur[cur["is_real_bool"]].copy()

    # Load base chain_table (788 chains, 351 PDB)
    ct = cc.read_csv_safe(cc.DATA / "structures/chain_table.csv")
    ct["pdb_norm"] = ct["entry_id"].astype(str).str.strip().str.upper()
    ct["lig_norm"] = ct["ligand_id_primary"].fillna("").astype(str).str.strip().str.upper()
    assert len(ct) == 788 and ct["pdb_norm"].nunique() == 351

    # Verify groups against results/02_lid_groups/groups.csv
    gr = cc.read_csv_safe(cc.RESULTS_BASE / "02_lid_groups/groups.csv")
    ct = ct.merge(gr[["chain_key", "cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5"]], on="chain_key", how="left", validate="one_to_one")
    assert (ct["cut_1.5"].astype(int) == ct["cluster_1p5"].astype(int)).all(), "cluster_1p5 != cut_1.5"
    ct["group"] = ct["cut_1.5"].map(cc.group_of)

    # Historical ChEMBL action_class (labelled explicitly as old_chembl_class)
    old_cc = cc.read_csv_safe(cc.DATA / "ligands/chain_class.csv")[["chain_key", "action_class"]].rename(
        columns={"action_class": "old_chembl_class"}
    )
    ct = ct.merge(old_cc, on="chain_key", how="left", validate="one_to_one")
    ct["old_chembl_class"] = ct["old_chembl_class"].fillna("unlabelled")

    # Entries metadata
    ent = cc.read_csv_safe(cc.DATA / "structures/entries.csv")[
        ["entry_id", "deposit_date", "qc_status", "qc_res", "qc_rfree", "qc_b", "qc_pass"]
    ].rename(columns={"entry_id": "pdb_norm"})
    ct = ct.merge(ent, on="pdb_norm", how="left", validate="many_to_one")
    assert len(ct) == 788

    # Geometry, mobility, lattice, outliers, pocket filling, cl4 subgroups
    geom = cc.read_csv_safe(cc.DATA / "lid/geometry_per_chain.csv")
    geom_cols = [
        "chain_key", "dist_H12c_H11axis_A", "dist_H12c_H3axis_A",
        "angle_H12_H3_deg", "angle_H12_H11_deg", "angle_H12_H45_deg",
        "hinge_H11C_H12_deg", "kink_H11_deg", "len_H12", "len_H12_w544",
        "dist_H12c_H11axis_A_w544", "angle_H12_H3_deg_w544", "angle_H12_H11_deg_w544",
    ]
    ct = ct.merge(geom[geom_cols], on="chain_key", how="left", validate="one_to_one")

    mob = cc.read_csv_safe(cc.DATA / "lid/mobility_per_chain.csv")[
        ["chain_key", "h12_status", "mobility_h12_z", "mobility_z", "n_obs_538_547", "n_obs_528_537", "core_mean_b"]
    ]
    ct = ct.merge(mob, on="chain_key", how="left", validate="one_to_one")

    # Also compute H12 538-544 completeness from aligned_ca.npz
    z_ca = np.load(cc.RESULTS_BASE / "02_lid_groups/aligned_ca.npz", allow_pickle=True)
    ca_keys = [str(k) for k in z_ca["chain_keys"]]
    ca_pos = z_ca["positions"].astype(int)
    ca_xyz = z_ca["coords"].astype(float)
    m544 = (ca_pos >= 538) & (ca_pos <= 544)
    m547 = (ca_pos >= 538) & (ca_pos <= 547)
    n544 = np.isfinite(ca_xyz[:, m544, 0]).sum(axis=1)
    n547 = np.isfinite(ca_xyz[:, m547, 0]).sum(axis=1)
    h12_win = pd.DataFrame({
        "chain_key": ca_keys,
        "n_obs_538_544": n544,
        "h12_complete_544": n544 == 7,
        "h12_complete_547": n547 == 10,
    })
    ct = ct.merge(h12_win, on="chain_key", how="left", validate="one_to_one")

    lat = cc.read_csv_safe(cc.DATA / "lid/lattice_contacts_per_chain.csv")[
        ["chain_key", "sym_atoms", "asu_atoms", "sym_receptor_atoms", "n_residues_sym_contact"]
    ]
    ct = ct.merge(lat, on="chain_key", how="left", validate="one_to_one")

    of = cc.read_csv_safe(cc.DATA / "outliers/outlier_flags.csv")[["chain_key", "n_methods", "consensus_outlier"]]
    ct = ct.merge(of, on="chain_key", how="left", validate="one_to_one")
    ov = cc.read_csv_safe(cc.DATA / "outliers/outlier_verdicts.csv")[["chain_key", "verdict", "reasons"]].rename(
        columns={"verdict": "outlier_verdict", "reasons": "outlier_reasons"}
    )
    ct = ct.merge(ov, on="chain_key", how="left", validate="one_to_one")

    pf = cc.read_csv_safe(cc.DATA / "ligands/pocket_filling_per_chain.csv")
    pf_cols = [
        "chain_key", "ligand_ids_used", "ligand_incomplete", "in_primary_analysis",
        "pocket_volume_A3", "ligand_volume_A3", "filling_pct",
        "ligand_pct_of_E2", "pocket_pct_of_E2_pocket", "pocket_open",
    ]
    pf_ren = pf[pf_cols].rename(columns={
        "in_primary_analysis": "pocket_in_primary",
        "ligand_incomplete": "pocket_ligand_incomplete",
    })
    ct = ct.merge(pf_ren, on="chain_key", how="left", validate="one_to_one")

    cl4s = cc.read_csv_safe(cc.DATA / "lid/cl4_subgroups.csv")[
        ["chain_key", "sub_cluster", "h12_shift_A", "h12_h3_angle_deg", "k362_e542_A", "h12_contact_atoms"]
    ].rename(columns={
        "sub_cluster": "cl4_sub_cluster",
        "h12_h3_angle_deg": "cl4_h12_h3_angle_deg",
        "h12_contact_atoms": "cl4_h12_contact_atoms",
    })
    ct = ct.merge(cl4s, on="chain_key", how="left", validate="one_to_one")

    # Ligand descriptors and IDs by primary CCD
    ld = cc.read_csv_safe(cc.DATA / "ligands/ligand_descriptors.csv").rename(columns={
        "ligand_id": "lig_norm",
        "n_chains": "ld_n_chains",
        "basic_amine": "ld_basic_amine",
        "heavy_atoms": "ld_heavy_atoms",
    })
    ct = ct.merge(ld, on="lig_norm", how="left", validate="many_to_one")

    lids = cc.read_csv_safe(cc.DATA / "ligands/ligand_ids.csv")[
        ["ccd", "smiles_stereo", "inchikey", "parent_pref_name", "chembl_parent", "is_e2_reference"]
    ].drop_duplicates("ccd").rename(columns={
        "ccd": "lig_norm",
        "inchikey": "repo_inchikey",
        "parent_pref_name": "repo_ligand_name",
    })
    ct = ct.merge(lids, on="lig_norm", how="left", validate="many_to_one")

    aff = cc.read_csv_safe(cc.DATA / "ligands/affinity_per_ligand.csv")[
        ["ccd", "affinity_basis", "dpx_median", "rba_median", "n_paired_assays", "px_abs", "px_abs_type", "dpx_plot"]
    ].drop_duplicates("ccd").rename(columns={"ccd": "lig_norm"})
    ct = ct.merge(aff, on="lig_norm", how="left", validate="many_to_one")

    assert len(ct) == 788, f"Row count changed after merges: {len(ct)}"

    # Exact primary-CCD join with curated dataset (real ligands)
    cur_join_cols = [
        "pdb_norm", "lig_norm", "doi", "action_class", "drug_classes", "compound_id",
        "mutations", "mutations_norm", "mutation_group", "is_mutant_bool", "is_real_bool",
        "has_peptide_bool", "inchikey", "evidence_sources", "primary_source",
        "luc_efficacy_pct", "proliferation_pct", "pdb_title", "note",
    ]
    merged = ct.merge(cur_real[cur_join_cols], on=["pdb_norm", "lig_norm"], how="left", validate="many_to_one")
    assert len(merged) == 788

    merged["curated_matched_primary"] = merged["action_class"].notna()
    merged["in_strict_single"] = merged["curated_matched_primary"] & (merged["ligand_status"] == "single")
    merged["in_core3_strict"] = merged["in_strict_single"] & merged["action_class"].isin(cc.CLASS_ORDER_CORE)
    merged["in_core3_with_multi"] = merged["curated_matched_primary"] & merged["action_class"].isin(cc.CLASS_ORDER_CORE)
    merged["doi_block"] = cc.doi_block_id(merged)

    # Numeric coercions
    for col in [
        "dist_H12c_H11axis_A", "dist_H12c_H3axis_A", "angle_H12_H3_deg", "angle_H12_H11_deg",
        "len_H12", "mobility_h12_z", "resolution_A", "r_free", "mean_b_entry_A2",
        "sym_atoms", "asu_atoms", "pocket_volume_A3", "ligand_volume_A3", "filling_pct",
        "pocket_pct_of_E2_pocket", "ligand_pct_of_E2", "ld_heavy_atoms", "mw", "clogp", "tpsa",
        "h12_shift_A", "k362_e542_A", "dpx_median",
    ]:
        if col in merged.columns:
            merged[col] = pd.to_numeric(merged[col], errors="coerce")

    merged["peptide_bound_bool"] = merged["peptide_bound"].astype(str).eq("True")
    merged["basic_amine_bool"] = merged["ligand_basic_amine"].astype(str).eq("True")
    merged["tls_refined_bool"] = merged["tls_refined"].astype(str).eq("True")

    # Audit unmatched curated real pairs & unmatched chains
    matched_pairs = set(zip(
        merged.loc[merged["curated_matched_primary"], "pdb_norm"],
        merged.loc[merged["curated_matched_primary"], "lig_norm"],
    ))
    repo_pdbs = set(ct["pdb_norm"])
    # Also check all ligand_codes (semicolon separated)
    all_expanded_pairs = set()
    for _, r in ct.iterrows():
        for code in str(r["ligand_codes"] or "").split(";"):
            code = code.strip().upper()
            if code:
                all_expanded_pairs.add((r["pdb_norm"], code))

    unmatched_rows = []
    for _, r in cur.iterrows():
        key = (r["pdb_norm"], r["lig_norm"])
        if not r["is_real_bool"]:
            reason = f"excluded_non_ligand ({r['action_class']})"
            status = "excluded_non_ligand"
        elif key in matched_pairs:
            sub = merged[(merged["pdb_norm"] == key[0]) & (merged["lig_norm"] == key[1])]
            if (sub["ligand_status"] == "single").all():
                reason = "included_strict_single"
                status = "included_strict"
            else:
                reason = f"matched_primary_but_multi_ligand (statuses={sorted(sub['ligand_status'].unique())})"
                status = "matched_multi"
        else:
            if key[0] not in repo_pdbs:
                reason = "pdb_not_in_788_chain_table"
                status = "unmatched_missing_pdb"
            elif key in all_expanded_pairs:
                reason = "secondary_ccd_in_multi_ligand_chain"
                status = "unmatched_secondary_ccd"
            else:
                reason = "pdb_present_but_ccd_not_in_chain_table"
                status = "unmatched_ccd_mismatch"
        unmatched_rows.append({
            "pdb_id": r["pdb_norm"], "ligand_id": r["lig_norm"],
            "action_class": r["action_class"], "is_real_ligand": r["is_real_bool"],
            "doi": r["doi"], "status": status, "reason": reason,
        })
    # Also chains in chain_table with no curated match
    unmatched_chains = merged[~merged["curated_matched_primary"]][
        ["chain_key", "pdb_norm", "lig_norm", "ligand_status", "ligand_codes", "group"]
    ].copy()
    unmatched_chains["reason"] = np.where(
        unmatched_chains["lig_norm"] == "", "chain_without_primary_ligand", "primary_ccd_not_in_curated_real"
    )
    audit_df = pd.DataFrame(unmatched_rows)
    audit_df.to_csv(cc.TABLES_DIR / "pair_inclusion_audit.csv", index=False)
    unmatched_chains.to_csv(cc.TABLES_DIR / "unmatched_chains_audit.csv", index=False)

    # Descriptive audit of expanded ligand_codes (';')
    exp_rows = []
    for _, r in ct[ct["ligand_status"] == "multi"].iterrows():
        codes = [x.strip().upper() for x in str(r["ligand_codes"]).split(";") if x.strip()]
        for idx_c, code in enumerate(codes):
            hit = cur_real[(cur_real["pdb_norm"] == r["pdb_norm"]) & (cur_real["lig_norm"] == code)]
            exp_rows.append({
                "chain_key": r["chain_key"], "pdb_norm": r["pdb_norm"], "group": r["group"],
                "ligand_id_primary": r["lig_norm"], "expanded_ccd": code, "is_primary": code == r["lig_norm"],
                "n_codes_in_chain": len(codes),
                "curated_action_class": hit["action_class"].iloc[0] if len(hit) else "unmatched",
            })
    pd.DataFrame(exp_rows).to_csv(cc.TABLES_DIR / "multi_ligand_expanded_audit.csv", index=False)

    # Build pair-level (PDB x CCD) table on strict single matched chains
    strict_df = merged[merged["in_strict_single"]].copy()
    pair_rows = []
    for (pdb, lig), g in strict_df.groupby(["pdb_norm", "lig_norm"], sort=True):
        groups_present = sorted(g["group"].unique())
        if len(groups_present) == 1:
            g_state = f"{groups_present[0]} only"
            g_excl = groups_present[0]
        else:
            g_state = "multiple (" + "+".join(groups_present) + ")"
            g_excl = "multiple"
        r0 = g.iloc[0]
        pair_rows.append({
            "pdb_norm": pdb,
            "lig_norm": lig,
            "pair_key": f"{pdb}|{lig}",
            "action_class": r0["action_class"],
            "drug_classes": r0["drug_classes"],
            "old_chembl_class": r0["old_chembl_class"],
            "doi": r0["doi"],
            "doi_block": r0["doi_block"],
            "inchikey": r0["inchikey"],
            "repo_inchikey": r0["repo_inchikey"],
            "primary_source": r0["primary_source"],
            "evidence_sources": r0["evidence_sources"],
            "luc_efficacy_pct": r0["luc_efficacy_pct"],
            "proliferation_pct": r0["proliferation_pct"],
            "n_chains": int(len(g)),
            "groups_present": ";".join(groups_present),
            "h12_group_state": g_state,
            "h12_group_excl": g_excl,
            "has_cl7": "cl7" in groups_present,
            "has_cl4": "cl4" in groups_present,
            "has_minor": "minor" in groups_present,
            "construct_mutations": r0["construct_mutations"],
            "mutation_status": r0["mutation_status"],
            "mutations_norm_csv": r0["mutations_norm"],
            "mutation_group_csv": r0["mutation_group"],
            "is_mutant_csv": bool(r0["is_mutant_bool"]),
            "peptide_bound_any": bool(g["peptide_bound_bool"].any()),
            "peptide_bound_all": bool(g["peptide_bound_bool"].all()),
            "has_peptide_csv": bool(r0["has_peptide_bool"]),
            "qc_label": r0["qc_label"],
            "qc_pass": r0["qc_pass"],
            "resolution_A": float(g["resolution_A"].median()),
            "r_free": float(g["r_free"].median()),
            "dist_H12c_H11axis_A": float(g["dist_H12c_H11axis_A"].median()),
            "angle_H12_H3_deg": float(g["angle_H12_H3_deg"].median()),
            "angle_H12_H11_deg": float(g["angle_H12_H11_deg"].median()),
            "len_H12": float(g["len_H12"].median()),
            "mobility_h12_z": float(g["mobility_h12_z"].median()),
            "h12_complete_547_share": float(g["h12_complete_547"].mean()),
            "h12_complete_544_share": float(g["h12_complete_544"].mean()),
            "sym_atoms_median": float(g["sym_atoms"].median()),
            "ld_heavy_atoms": r0["ld_heavy_atoms"],
            "mw": r0["mw"],
            "clogp": r0["clogp"],
            "tpsa": r0["tpsa"],
            "basic_amine": bool(r0["basic_amine_bool"]),
            "murcko": r0["murcko"],
            "smiles_stereo": r0["smiles_stereo"],
            "repo_ligand_name": r0["repo_ligand_name"],
            "note": r0["note"],
            "pdb_title": r0["pdb_title"],
        })
    pairs_df = pd.DataFrame(pair_rows)

    # Pair x H12-group level table (for within-conformation analyses, median over chains)
    pg_rows = []
    for (pdb, lig, grp), g in strict_df.groupby(["pdb_norm", "lig_norm", "group"], sort=True):
        r0 = g.iloc[0]
        # Note-based check for WT assay on mutant structure
        note_str = str(r0["note"] or "")
        wt_assay_on_mut = bool(
            ("WT" in note_str or "wild-type" in note_str.lower() or "wild type" in note_str.lower())
            and (
                "mutant" in note_str.lower()
                or "Y537S" in note_str
                or "L536S" in note_str
                or "D538G" in note_str
                or "mismatch" in note_str.lower()
            )
        )
        pg_rows.append({
            "pdb_norm": pdb,
            "lig_norm": lig,
            "pair_key": f"{pdb}|{lig}",
            "group": grp,
            "action_class": r0["action_class"],
            "drug_classes": r0["drug_classes"],
            "old_chembl_class": r0["old_chembl_class"],
            "doi": r0["doi"],
            "doi_block": r0["doi_block"],
            "inchikey": r0["inchikey"],
            "primary_source": r0["primary_source"],
            "luc_efficacy_pct": r0["luc_efficacy_pct"],
            "proliferation_pct": r0["proliferation_pct"],
            "n_chains": int(len(g)),
            "construct_mutations": r0["construct_mutations"],
            "mutation_status": r0["mutation_status"],
            "mutations_norm_csv": r0["mutations_norm"],
            "mutation_group_csv": r0["mutation_group"],
            "peptide_bound": bool(g["peptide_bound_bool"].any()),
            "has_peptide_csv": bool(r0["has_peptide_bool"]),
            "qc_label": r0["qc_label"],
            "qc_pass": r0["qc_pass"],
            "resolution_A": float(g["resolution_A"].median()),
            "r_free": float(g["r_free"].median()),
            "tls_refined": bool(g["tls_refined_bool"].any()),
            "dist_H12c_H11axis_A": float(g["dist_H12c_H11axis_A"].median()),
            "dist_H12c_H11axis_A_w544": float(g["dist_H12c_H11axis_A_w544"].median()),
            "angle_H12_H3_deg": float(g["angle_H12_H3_deg"].median()),
            "angle_H12_H11_deg": float(g["angle_H12_H11_deg"].median()),
            "len_H12": float(g["len_H12"].median()),
            "mobility_h12_z": float(g["mobility_h12_z"].median()),
            "n_obs_538_547": float(g["n_obs_538_547"].median()),
            "n_obs_538_544": float(g["n_obs_538_544"].median()),
            "h12_complete_547": bool(g["h12_complete_547"].all()),
            "h12_complete_544": bool(g["h12_complete_544"].all()),
            "sym_atoms": float(g["sym_atoms"].median()),
            "pocket_volume_A3": float(g["pocket_volume_A3"].median()),
            "ligand_volume_A3": float(g["ligand_volume_A3"].median()),
            "filling_pct": float(g["filling_pct"].median()),
            "pocket_pct_of_E2_pocket": float(g["pocket_pct_of_E2_pocket"].median()),
            "pocket_complete_primary": bool(
                (g["pocket_in_primary"].astype(str) == "True").all()
                & (g["pocket_ligand_incomplete"].astype(str) == "False").all()
            ),
            "cl4_sub_cluster": g["cl4_sub_cluster"].dropna().iloc[0] if g["cl4_sub_cluster"].notna().any() else np.nan,
            "h12_shift_A": float(g["h12_shift_A"].median()),
            "k362_e542_A": float(g["k362_e542_A"].median()),
            "wt_assay_on_mutant_note": wt_assay_on_mut,
            "note": note_str,
        })
    pg_df = pd.DataFrame(pg_rows)

    # Construct & peptide reconciliation
    recon = pairs_df[[
        "pdb_norm", "lig_norm", "action_class", "construct_mutations", "mutation_status",
        "mutations_norm_csv", "is_mutant_csv", "peptide_bound_any", "has_peptide_csv",
    ]].copy()
    recon["struct_is_mutant"] = recon["mutation_status"] != "wild type"
    recon["mutant_flag_mismatch"] = recon["struct_is_mutant"] != recon["is_mutant_csv"]
    recon["peptide_flag_mismatch"] = recon["peptide_bound_any"] != recon["has_peptide_csv"]
    recon.to_csv(cc.TABLES_DIR / "construct_peptide_reconciliation.csv", index=False)

    # InChIKey conflicts (same InChIKey with >1 action_class in curated real)
    ik_cls = cur_real.groupby("inchikey")["action_class"].nunique()
    conflict_iks = ik_cls[ik_cls > 1].index.tolist()
    conflicts = cur_real[cur_real["inchikey"].isin(conflict_iks)][
        ["inchikey", "pdb_norm", "lig_norm", "action_class", "primary_source", "mutations_norm", "has_peptide", "doi", "note"]
    ].sort_values(["inchikey", "action_class", "pdb_norm"])
    conflicts.to_csv(cc.TABLES_DIR / "inchikey_conflicts.csv", index=False)

    # WT-assay on mutant-structure flags
    wt_flags = pg_df[pg_df["wt_assay_on_mutant_note"]][
        ["pdb_norm", "lig_norm", "group", "action_class", "mutation_status", "primary_source", "luc_efficacy_pct", "proliferation_pct", "note"]
    ]
    wt_flags.to_csv(cc.TABLES_DIR / "wt_assay_mutant_flags.csv", index=False)

    # Save master tables
    merged.to_csv(cc.TABLES_DIR / "chains_curated.csv", index=False)
    pairs_df.to_csv(cc.TABLES_DIR / "pairs_curated.csv", index=False)
    pg_df.to_csv(cc.TABLES_DIR / "pairs_by_group.csv", index=False)

    core3_pairs = pairs_df[pairs_df["action_class"].isin(cc.CLASS_ORDER_CORE)]
    core3_chains = merged[merged["in_core3_strict"]]
    summary = {
        "git_sha": get_git_sha(),
        "dataset_path": str(dataset_path),
        "curated_total_rows": int(len(cur)),
        "curated_real_rows": int(len(cur_real)),
        "curated_real_action_counts": cur_real["action_class"].value_counts().to_dict(),
        "curated_core3_rows": int(cur_real["action_class"].isin(cc.CLASS_ORDER_CORE).sum()),
        "curated_drug_classes_filled_real": int(cur_real["drug_classes"].notna().sum()),
        "curated_unknown_serd": int(
            ((cur_real["action_class"] == "unknown") & cur_real["drug_classes"].fillna("").str.contains("SERD")).sum()
        ),
        "chain_table_chains": int(len(ct)),
        "chain_table_pdbs": int(ct["pdb_norm"].nunique()),
        "chain_table_groups": ct["group"].value_counts().to_dict(),
        "exact_primary_match_chains": int(merged["curated_matched_primary"].sum()),
        "exact_primary_match_pairs": int(merged.loc[merged["curated_matched_primary"], ["pdb_norm", "lig_norm"]].drop_duplicates().shape[0]),
        "strict_single_chains": int(merged["in_strict_single"].sum()),
        "strict_single_pairs": int(len(pairs_df)),
        "strict_single_pdbs": int(pairs_df["pdb_norm"].nunique()),
        "core3_strict_chains": int(len(core3_chains)),
        "core3_strict_pairs": int(len(core3_pairs)),
        "core3_strict_pdbs": int(core3_pairs["pdb_norm"].nunique()),
        "core3_strict_inchikeys": int(core3_pairs["inchikey"].nunique()),
        "core3_strict_by_class": core3_pairs["action_class"].value_counts().to_dict(),
        "pair_audit_status_counts": audit_df["status"].value_counts().to_dict(),
        "unmatched_chains_count": int(len(unmatched_chains)),
        "mutant_flag_mismatches": int(recon["mutant_flag_mismatch"].sum()),
        "peptide_flag_mismatches": int(recon["peptide_flag_mismatch"].sum()),
        "inchikey_conflict_count": int(len(conflict_iks)),
        "wt_assay_on_mutant_flagged_rows": int(len(wt_flags)),
    }
    cc.write_json(cc.TABLES_DIR / "join_audit_summary.json", summary)
    print("=== Stage 12.1 Join & Audit Summary ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
