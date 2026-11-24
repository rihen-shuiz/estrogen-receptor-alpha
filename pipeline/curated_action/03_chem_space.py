"""Stage 12.3 - Chemical space (RDKit ECFP4 & MACCS), Tanimoto matrices, clustering stability, structural overlay, and discordant close pairs."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.manifold import MDS
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402


def tanimoto_matrix(fps) -> np.ndarray:
    n = len(fps)
    m = np.eye(n, dtype=float)
    for i in range(n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i + 1:])
        m[i, i + 1:] = sims
        m[i + 1:, i] = sims
    return m


def cramers_v_perm(a, b, n_perm: int = 1000, seed: int = cc.SEED) -> tuple[float, float]:
    tab = pd.crosstab(a, b)
    if min(tab.shape) < 2:
        return (np.nan, np.nan)
    chi2 = chi2_contingency(tab, correction=False)[0]
    n = tab.values.sum()
    v = float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))
    rng = np.random.default_rng(seed)
    a_arr, b_arr = np.asarray(a), np.asarray(b)
    ge = 0
    for _ in range(n_perm):
        bp = rng.permutation(b_arr)
        tp = pd.crosstab(a_arr, bp)
        c2 = chi2_contingency(tp, correction=False)[0]
        if c2 >= chi2 - 1e-12:
            ge += 1
    return (v, (ge + 1) / (n_perm + 1))


def main():
    cc.ensure_dirs()
    pairs = cc.read_csv_safe(cc.TABLES_DIR / "pairs_curated.csv")
    core3 = pairs[pairs["action_class"].isin(cc.CLASS_ORDER_CORE)].copy()

    # Check SMILES availability & InChIKey match
    n_smiles_found = int(core3["smiles_stereo"].notna().sum())
    ik_match = (core3["inchikey"].fillna("") == core3["repo_inchikey"].fillna("")).sum()

    # Aggregate to unique CCD ligands across core3 strict pairs (preserving context for repeats)
    lig_rows = []
    for ccd, g in core3.groupby("lig_norm", sort=True):
        smi = g["smiles_stereo"].dropna().iloc[0] if g["smiles_stereo"].notna().any() else ""
        cls_counts = g["action_class"].value_counts()
        grp_counts = g["h12_group_excl"].value_counts()
        cls_list = sorted(g["action_class"].unique())
        grp_list = sorted(g["h12_group_excl"].unique())
        lig_rows.append({
            "lig_norm": ccd,
            "inchikey": g["inchikey"].iloc[0],
            "repo_inchikey": g["repo_inchikey"].iloc[0],
            "inchikey_exact_match": g["inchikey"].iloc[0] == g["repo_inchikey"].iloc[0],
            "smiles_stereo": smi,
            "ligand_name": g["repo_ligand_name"].dropna().iloc[0] if g["repo_ligand_name"].notna().any() else ccd,
            "n_pairs": int(len(g)),
            "n_pdbs": int(g["pdb_norm"].nunique()),
            "n_dois": int(g["doi_block"].nunique()),
            "pdbs": ";".join(sorted(g["pdb_norm"].unique())),
            "action_class_majority": cls_counts.index[0],
            "action_class_all": ";".join(cls_list),
            "action_conflict": len(cls_list) > 1,
            "h12_group_majority": grp_counts.index[0],
            "h12_group_all": ";".join(grp_list),
            "h12_conflict": len(grp_list) > 1,
            "heavy_atoms": float(pd.to_numeric(g["ld_heavy_atoms"], errors="coerce").median()),
            "basic_amine": bool(g["basic_amine"].astype(str).eq("True").any()),
            "murcko": g["murcko"].dropna().iloc[0] if g["murcko"].notna().any() else "",
        })
    lig_df = pd.DataFrame(lig_rows).sort_values("lig_norm").reset_index(drop=True)

    # Compute RDKit fingerprints
    ecfp4_fps, maccs_fps, valid_mask = [], [], []
    for smi in lig_df["smiles_stereo"]:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) and smi else None
        if mol is None:
            valid_mask.append(False)
            ecfp4_fps.append(None)
            maccs_fps.append(None)
        else:
            valid_mask.append(True)
            ecfp4_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
            maccs_fps.append(MACCSkeys.GenMACCSKeys(mol))

    assert all(valid_mask), f"Invalid SMILES for {lig_df.loc[~np.array(valid_mask), 'lig_norm'].tolist()}"
    sim_ecfp4 = tanimoto_matrix(ecfp4_fps)
    sim_maccs = tanimoto_matrix(maccs_fps)

    # Save Tanimoto matrices (both npy and csv)
    np.save(cc.TABLES_DIR / "tanimoto_ecfp4.npy", sim_ecfp4)
    np.save(cc.TABLES_DIR / "tanimoto_maccs.npy", sim_maccs)
    pd.DataFrame(sim_ecfp4, index=lig_df["lig_norm"], columns=lig_df["lig_norm"]).round(4).to_csv(
        cc.TABLES_DIR / "tanimoto_ecfp4.csv"
    )
    pd.DataFrame(sim_maccs, index=lig_df["lig_norm"], columns=lig_df["lig_norm"]).round(4).to_csv(
        cc.TABLES_DIR / "tanimoto_maccs.csv"
    )

    # Hierarchical clustering (average linkage on 1 - Tanimoto) at multiple cutoffs
    dist_ecfp4 = np.clip(1.0 - sim_ecfp4, 0.0, 1.0)
    dist_maccs = np.clip(1.0 - sim_maccs, 0.0, 1.0)
    z_ecfp4 = linkage(squareform(dist_ecfp4, checks=False), method="average")
    z_maccs = linkage(squareform(dist_maccs, checks=False), method="average")

    cut_cols = []
    for name, z in (("ecfp4", z_ecfp4), ("maccs", z_maccs)):
        for cut in (0.4, 0.5, 0.6):
            col = f"cl_{name}_{str(cut).replace('.', 'p')}"
            lab = fcluster(z, t=cut, criterion="distance")
            # Rank cluster IDs by descending size
            order = pd.Series(lab).value_counts().index.tolist()
            rmap = {o: i + 1 for i, o in enumerate(order)}
            lig_df[col] = [rmap[x] for x in lab]
            cut_cols.append(col)

    # Stability matrix (ARI between cutoffs and fingerprints)
    stab_rows = []
    for c1 in cut_cols:
        for c2 in cut_cols:
            stab_rows.append({
                "clustering_a": c1,
                "clustering_b": c2,
                "n_clusters_a": int(lig_df[c1].nunique()),
                "n_clusters_b": int(lig_df[c2].nunique()),
                "ari": float(adjusted_rand_score(lig_df[c1], lig_df[c2])),
            })
    stab_df = pd.DataFrame(stab_rows)
    stab_df.to_csv(cc.TABLES_DIR / "chem_cluster_stability.csv", index=False)

    # 2D MDS embedding on ECFP4 distance (seed 42)
    mds = MDS(n_components=2, dissimilarity="precomputed", random_state=cc.SEED, normalized_stress="auto", n_init=8)
    xy = mds.fit_transform(dist_ecfp4)
    lig_df["mds_x"] = xy[:, 0]
    lig_df["mds_y"] = xy[:, 1]
    lig_df.to_csv(cc.TABLES_DIR / "chem_ligands_order.csv", index=False)

    # Merge primary chemical cluster (cl_ecfp4_0p6) back to pair level for contingency & association tests
    primary_cl_col = "cl_ecfp4_0p6"
    pairs_chem = core3.merge(lig_df[["lig_norm", primary_cl_col, "cl_ecfp4_0p5", "cl_maccs_0p4", "mds_x", "mds_y"]],
                             on="lig_norm", how="left", validate="many_to_one")
    pairs_chem.to_csv(cc.TABLES_DIR / "pairs_core3_with_chem.csv", index=False)

    # Association strength (Cramer's V + permutation p) at 3 levels: pairs, unique ligands, DOI-deduplicated
    assoc_rows = []
    # Restrict to clusters with >= 3 items for clean contingency or pool small clusters as "small"
    for level_name, df_lvl, grp_col, cls_col in [
        ("pairs_core3", pairs_chem, "h12_group_excl", "action_class"),
        ("unique_ligands", lig_df, "h12_group_majority", "action_class_majority"),
        ("doi_dedup_pairs", pairs_chem.sort_values("pdb_norm").drop_duplicates(["doi_block", primary_cl_col]), "h12_group_excl", "action_class"),
    ]:
        # Pool clusters with < 4 members into "other_small"
        vc = df_lvl[primary_cl_col].value_counts()
        big = set(vc[vc >= 4].index)
        cl_pooled = df_lvl[primary_cl_col].map(lambda x: f"C{x}" if x in big else "C_small")
        for target_name, tcol in [("h12_group", grp_col), ("action_class", cls_col)]:
            v, pv = cramers_v_perm(cl_pooled, df_lvl[tcol], n_perm=1000, seed=cc.SEED)
            assoc_rows.append({
                "level": level_name,
                "chem_clustering": primary_cl_col,
                "target": target_name,
                "n_obs": int(len(df_lvl)),
                "n_big_clusters": int(len(big)),
                "cramers_v": v,
                "perm_p": pv,
            })
    assoc_df = pd.DataFrame(assoc_rows)
    assoc_df.to_csv(cc.TABLES_DIR / "chem_struct_class_associations.csv", index=False)

    # Find close ligand pairs (ECFP4 Tanimoto >= 0.55) with discordant action_class or H12 group
    disc_rows = []
    n_lig = len(lig_df)
    for i in range(n_lig):
        for j in range(i + 1, n_lig):
            tc = float(sim_ecfp4[i, j])
            if tc >= 0.55:
                ri, rj = lig_df.iloc[i], lig_df.iloc[j]
                diff_cls = ri["action_class_majority"] != rj["action_class_majority"]
                diff_grp = ri["h12_group_majority"] != rj["h12_group_majority"]
                if diff_cls or diff_grp:
                    pi = core3[core3["lig_norm"] == ri["lig_norm"]].iloc[0]
                    pj = core3[core3["lig_norm"] == rj["lig_norm"]].iloc[0]
                    disc_rows.append({
                        "lig_a": ri["lig_norm"], "lig_b": rj["lig_norm"],
                        "tanimoto_ecfp4": round(tc, 4),
                        "tanimoto_maccs": round(float(sim_maccs[i, j]), 4),
                        "diff_class": diff_cls, "diff_h12_group": diff_grp,
                        "class_a": ri["action_class_majority"], "class_b": rj["action_class_majority"],
                        "group_a": ri["h12_group_majority"], "group_b": rj["h12_group_majority"],
                        "pdb_a": pi["pdb_norm"], "pdb_b": pj["pdb_norm"],
                        "construct_a": pi["construct_mutations"], "construct_b": pj["construct_mutations"],
                        "peptide_a": pi["peptide_bound_any"], "peptide_b": pj["peptide_bound_any"],
                        "source_a": pi["primary_source"], "source_b": pj["primary_source"],
                        "doi_a": pi["doi"], "doi_b": pj["doi"],
                    })
    disc_df = pd.DataFrame(disc_rows).sort_values("tanimoto_ecfp4", ascending=False)
    disc_df.to_csv(cc.TABLES_DIR / "chem_discordant_pairs.csv", index=False)

    summary = {
        "core3_pairs": int(len(core3)),
        "smiles_found_pairs": n_smiles_found,
        "inchikey_exact_match_pairs": int(ik_match),
        "unique_ccd_ligands": int(len(lig_df)),
        "unique_inchikeys": int(lig_df["inchikey"].nunique()),
        "action_conflicts_same_ccd": int(lig_df["action_conflict"].sum()),
        "h12_conflicts_same_ccd": int(lig_df["h12_conflict"].sum()),
        "mds_stress": float(mds.stress_),
        "ecfp4_0p6_n_clusters": int(lig_df["cl_ecfp4_0p6"].nunique()),
        "ari_ecfp4_0p6_vs_maccs_0p4": float(adjusted_rand_score(lig_df["cl_ecfp4_0p6"], lig_df["cl_maccs_0p4"])),
        "n_discordant_close_pairs_tc_ge_0p55": int(len(disc_df)),
        "top_discordant_pairs": disc_df.head(6).to_dict("records"),
        "associations": assoc_df.to_dict("records"),
    }
    cc.write_json(cc.TABLES_DIR / "chem_space_summary.json", summary)
    print("=== Stage 12.3 Chemical Space Summary ===")
    print(f"Unique CCDs: {len(lig_df)}, SMILES found: {n_smiles_found}/{len(core3)}, InChIKey exact match: {ik_match}/{len(core3)}")
    print(assoc_df.to_string(index=False))
    print("Top discordant pairs:\n", disc_df.head(5)[["lig_a", "lig_b", "tanimoto_ecfp4", "class_a", "class_b", "group_a", "group_b", "pdb_a", "pdb_b"]].to_string(index=False))


if __name__ == "__main__":
    main()
