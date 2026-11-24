"""Stage 12.2 - Redundancy accumulation curves by curated class, representatives/medoids, and QC clustering sensitivity."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common_curated as cc  # noqa: E402

CUTS = ("1.0", "1.5", "2.0", "2.5")
N_PERM = 200


def pdb_cluster_sets(df: pd.DataFrame, cut_col: str) -> dict[str, set]:
    """Map each PDB to the set of cluster IDs observed among its eligible chains."""
    return {pdb: set(g[cut_col].astype(int).tolist()) for pdb, g in df.groupby("pdb_norm")}


def accumulation_curve_pdb(pdb_to_clusters: dict[str, set], n_perm: int = N_PERM, seed: int = cc.SEED) -> pd.DataFrame:
    """Random PDB addition without replacement: unique clusters discovered after k PDBs."""
    pdbs = sorted(pdb_to_clusters.keys())
    n = len(pdbs)
    if n == 0:
        return pd.DataFrame()
    all_clusters = set().union(*pdb_to_clusters.values())
    total_c = len(all_clusters)
    sets_list = [pdb_to_clusters[p] for p in pdbs]
    rng = np.random.default_rng(seed)
    mat = np.zeros((n_perm, n), dtype=float)
    for r in range(n_perm):
        order = rng.permutation(n)
        seen = set()
        for k, idx in enumerate(order):
            seen |= sets_list[idx]
            mat[r, k] = len(seen)
    mean = mat.mean(axis=0)
    sd = mat.std(axis=0, ddof=1)
    p10 = np.percentile(mat, 10, axis=0)
    p90 = np.percentile(mat, 90, axis=0)
    n90 = int(np.argmax(mean >= 0.9 * total_c)) + 1 if total_c > 0 else np.nan
    return pd.DataFrame({
        "n_pdbs": np.arange(1, n + 1),
        "mean_clusters": mean,
        "sd_clusters": sd,
        "p10_clusters": p10,
        "p90_clusters": p90,
        "mean_frac_total": mean / total_c if total_c > 0 else np.nan,
        "total_clusters_denom": total_c,
        "total_pdbs": n,
        "n90_pdbs": n90,
    })


def pdb_weighted_hill(pdb_to_clusters: dict[str, set]) -> dict:
    """Equal weight (1.0) per PDB; split equally across clusters if a PDB spans multiple clusters."""
    weights: dict[int, float] = {}
    for cl_set in pdb_to_clusters.values():
        w = 1.0 / len(cl_set)
        for c_id in cl_set:
            weights[c_id] = weights.get(c_id, 0.0) + w
    return cc.hill_numbers(list(weights.values()))


def greedy_cover_detailed(d: np.ndarray, chains_df: pd.DataFrame, radius: float) -> pd.DataFrame:
    cand = (chains_df["h12_label"] == "complete").to_numpy()
    cand_idx = np.where(cand)[0]
    within = d[:, cand_idx] <= radius
    covered = np.zeros(len(d), dtype=bool)
    rows = []
    step = 0
    while covered.mean() < 1.0:
        gain = within[~covered].sum(axis=0)
        k = int(np.argmax(gain))
        if gain[k] == 0:
            break
        step += 1
        covered |= within[:, k]
        ci = cand_idx[k]
        r = chains_df.iloc[ci]
        rows.append({
            "radius_A": radius,
            "step": step,
            "chain_key": r["chain_key"],
            "pdb_norm": r["pdb_norm"],
            "qc_label": r["qc_label"],
            "cluster_1p5": int(r["cluster_1p5"]),
            "group": r["group"],
            "action_class": r["action_class"] if pd.notna(r["action_class"]) else "unmatched",
            "gain_chains": int(gain[k]),
            "cum_covered_chains": int(covered.sum()),
            "cum_coverage": float(covered.mean()),
        })
    return pd.DataFrame(rows)


def main():
    cc.ensure_dirs()
    ch = cc.read_csv_safe(cc.TABLES_DIR / "chains_curated.csv")
    d_rmsd = np.load(cc.RESULTS_BASE / "02_lid_groups/rmsd.npy")
    idx_rmsd = cc.read_csv_safe(cc.RESULTS_BASE / "02_lid_groups/chains.csv")
    assert (ch["chain_key"].values == idx_rmsd["chain_key"].values).all(), "Chain order mismatch vs RMSD matrix"

    # ------------------------------------------------------------------ 1. Accumulation curves & Hill numbers
    subsets_def = {
        "all_351_pdbs": ch,
        "strict_single_all": ch[ch["in_strict_single"].astype(str) == "True"],
        "agonist": ch[(ch["in_core3_strict"].astype(str) == "True") & (ch["action_class"] == "agonist")],
        "partial_agonist": ch[(ch["in_core3_strict"].astype(str) == "True") & (ch["action_class"] == "partial_agonist")],
        "antagonist": ch[(ch["in_core3_strict"].astype(str) == "True") & (ch["action_class"] == "antagonist")],
    }
    # Sensitivity deduplicated subsets (on core3 strict)
    core3 = ch[ch["in_core3_strict"].astype(str) == "True"].copy()
    # One PDB per unique InChIKey (deterministic: first sorted pdb_norm)
    ik_pdbs = core3.sort_values("pdb_norm").drop_duplicates(["action_class", "inchikey"])["pdb_norm"].unique()
    doi_pdbs = core3.sort_values("pdb_norm").drop_duplicates(["action_class", "doi_block"])["pdb_norm"].unique()
    tierA_core3 = core3[core3["qc_label"] == "A"]
    sens_def = {
        "core3_all": core3,
        "core3_dedup_inchikey": core3[core3["pdb_norm"].isin(ik_pdbs)],
        "core3_dedup_doi": core3[core3["pdb_norm"].isin(doi_pdbs)],
        "core3_tierA": tierA_core3,
    }

    curve_rows = []
    hill_rows = []
    for sub_name, sub_df in {**subsets_def, **sens_def}.items():
        for cut in CUTS:
            cut_col = f"cut_{cut}"
            p2c = pdb_cluster_sets(sub_df, cut_col)
            cdf = accumulation_curve_pdb(p2c, n_perm=N_PERM, seed=cc.SEED)
            cdf["subset"] = sub_name
            cdf["cut_A"] = float(cut)
            curve_rows.append(cdf)
            hn = pdb_weighted_hill(p2c)
            n90 = int(cdf["n90_pdbs"].iloc[0]) if len(cdf) else np.nan
            hill_rows.append({
                "subset": sub_name,
                "cut_A": float(cut),
                "n_pdbs": len(p2c),
                "n_chains": len(sub_df),
                "n_inchikeys": int(sub_df["inchikey"].dropna().nunique()),
                "n_doi_blocks": int(sub_df["doi_block"].nunique()),
                "q0_clusters": hn["q0"],
                "shannon_bits": hn["shannon_bits"],
                "hill_q1": hn["q1"],
                "hill_q2": hn["q2"],
                "n90_pdbs": n90,
                "n90_frac_pdbs": n90 / len(p2c) if len(p2c) else np.nan,
            })

    curves_all = pd.concat(curve_rows, ignore_index=True)
    hill_all = pd.DataFrame(hill_rows)
    curves_all.to_csv(cc.TABLES_DIR / "redundancy_accumulation_curves.csv", index=False)
    hill_all.to_csv(cc.TABLES_DIR / "redundancy_hill_summary.csv", index=False)

    # ------------------------------------------------------------------ 2. Greedy set cover & cluster medoids
    sc_dfs = []
    sc_summary = []
    for r in (0.5, 0.75, 1.0, 1.5):
        gdf = greedy_cover_detailed(d_rmsd, ch, r)
        sc_dfs.append(gdf)
        hit90 = gdf[gdf["cum_coverage"] >= 0.90 - 1e-12]
        hit95 = gdf[gdf["cum_coverage"] >= 0.95 - 1e-12]
        sc_summary.append({
            "radius_A": r,
            "reps_90pct": int(hit90["step"].iloc[0]) if len(hit90) else np.nan,
            "reps_95pct": int(hit95["step"].iloc[0]) if len(hit95) else np.nan,
            "max_coverage": float(gdf["cum_coverage"].iloc[-1]) if len(gdf) else 0.0,
        })
    sc_all = pd.concat(sc_dfs, ignore_index=True)
    sc_all.to_csv(cc.TABLES_DIR / "greedy_set_cover_steps.csv", index=False)
    sc_sum_df = pd.DataFrame(sc_summary)
    sc_sum_df.to_csv(cc.TABLES_DIR / "greedy_set_cover_summary.csv", index=False)

    # Cluster medoids at 1.5 A (min sum of distances within cluster)
    medoid_rows = []
    for cl_id, g in ch.groupby("cluster_1p5"):
        idxs = g.index.to_numpy()
        sub_d = d_rmsd[np.ix_(idxs, idxs)]
        sums = sub_d.sum(axis=1)
        best_local = int(np.argmin(sums))
        best_idx = idxs[best_local]
        r = ch.iloc[best_idx]
        mean_d = float(sums[best_local] / max(len(idxs) - 1, 1))
        medoid_rows.append({
            "role": "cluster_medoid_1p5",
            "cluster_1p5": int(cl_id),
            "group": r["group"],
            "cluster_size_chains": int(len(idxs)),
            "cluster_size_pdbs": int(g["pdb_norm"].nunique()),
            "chain_key": r["chain_key"],
            "pdb_norm": r["pdb_norm"],
            "lig_norm": r["lig_norm"],
            "action_class": r["action_class"] if pd.notna(r["action_class"]) else "unmatched",
            "qc_label": r["qc_label"],
            "h12_label": r["h12_label"],
            "mean_rmsd_to_cluster_mates_A": mean_d,
        })
    medoids_df = pd.DataFrame(medoid_rows).sort_values("cluster_size_chains", ascending=False)
    medoids_df.to_csv(cc.TABLES_DIR / "cluster_medoids_1p5.csv", index=False)

    # ------------------------------------------------------------------ 3. QC clustering sensitivity
    rfree_num = pd.to_numeric(ch["r_free"], errors="coerce")
    res_num = pd.to_numeric(ch["resolution_A"], errors="coerce")
    b_num = pd.to_numeric(ch["mean_b_entry_A2"], errors="coerce")

    qc_masks = {
        "full_788": np.ones(len(ch), dtype=bool),
        "rfree_le_0p30": (rfree_num.notna() & (rfree_num <= 0.30)).to_numpy(),
        "qc_all3_tierA": (ch["qc_label"] == "A").to_numpy(),
    }
    qc_missing_counts = {
        "missing_rfree_chains": int(rfree_num.isna().sum()),
        "missing_rfree_pdbs": int(ch.loc[rfree_num.isna(), "pdb_norm"].nunique()),
        "missing_resolution_chains": int(res_num.isna().sum()),
        "missing_mean_b_chains": int(b_num.isna().sum()),
        "qc_pass_all3_chains": int(qc_masks["qc_all3_tierA"].sum()),
        "qc_pass_all3_pdbs": int(ch.loc[qc_masks["qc_all3_tierA"], "pdb_norm"].nunique()),
    }

    qc_rows = []
    trans_rows = []
    for sub_name, mask in qc_masks.items():
        sub_ch = ch[mask].copy().reset_index(drop=True)
        sub_d = d_rmsd[np.ix_(mask, mask)]
        z = linkage(squareform(sub_d, checks=False), method="average")
        for cut in CUTS:
            new_lab = fcluster(z, t=float(cut), criterion="distance")
            orig_lab = sub_ch[f"cut_{cut}"].astype(int).to_numpy()
            ari = float(adjusted_rand_score(orig_lab, new_lab))
            # PDB-weighted Hill & n90
            sub_tmp = sub_ch.assign(new_cl=new_lab)
            p2c = pdb_cluster_sets(sub_tmp, "new_cl")
            hn = pdb_weighted_hill(p2c)
            cdf = accumulation_curve_pdb(p2c, n_perm=100, seed=cc.SEED)
            n90 = int(cdf["n90_pdbs"].iloc[0]) if len(cdf) else np.nan
            qc_rows.append({
                "subset": sub_name,
                "cut_A": float(cut),
                "n_chains": int(mask.sum()),
                "n_chains_excluded": int((~mask).sum()),
                "n_pdbs": int(sub_ch["pdb_norm"].nunique()),
                "n_clusters": int(len(np.unique(new_lab))),
                "ari_vs_full": ari,
                "shannon_bits_pdb": hn["shannon_bits"],
                "hill_q1_pdb": hn["q1"],
                "n90_pdbs": n90,
                "n90_frac_pdbs": n90 / len(p2c) if len(p2c) else np.nan,
            })
            if cut == "1.5":
                # Transition matrix between original group (cl7/cl4/minor) and new 1.5A clusters (ranked by size)
                order_new = pd.Series(new_lab).value_counts().index.tolist()
                rank_map = {old: f"new_rank{i+1} (n={int((new_lab==old).sum())})" for i, old in enumerate(order_new)}
                ctab = pd.crosstab(sub_ch["group"], pd.Series(new_lab).map(rank_map))
                for og in ctab.index:
                    for ng in ctab.columns:
                        if ctab.loc[og, ng] > 0:
                            trans_rows.append({
                                "subset": sub_name,
                                "orig_group": og,
                                "new_cluster_ranked": ng,
                                "n_chains": int(ctab.loc[og, ng]),
                            })

    qc_df = pd.DataFrame(qc_rows)
    qc_df.to_csv(cc.TABLES_DIR / "qc_clustering_sensitivity.csv", index=False)
    pd.DataFrame(trans_rows).to_csv(cc.TABLES_DIR / "qc_clustering_transitions_1p5.csv", index=False)

    summary = {
        "n_perm": N_PERM,
        "set_cover_check": sc_sum_df.to_dict("records"),
        "qc_missing_counts": qc_missing_counts,
        "hill_1p5_core_classes": hill_all[
            (hill_all["cut_A"] == 1.5) & (hill_all["subset"].isin(["all_351_pdbs", "agonist", "partial_agonist", "antagonist"]))
        ].to_dict("records"),
        "qc_ari_1p5": qc_df[qc_df["cut_A"] == 1.5][["subset", "n_chains", "n_clusters", "ari_vs_full", "hill_q1_pdb"]].to_dict("records"),
    }
    cc.write_json(cc.TABLES_DIR / "redundancy_qc_summary.json", summary)
    print("=== Stage 12.2 Redundancy & QC Sensitivity ===")
    print(sc_sum_df.to_string(index=False))
    print(qc_df[qc_df["cut_A"] == 1.5].to_string(index=False))


if __name__ == "__main__":
    main()
