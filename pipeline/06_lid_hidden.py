"""Stage 06 - hide the lid (H12) and cluster again: which groups appear and what separates them.

Input: the saved RMSD matrix without residues 528-550 (H11-H12 loop, H12, C-terminal tail; results/02_lid_groups/rmsd_without_528_550.npy,
same core-v3 superposition, positions observed in both chains of a pair).
Declared before running:
  * Average linkage (as for the 1.5 A groups). K = the value in 2..10 with the highest silhouette (precomputed
    distances); K = 2 is also reported for a direct comparison with the lid split.
  * Robustness: spectral clustering (affinity exp(-d^2 / 2 s^2), s = median distance, seed 42) at the same K;
    ARI between the two methods.
  * Which groups: composition by lid group, crystal form, construct, peptide, ligand basic amine, action class
    (uncurated ChEMBL); association = ARI of the new labels with each covariate at chain level and Cramer's V at entry level
    (one row per entry and new group).
  * What separates them: per-residue distance between group-mean C-alpha positions (superposed coordinates,
    528-550 excluded), top residues and per-helix sums of the squared difference.
  * Real only if the median between-group distance exceeds the median within-group distance by more than the 95th
    percentile of copy-to-copy distances (chains of one entry, same matrix).
Outputs: results/06_lid_hidden/no_lid_clusters.csv, no_lid_k_scan.csv, no_lid_composition.csv, no_lid_association.csv,
no_lid_residue_profile.csv, numbers_no_lid.json.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import chi2_contingency
from sklearn.cluster import SpectralClustering
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("06_lid_hidden")

HELICES = {"H1": (312, 321), "H3": (342, 361), "H4/5": (372, 394), "H7": (422, 437), "H8": (442, 455),
           "H9": (473, 492), "H10/11": (497, 526)}


def cramers_v(a, b) -> float:
    tab = pd.crosstab(a, b)
    if min(tab.shape) < 2:
        return float("nan")
    chi2 = chi2_contingency(tab, correction=False)[0]
    n = tab.values.sum()
    return float(np.sqrt(chi2 / (n * (min(tab.shape) - 1))))


def relabel_by_size(lab):
    order = pd.Series(lab).value_counts().index.tolist()
    m = {old: i + 1 for i, old in enumerate(order)}
    return np.array([m[x] for x in lab])


def main():
    t = c.chains()
    d, idx = c.rmsd("rmsd_without_lid")
    d = np.nan_to_num(d, nan=np.nanmedian(d))
    np.fill_diagonal(d, 0)
    z = linkage(squareform(d, checks=False), "average")
    scan = []
    for k in range(2, 11):
        lab = fcluster(z, k, "maxclust")
        if len(set(lab)) < 2:
            continue
        scan.append(dict(k=k, silhouette=silhouette_score(d, lab, metric="precomputed"),
                         sizes=";".join(map(str, sorted(pd.Series(lab).value_counts().tolist(), reverse=True))),
                         ari_vs_lid_groups=adjusted_rand_score(t.cluster_1p5, lab)))
    scan = pd.DataFrame(scan)
    scan.to_csv(OUT / "no_lid_k_scan.csv", index=False)
    kbest = int(scan.loc[scan.silhouette.idxmax(), "k"])
    # first balanced average-linkage split: smallest K whose second-largest group holds >= 5 % of the chains
    bal = scan[scan.sizes.map(lambda x: int(x.split(";")[1]) >= 0.05 * len(d))]
    kbal = int(bal.k.min()) if len(bal) else None
    s = np.median(d[np.triu_indices(len(d), 1)])
    aff = np.exp(-d ** 2 / (2 * s * s))
    labels = {
        "average_silhouette_best": relabel_by_size(fcluster(z, kbest, "maxclust")),
        "average_first_balanced": relabel_by_size(fcluster(z, kbal, "maxclust")) if kbal else None,
        "spectral_k2": relabel_by_size(SpectralClustering(2, affinity="precomputed",
                                                          random_state=c.SEED).fit_predict(aff)),
    }
    labels = {k: v for k, v in labels.items() if v is not None}
    out = t[["chain_key", "entry_id", "group", "cluster_1p5", "space_group", "mutation_status", "peptide_bound",
             "ligand_basic_amine", "action_class", "h12_label", "ligand_chemotype"]].copy()
    for k, v in labels.items():
        out[k] = v
    out.to_csv(OUT / "no_lid_clusters.csv", index=False)

    same_entry = t.entry_id.to_numpy()[:, None] == t.entry_id.to_numpy()[None, :]
    copy = d[np.triu(same_entry, 1)]
    copy_p95 = float(np.percentile(copy, 95))
    zc = np.load(c.p("aligned"), allow_pickle=True)
    x, pos = zc["coords"].astype(float), zc["positions"].astype(int)
    keep_pos = ~((pos >= 528) & (pos <= 550))
    # per-residue copy-to-copy distance (chains of one entry, superposed frame): 95th percentile per position
    pairs = np.argwhere(np.triu(same_entry, 1))
    dcopy = np.linalg.norm(x[pairs[:, 0]] - x[pairs[:, 1]], axis=-1)       # (pairs, positions), NaN if unobserved
    copy_res_p95 = np.array([np.nanpercentile(col, 95) if np.isfinite(col).sum() >= 20 else np.nan
                             for col in dcopy.T])
    comp_rows, assoc_rows, sep_rows, prof_rows, summary = [], [], [], [], {}
    for name, lab in labels.items():
        # composition
        for g, sub in out.groupby(name):
            sg = sub.space_group.value_counts(normalize=True)
            mu = sub.mutation_status.value_counts(normalize=True)
            comp_rows.append(dict(split=name, nolid_group=int(g), n_chains=len(sub), n_entries=sub.entry_id.nunique(),
                                  cl7=int((sub.group == "cl7").sum()), cl4=int((sub.group == "cl4").sum()),
                                  minor=int((sub.group == "minor").sum()),
                                  top_space_group=sg.index[0], top_space_group_share=float(sg.iloc[0]),
                                  top_mutation=mu.index[0], top_mutation_share=float(mu.iloc[0]),
                                  peptide_share=float(sub.peptide_bound.mean()),
                                  basic_amine_share=float(sub.ligand_basic_amine.mean()),
                                  h12_complete_share=float((sub.h12_label == "complete").mean()),
                                  chains=";".join(sub.chain_key) if len(sub) <= 10 else "",
                                  classes=";".join(f"{k}:{v}" for k, v in sub.action_class.value_counts().items())))
        # association
        ent = out.drop_duplicates(["entry_id", name])
        for cov in ("group", "space_group", "mutation_status", "peptide_bound", "ligand_basic_amine",
                    "ligand_chemotype", "action_class", "h12_label", "entry_id"):
            assoc_rows.append(dict(split=name, covariate=cov,
                                   ari_chain=adjusted_rand_score(out[cov].map(str), out[name]),
                                   cramers_v_entry=cramers_v(ent[cov].map(str), ent[name]) if cov != "entry_id"
                                   else np.nan))
        # separation of the two largest groups
        ia, ib = np.where(lab == 1)[0], np.where(lab == 2)[0]
        wa = d[np.ix_(ia, ia)][np.triu_indices(len(ia), 1)]
        wb = d[np.ix_(ib, ib)][np.triu_indices(len(ib), 1)]
        within = float(np.median(np.concatenate([wa, wb])))
        between = float(np.median(d[np.ix_(ia, ib)]))
        sep_rows.append(dict(split=name, n_group1=len(ia), n_group2=len(ib), within_median=within,
                             between_median=between, gap=between - within, copy_p95=copy_p95,
                             real=bool(between - within > copy_p95),
                             ari_vs_lid_groups=adjusted_rand_score(t.cluster_1p5, lab)))
        # per-residue distance between the two largest groups' mean positions (528-550 excluded)
        m1, m2 = np.nanmean(x[ia], 0), np.nanmean(x[ib], 0)
        n1, n2 = np.isfinite(x[ia][..., 0]).sum(0), np.isfinite(x[ib][..., 0]).sum(0)
        spread = np.sqrt(np.nanmean(np.concatenate([((x[ia] - m1) ** 2).sum(-1), ((x[ib] - m2) ** 2).sum(-1)]), 0))
        pr = pd.DataFrame({"split": name, "position": pos, "dist_A": np.linalg.norm(m1 - m2, axis=1),
                           "within_rms_A": spread, "copy_p95_A": copy_res_p95, "n_obs_group1": n1,
                           "n_obs_group2": n2, "in_lid_region": ~keep_pos})
        pr["beyond_ordinary"] = pr.dist_A > np.fmax(pr.within_rms_A, pr.copy_p95_A)
        prof_rows.append(pr)
        pk = pr[keep_pos & (pr.n_obs_group1 >= 5) & (pr.n_obs_group2 >= 5)]
        sq = pk.dist_A ** 2
        top = pk.sort_values("dist_A", ascending=False).head(8)
        sq = sq.fillna(0)
        summary[name] = dict(
            sizes=";".join(map(str, pd.Series(lab).value_counts().tolist())),
            separation=sep_rows[-1],
            top_residues=top[["position", "dist_A", "within_rms_A"]].round(3).to_dict("records"),
            n_residues_dist_gt_within=int((pk.dist_A > pk.within_rms_A).sum()),
            residues_beyond_ordinary=pk[pk.beyond_ordinary].position.tolist(),
            helix_share_sq_diff={h: float(sq[(pk.position >= a) & (pk.position <= b)].sum() / sq.sum())
                                 for h, (a, b) in HELICES.items()},
            loop_331_340_share=float(sq[(pk.position >= 331) & (pk.position <= 340)].sum() / sq.sum()),
            loop_457_472_share=float(sq[(pk.position >= 457) & (pk.position <= 472)].sum() / sq.sum()),
        )
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "no_lid_composition.csv", index=False)
    assoc = pd.DataFrame(assoc_rows)
    assoc.to_csv(OUT / "no_lid_association.csv", index=False)
    sep = pd.DataFrame(sep_rows)
    sep.to_csv(OUT / "no_lid_separation.csv", index=False)
    prof = pd.concat(prof_rows, ignore_index=True)
    prof.to_csv(OUT / "no_lid_residue_profile.csv", index=False)

    sc = scan.set_index("k")
    c.write_numbers(OUT, {
        "k_best": kbest, "silhouette_best": float(sc.loc[kbest, "silhouette"]), "sizes_best": sc.loc[kbest, "sizes"],
        "k_first_balanced": kbal, "silhouette_first_balanced": float(sc.loc[kbal, "silhouette"]) if kbal else None,
        "copy_p95_A": copy_p95, "median_copy_A": float(np.median(copy)),
        "ari_spectral_vs_balanced": adjusted_rand_score(labels["spectral_k2"], labels["average_first_balanced"])
        if "average_first_balanced" in labels else None,
        "splits": summary,
        "association_top": assoc.sort_values("ari_chain", ascending=False).groupby("split").head(3)
        .to_dict("records"),
    })
    print(scan.to_string())
    print(comp.drop(columns=["chains"]).to_string())
    print(assoc.sort_values(["split", "ari_chain"], ascending=[True, False]).to_string())
    print(sep.to_string())
    for k, v in summary.items():
        print(k, v["top_residues"], v["helix_share_sq_diff"], v["loop_331_340_share"], v["loop_457_472_share"],
              v["n_residues_dist_gt_within"], v["residues_beyond_ordinary"])
    print(comp[comp.chains != ""][["split", "nolid_group", "chains"]].to_string())


if __name__ == "__main__":
    main()
