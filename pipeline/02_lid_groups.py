"""Stage 02 - superposition, RMSD and the groups defined by the lid (H12) position.

1. Superposition: every chain's C-alpha atoms (deposited coordinates, data/coordinates/backbone.npz) are fitted onto
   the reference chain 7NFB|B|1 by Kabsch on the core-v3 positions both chains share (123 positions; H3, H4/5, H7,
   H8, H10/11), one fit per chain.
2. RMSD: C-alpha RMSD over the positions observed in both chains of a pair, all 246 LBD positions (305-550); the same
   without 528-550 (H11-H12 loop, H12, C-terminal tail).
3. Groups: average linkage, cuts at 1.0 / 1.5 / 2.0 / 2.5 A (1.5 A primary). cl7 and cl4 are the two main groups.
4. Where the groups differ: per-residue distance between the cl7 and cl4 mean C-alpha positions.
Checks: the labels must equal the frozen October groups (data/structures/frozen_groups.csv).
Outputs (results/02_lid_groups/): aligned_ca.npz, rmsd.npy, rmsd_without_528_550.npy, chains.csv, groups.csv,
displacement_profile.csv, numbers.json.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("02_lid_groups")
CUTS = (1.0, 1.5, 2.0, 2.5)


def kabsch(mobile, ref):
    pc, qc = mobile.mean(0), ref.mean(0)
    u, _, vt = np.linalg.svd((mobile - pc).T @ (ref - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    return rot, pc, qc


def main():
    z = np.load(c.p("backbone"), allow_pickle=True)
    keys = [str(k) for k in z["chain_keys"]]
    pos = z["positions"].astype(int)
    ca = z["coords"][:, :, 1, :].astype(float)
    core = [q for a, b in json.loads(c.p("core").read_text())["core_ranges"] for q in range(a, b + 1)]
    ref_key = json.loads(c.p("reference").read_text())["chain_key"]
    ref = ca[keys.index(ref_key)]
    in_core = np.isin(pos, core)
    aligned = np.full_like(ca, np.nan)
    fit = []
    for i, k in enumerate(keys):
        m = in_core & np.isfinite(ca[i, :, 0]) & np.isfinite(ref[:, 0])
        rot, pc, qc = kabsch(ca[i][m], ref[m])
        aligned[i] = (ca[i] - pc) @ rot.T + qc
        fit.append(dict(row=i, chain_key=k, entry_id=k.split("|")[0], n_core_ca=int(m.sum()),
                        core_rmsd_A=float(np.sqrt(((aligned[i][m] - ref[m]) ** 2).sum(1).mean()))))
    idx = pd.DataFrame(fit)
    idx.to_csv(OUT / "chains.csv", index=False)
    np.savez_compressed(OUT / "aligned_ca.npz", chain_keys=np.array(keys), positions=pos, coords=aligned)

    d = c.masked_rmsd(aligned)
    lid = (pos >= 528) & (pos <= 550)
    d_nolid = c.masked_rmsd(aligned, ~lid)
    np.save(OUT / "rmsd.npy", d)
    np.save(OUT / "rmsd_without_528_550.npy", d_nolid)

    zl = linkage(squareform(d, checks=False), "average")
    gr = idx[["chain_key", "entry_id"]].copy()
    for cut in CUTS:
        gr[f"cut_{cut}"] = fcluster(zl, t=cut, criterion="distance")
    gr.to_csv(OUT / "groups.csv", index=False)
    frozen = pd.read_csv(c.p("frozen_groups")).set_index("chain_key").loc[gr.chain_key]
    same = {cut: bool((frozen[f"cut_{cut}"].to_numpy() == gr[f"cut_{cut}"].to_numpy()).all()) for cut in CUTS}
    ari = {cut: adjusted_rand_score(frozen[f"cut_{cut}"], gr[f"cut_{cut}"]) for cut in CUTS}
    assert all(same.values()), f"groups differ from the frozen October groups: {ari}"

    grp = gr["cut_1.5"].map(c.group_of).to_numpy()
    m7, m4 = grp == "cl7", grp == "cl4"
    iu = np.triu_indices(len(d), 1)

    def within(m):
        sub = d[np.ix_(m, m)]
        return sub[np.triu_indices(m.sum(), 1)]

    same_entry = idx.entry_id.to_numpy()[:, None] == idx.entry_id.to_numpy()[None, :]
    same_group = grp[:, None] == grp[None, :]
    copies = d[np.triu(same_entry & same_group, 1)]

    mean7, mean4 = np.nanmean(aligned[m7], 0), np.nanmean(aligned[m4], 0)
    prof = pd.DataFrame({"position": pos, "distance_cl7_cl4_A": np.linalg.norm(mean7 - mean4, axis=1),
                         "n_cl7": np.isfinite(aligned[m7][..., 0]).sum(0), "n_cl4": np.isfinite(aligned[m4][..., 0]).sum(0)})
    prof.to_csv(OUT / "displacement_profile.csv", index=False)
    sq = prof.distance_cl7_cl4_A ** 2
    n = {
        "reference": ref_key, "core_positions": len(core), "chains": len(keys),
        "core_fit_rmsd_median_A": float(idx.core_rmsd_A.median()),
        "groups_equal_frozen": same, "n_groups_by_cut": {cut: int(gr[f"cut_{cut}"].nunique()) for cut in CUTS},
        "group_sizes_1p5": pd.Series(grp).value_counts().to_dict(),
        "entries_by_group": pd.Series(grp).groupby(idx.entry_id).first().value_counts().to_dict(),
        "median_rmsd_within_cl7": float(np.median(within(m7))), "median_rmsd_within_cl4": float(np.median(within(m4))),
        "mean_rmsd_within_cl7": float(np.mean(within(m7))), "mean_rmsd_within_cl4": float(np.mean(within(m4))),
        "median_rmsd_cl7_cl4": float(np.median(d[np.ix_(m7, m4)])),
        "median_rmsd_cl7_cl4_without_528_550": float(np.median(d_nolid[np.ix_(m7, m4)])),
        "median_rmsd_all_pairs": float(np.median(d[iu])),
        "share_sq_difference_528_550": float(sq[(prof.position >= 528) & (prof.position <= 550)].sum() / sq.sum()),
        "max_distance_outside_528_550_A": float(prof[(prof.position < 528)].distance_cl7_cl4_A.max()),
        "copy_rmsd_median": float(np.median(copies)), "copy_rmsd_p95": float(np.percentile(copies, 95)),
    }
    ct = c.read("chain_table").set_index("chain_key").loc[idx.chain_key]
    for g in ("cl7", "cl4"):
        sub = ct[grp == g]
        n[f"peptide_bound_{g}"] = int(sub.peptide_bound.astype(str).eq("True").sum())
        n[f"constructs_{g}"] = sub.mutation_status.value_counts().to_dict()
        n[f"h12_complete_{g}"] = int(sub.h12_label.eq("complete").sum())
    c.write_numbers(OUT, n)
    print(json.dumps({k: v for k, v in n.items()}, indent=1, default=str))


if __name__ == "__main__":
    main()
