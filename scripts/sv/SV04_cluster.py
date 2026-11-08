"""SV04 steps 4-5 (RUN_SV1028): cluster chains on helix midpoints, ARI against S08, landing of partial-H12 chains.

Inputs: results/sv/SV04/midpoints.csv (SV04_midpoints.py), results/sv/SV00/chain_table.csv (s08_cluster = S08 cut
1.5 A, the cut that gives 7 clusters; equal to results/entry_clusters.csv cut_1.5), results/entry_level_clusters.csv
(K07: one chain per entry by its stated rule), results/sv/SV03/geometry_per_chain.csv (SV03 H12 map coordinates).

Rules fixed before any output was looked at:
- Helices = the 8 SV02 consensus helices (H12 = 538-547); a chain x helix enters only if eligible (SV04 [DEFAULT]).
- Step 4 chain distance = sqrt(mean over the helices both chains have of |midpoint_i - midpoint_j|^2) (A).
  Axis angle reported separately: RMS over the same shared helices of the angle between the unit axes (deg); it is
  clustered the same way and gets its own ARI rows. A pair with no shared helix would be undefined: then the chain
  with fewest eligible helices is dropped until none is left (recorded).
- Average linkage (scipy). ARI (sklearn adjusted_rand_score) against S08:
    K = 2: cl7 and cl4 chains only, linkage recomputed on their sub-matrix, cut into 2 clusters, vs cl7/cl4;
    K = 7: all chains, cut into 7 clusters (maxclust), vs the S08 1.5 A labels (7 clusters).
  Chain level = all chains; entry level = the K07 chain per entry (results/entry_level_clusters.csv), sub-matrix
  sliced, linkage recomputed, same two cuts, vs the S08 labels of those chains.
- Step 5, partial-H12 chains = SV01a h12_status 'partial' (1-9 of the 10 C-alpha 538-547 observed). Where they land:
    landing_h12 (primary): the chain's estimated H12 midpoint is nearer to the cl7 or to the cl4 reference H12
      midpoint; reference = coordinate-wise median of the H12 midpoints of SV01a-complete-H12 chains of that S08
      group with an eligible H12. Needs an eligible H12 midpoint.
    landing_k2: the chain's SV04 cluster in the K = 2 clustering (cl7 and cl4 chains only), labelled by the majority
      S08 group of that cluster's SV01a-complete-H12 members.
    landing_k7: the chain's SV04 cluster at K = 7 (all chains), labelled the same way (cl7, cl4 or minor; 'none' if
      the cluster has no complete-H12 member).
  Agreement = landing label equals the chain's S08 group; for minor-cluster chains landing_h12 (binary cl7 / cl4)
  cannot agree by construction, so their agreement is n/a (their landing is still reported). Entry level: per entry with >= 1 landed partial chain,
  the majority landing and the majority S08 group of those chains (a tie is reported as a tie, not as agreement).
- Figure coordinates (SV03 H12 map axes): x = distance from the SV04 H12 midpoint to the SV04 H11 axis line (the
  H10/11 axis through its midpoint), y = angle between the SV04 H12 and H3 axes; a chain is placed only if H12, H3
  and H10/11 are all eligible. Filled points = SV01a-complete-H12 chains at their SV03 coordinates; hollow =
  partial-H12 chains at their SV04 estimated coordinates.
Outputs: results/sv/SV04/ari.csv, partial_h12_assignment.csv, h12_map_coords.csv, cluster_labels.csv,
  cluster_crosstab.csv, distance_summary.csv, numbers_cluster.json; data/sv/SV04/{midpoint_rms,axis_angle_rms,
  n_shared}.npy + matrix_index.csv.
Run: .venv/bin/python -B -I scripts/sv/SV04_cluster.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.helix_geometry import axis_angle, point_axis_distance  # noqa: E402

OUT = ROOT / "results/sv/SV04"
DATA = ROOT / "data/sv/SV04"
HELICES = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]


def matrices(mp, keys):
    n, h = len(keys), len(HELICES)
    M = np.full((n, h, 3), np.nan)
    U = np.full((n, h, 3), np.nan)
    E = np.zeros((n, h), bool)
    idx = {k: i for i, k in enumerate(keys)}
    for r in mp.itertuples():
        if r.helix not in HELICES:
            continue
        i, j = idx[r.chain_key], HELICES.index(r.helix)
        if r.eligible:
            M[i, j] = (r.mid_x, r.mid_y, r.mid_z)
            U[i, j] = (r.u_x, r.u_y, r.u_z)
            E[i, j] = True
    d2 = np.zeros((n, n))
    a2 = np.zeros((n, n))
    cnt = np.zeros((n, n), int)
    for j in range(h):
        e = E[:, j]
        both = np.outer(e, e)
        x = np.where(e[:, None], M[:, j], 0.0)
        u = np.where(e[:, None], U[:, j], 0.0)
        sq = ((x[:, None, :] - x[None, :, :]) ** 2).sum(-1)
        ang = np.degrees(np.arccos(np.clip(u @ u.T, -1.0, 1.0)))
        d2 += np.where(both, sq, 0.0)
        a2 += np.where(both, ang ** 2, 0.0)
        cnt += both
    with np.errstate(invalid="ignore", divide="ignore"):
        D = np.sqrt(d2 / cnt)
        A = np.sqrt(a2 / cnt)
    np.fill_diagonal(D, 0.0)
    np.fill_diagonal(A, 0.0)
    return M, U, E, D, A, cnt


def cut(dm, k):
    Z = linkage(squareform(dm, checks=False), method="average")
    return fcluster(Z, t=k, criterion="maxclust")


def majority(series):
    vc = series.value_counts()
    if not len(vc):
        return "none"
    if len(vc) > 1 and vc.iloc[0] == vc.iloc[1]:
        return "tie"
    return vc.index[0]


def main():
    mp = pd.read_csv(OUT / "midpoints.csv", keep_default_na=False, na_values=[""])
    mp["eligible"] = mp.eligible.astype(str).eq("True")
    ct = pd.read_csv(ROOT / "results/sv/SV00/chain_table.csv", keep_default_na=False, na_values=[""])
    ec = pd.read_csv(ROOT / "results/entry_clusters.csv").set_index("chain_key")
    assert (ec.loc[ct.chain_key, "cut_1.5"].to_numpy() == ct.s08_cluster.to_numpy()).all()
    el = pd.read_csv(ROOT / "results/entry_level_clusters.csv")
    keys = ct.chain_key.tolist()
    M, U, E, D, A, cnt = matrices(mp, keys)
    # pairs with no shared helix (rule: drop the chain with fewest eligible helices until none)
    keep = np.ones(len(keys), bool)
    dropped = []
    while True:
        sub = cnt[np.ix_(keep, keep)]
        if (sub > 0).all():
            break
        kk = np.where(keep)[0]
        bad = kk[(sub == 0).any(1)]
        worst = bad[np.argmin(E[bad].sum(1))]
        keep[worst] = False
        dropped.append(keys[worst])
    DATA.mkdir(parents=True, exist_ok=True)
    np.save(DATA / "midpoint_rms.npy", D)
    np.save(DATA / "axis_angle_rms.npy", A)
    np.save(DATA / "n_shared.npy", cnt)
    pd.DataFrame({"row": range(len(keys)), "chain_key": keys, "in_clustering": keep}).to_csv(
        DATA / "matrix_index.csv", index=False)

    group = ct.s08_cluster.map({7: "cl7", 4: "cl4"}).fillna("minor").to_numpy()
    s08 = ct.s08_cluster.to_numpy()
    oc = pd.read_csv(ROOT / "results/sv/SV01a/chain_outcomes.csv").set_index("chain_key")
    h12_status = oc.loc[keys, "h12_status"].to_numpy()
    rep = set(el.chain_key)
    is_rep = np.array([k in rep for k in keys])
    labels = pd.DataFrame({"chain_key": keys, "entry_id": ct.entry_id, "s08_cluster": s08, "group": group,
                           "h12_status": h12_status, "k07_entry_chain": is_rep, "n_eligible_helices": E.sum(1)})

    ari_rows, xt_rows = [], []
    for dname, dm in (("midpoint_rms", D), ("axis_angle_rms", A)):
        for level, base in (("chain", keep), ("entry", keep & is_rep)):
            for K, sel in ((2, base & np.isin(group, ["cl7", "cl4"])), (7, base)):
                ii = np.where(sel)[0]
                lab = cut(dm[np.ix_(ii, ii)], K)
                ref = group[ii] if K == 2 else s08[ii]
                col = f"{dname}_{level}_K{K}"
                labels[col] = pd.Series(lab, index=ii).reindex(range(len(keys))).astype("Int64")
                sizes = pd.Series(lab).value_counts().sort_values(ascending=False)
                ari_rows.append({"distance": dname, "level": level, "K": K,
                                 "chains": "cl7 + cl4 chains" if K == 2 else "all chains",
                                 "reference": "S08 group cl7 / cl4" if K == 2 else "S08 cut 1.5 A (7 clusters)",
                                 "n": int(len(ii)), "n_clusters_sv04": int(len(sizes)),
                                 "n_clusters_reference": int(len(set(ref))),
                                 "sv04_cluster_sizes": ";".join(map(str, sizes.to_list())),
                                 "ari": adjusted_rand_score(ref, lab)})
                xt = pd.crosstab(pd.Series(lab, name="sv04_cluster"),
                                 pd.Series(ref, name="s08"))
                for c_ in xt.index:
                    for r_ in xt.columns:
                        if xt.loc[c_, r_]:
                            xt_rows.append({"distance": dname, "level": level, "K": K, "sv04_cluster": int(c_),
                                            "s08": str(r_), "n": int(xt.loc[c_, r_])})
    ari = pd.DataFrame(ari_rows)
    ari.to_csv(OUT / "ari.csv", index=False)
    pd.DataFrame(xt_rows).to_csv(OUT / "cluster_crosstab.csv", index=False)
    labels.to_csv(OUT / "cluster_labels.csv", index=False)

    # distance summary per helix (chain pairs): midpoint distance and axis angle within / between cl7 and cl4
    srows = []
    g7, g4 = np.where(group == "cl7")[0], np.where(group == "cl4")[0]
    for j, hname in enumerate(HELICES + ["RMS over shared helices"]):
        for what in ("midpoint_distance_A", "axis_angle_deg"):
            rec = {"helix": hname, "measure": what}
            for pname, (a_, b_) in (("within_cl7", (g7, g7)), ("within_cl4", (g4, g4)), ("cl7_vs_cl4", (g7, g4))):
                if hname.startswith("RMS"):
                    mat = (D if what.startswith("midpoint") else A)[np.ix_(a_, b_)]
                    ok = np.ones_like(mat, bool)
                else:
                    ea, eb = E[a_, j], E[b_, j]
                    if what.startswith("midpoint"):
                        mat = np.linalg.norm(M[a_, j][:, None] - M[b_, j][None], axis=-1)
                    else:
                        mat = np.degrees(np.arccos(np.clip(np.einsum("ik,jk->ij", np.nan_to_num(U[a_, j]),
                                                                     np.nan_to_num(U[b_, j])), -1, 1)))
                    ok = np.outer(ea, eb)
                if a_ is b_:
                    ok = ok & np.triu(np.ones_like(ok, bool), 1)
                v = mat[ok]
                rec[f"{pname}_n_pairs"] = int(v.size)
                rec[f"{pname}_median"] = float(np.median(v)) if v.size else np.nan
                rec[f"{pname}_q90"] = float(np.quantile(v, 0.9)) if v.size else np.nan
            srows.append(rec)
    pd.DataFrame(srows).to_csv(OUT / "distance_summary.csv", index=False)

    # step 5: partial-H12 chains
    h12 = HELICES.index("H12")
    h3, h11 = HELICES.index("H3"), HELICES.index("H10/11")
    complete = h12_status == "complete"
    refs = {g: np.median(M[complete & (group == g) & E[:, h12], h12], axis=0) for g in ("cl7", "cl4")}
    n_ref = {g: int((complete & (group == g) & E[:, h12]).sum()) for g in ("cl7", "cl4")}

    def clabel(col):
        lab = labels[col]
        out = {}
        for c_ in lab.dropna().unique():
            m_ = (lab == c_).fillna(False).to_numpy() & complete
            out[int(c_)] = majority(pd.Series(group[m_]))
        return out

    lab_k2 = clabel("midpoint_rms_chain_K2")
    lab_k7 = clabel("midpoint_rms_chain_K7")
    gm = pd.read_csv(ROOT / "results/sv/SV03/geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    gm = gm.set_index("chain_key").loc[keys]
    rows = []
    for i, k in enumerate(keys):
        placed = E[i, h12] and E[i, h3] and E[i, h11]
        x4 = point_axis_distance(M[i, h12], M[i, h11], U[i, h11]) if placed else np.nan
        y4 = axis_angle(U[i, h12], U[i, h3]) if placed else np.nan
        d7 = float(np.linalg.norm(M[i, h12] - refs["cl7"])) if E[i, h12] else np.nan
        d4 = float(np.linalg.norm(M[i, h12] - refs["cl4"])) if E[i, h12] else np.nan
        l_h12 = ("cl7" if d7 < d4 else "cl4") if E[i, h12] else "not placed"
        c2 = labels.loc[i, "midpoint_rms_chain_K2"]
        c7 = labels.loc[i, "midpoint_rms_chain_K7"]
        l2 = lab_k2.get(int(c2), "none") if pd.notna(c2) else "n/a (minor cluster)"
        l7 = lab_k7.get(int(c7), "none") if pd.notna(c7) else "not clustered"
        rows.append({"chain_key": k, "entry_id": ct.entry_id[i], "s08_cluster": int(s08[i]), "group": group[i],
                     "h12_status": h12_status[i], "n_obs_window_538_547": int(oc.loc[k, "n_obs_window"]),
                     "ligand_chemotype": ct.ligand_chemotype[i], "h12_eligible": bool(E[i, h12]),
                     "placed_on_map": bool(placed),
                     "mid_H12_x": M[i, h12, 0], "mid_H12_y": M[i, h12, 1], "mid_H12_z": M[i, h12, 2],
                     "d_ref_cl7_A": d7, "d_ref_cl4_A": d4,
                     "landing_h12": l_h12 if E[i, h12] else "not placed",
                     "agree_h12": (l_h12 == group[i]) if E[i, h12] and group[i] != "minor" else pd.NA,
                     "sv04_cluster_K2": c2, "landing_k2": l2,
                     "agree_k2": (l2 == group[i]) if pd.notna(c2) else pd.NA,
                     "sv04_cluster_K7": c7, "landing_k7": l7, "agree_k7": (l7 == group[i]) if pd.notna(c7) else pd.NA,
                     "map_x_sv04": x4, "map_y_sv04": y4,
                     "map_x_sv03": gm.loc[k, "dist_H12c_H11axis_A"], "map_y_sv03": gm.loc[k, "angle_H12_H3_deg"]})
    allc = pd.DataFrame(rows)
    allc["map_layer"] = np.where(allc.h12_status.eq("complete") & allc.map_x_sv03.notna() & allc.map_y_sv03.notna(),
                                 "filled (SV03 position)",
                                 np.where(allc.h12_status.eq("partial") & allc.placed_on_map,
                                          "hollow (SV04 estimate)", "not drawn"))
    allc.to_csv(OUT / "h12_map_coords.csv", index=False)
    part = allc[allc.h12_status == "partial"].drop(columns=["map_layer"]).copy()
    part.to_csv(OUT / "partial_h12_assignment.csv", index=False)

    # numbers
    def agree_tab(df, col):
        out = {}
        for g, x in df.groupby("group"):
            vc = x[col].value_counts().to_dict()
            out[g] = {"n": int(len(x)), **{str(a): int(b) for a, b in vc.items()}}
        return out

    ent = []
    for e_, x in part.groupby("entry_id"):
        for col in ("landing_h12", "landing_k2", "landing_k7"):
            y = x[~x[col].isin(["not placed", "n/a (minor cluster)", "not clustered"])]
            if not len(y):
                continue
            ml, mg = majority(y[col]), majority(y["group"])
            ent.append({"entry_id": e_, "landing": col, "n_chains": int(len(y)), "majority_landing": ml,
                        "majority_group": mg,
                        "agree": (ml == mg) if "tie" not in (ml, mg) and not (mg == "minor" and col != "landing_k7")
                        else pd.NA})
    ent = pd.DataFrame(ent)
    ent.to_csv(OUT / "partial_h12_entry_level.csv", index=False)
    comp = allc[allc.map_layer.str.startswith("filled") & allc.map_x_sv04.notna()]
    num = {"dropped_no_shared_helix": dropped, "n_clustered": int(keep.sum()),
           "n_shared_helices": {int(k): int(v) for k, v in
                                pd.Series(cnt[np.triu_indices(len(keys), 1)]).value_counts().sort_index().items()},
           "reference_h12_midpoints": {g: list(map(float, refs[g])) for g in refs}, "n_reference": n_ref,
           "reference_distance_cl7_cl4_A": float(np.linalg.norm(refs["cl7"] - refs["cl4"])),
           "cluster_label_k2": {str(k): v for k, v in lab_k2.items()},
           "cluster_label_k7": {str(k): v for k, v in lab_k7.items()},
           "partial": {"n": int(len(part)), "by_group": part.group.value_counts().to_dict(),
                       "h12_eligible": int(part.h12_eligible.sum()), "placed_on_map": int(part.placed_on_map.sum()),
                       "placed_by_group": part[part.placed_on_map].group.value_counts().to_dict(),
                       "landing_h12": agree_tab(part, "landing_h12"), "landing_k2": agree_tab(part, "landing_k2"),
                       "landing_k7": agree_tab(part, "landing_k7"),
                       "agree_h12": {g: [int(x.agree_h12.fillna(False).astype(bool).sum()), int(x.agree_h12.notna().sum())]
                                     for g, x in part.groupby("group")},
                       "agree_k2": {g: [int(x.agree_k2.fillna(False).astype(bool).sum()), int(x.agree_k2.notna().sum())]
                                    for g, x in part.groupby("group")},
                       "agree_k7": {g: [int(x.agree_k7.fillna(False).astype(bool).sum()), int(x.agree_k7.notna().sum())]
                                    for g, x in part.groupby("group")},
                       "agree_by_h12_eligible": {
                           f"{col}|{g}|h12_eligible={el_}": [int(x[col].fillna(False).astype(bool).sum()),
                                                             int(x[col].notna().sum())]
                           for col in ("agree_h12", "agree_k2", "agree_k7")
                           for (g, el_), x in part.groupby(["group", "h12_eligible"])},
                       "margin_d4_minus_d7_A": {g: {"median": float(x.eval("d_ref_cl4_A - d_ref_cl7_A").median()),
                                                    "min": float(x.eval("d_ref_cl4_A - d_ref_cl7_A").min()),
                                                    "max": float(x.eval("d_ref_cl4_A - d_ref_cl7_A").max())}
                                                for g, x in part[part.h12_eligible].groupby("group")},
                       "sv04_minus_sv03_map": {"n": int((part.map_x_sv04.notna() & part.map_x_sv03.notna()).sum()),
                                               "x_median": float((part.map_x_sv04 - part.map_x_sv03).median()),
                                               "x_median_abs": float((part.map_x_sv04 - part.map_x_sv03).abs().median()),
                                               "x_max_abs": float((part.map_x_sv04 - part.map_x_sv03).abs().max()),
                                               "y_max_abs": float((part.map_y_sv04 - part.map_y_sv03).abs().max()),
                                               "placed_without_sv03_point": int((part.placed_on_map &
                                                                                 part.map_x_sv03.isna()).sum())}},
           "entry_level": {col: {"n_entries": int((ent.landing == col).sum()),
                                 "agree": int(ent[(ent.landing == col)].agree.fillna(False).astype(bool).sum()),
                                 "disagree": int((ent[(ent.landing == col)].agree == False).sum()),  # noqa: E712
                                 "tie": int(((ent.landing == col) & ent.agree.isna() &
                                             (ent.majority_group == "tie")).sum() +
                                            ((ent.landing == col) & (ent.majority_landing == "tie") &
                                             (ent.majority_group != "tie")).sum()),
                                 "n_a_minor": int(((ent.landing == col) & ent.agree.isna() &
                                                   (ent.majority_group == "minor")).sum()),
                                 "by_group": {g: [int(x.agree.fillna(False).astype(bool).sum()), int(len(x))]
                                              for g, x in ent[ent.landing == col].groupby("majority_group")}}
                           for col in ("landing_h12", "landing_k2", "landing_k7")},
           "complete_sv04_minus_sv03_map": {"n": int(len(comp)),
                                            "x_median": float((comp.map_x_sv04 - comp.map_x_sv03).median()),
                                            "x_median_abs": float((comp.map_x_sv04 - comp.map_x_sv03).abs().median()),
                                            "x_max_abs": float((comp.map_x_sv04 - comp.map_x_sv03).abs().max()),
                                            "y_max_abs": float((comp.map_y_sv04 - comp.map_y_sv03).abs().max())},
           "sv04_minus_sv03_map_x_by_status_group": {
               f"{st}|{g}": {"n": int(len(x)), "median": float((x.map_x_sv04 - x.map_x_sv03).median()),
                             "min": float((x.map_x_sv04 - x.map_x_sv03).min()),
                             "max": float((x.map_x_sv04 - x.map_x_sv03).max())}
               for (st, g), x in allc[allc.map_x_sv04.notna() & allc.map_x_sv03.notna()].groupby(["h12_status", "group"])},
           "k2_chain_mismatch": (lambda x: {
               "n": int(len(x)), "n_without_usable_h12": int((x.n_eligible_helices < 8).sum()),
               "by_h12_status": x.h12_status.value_counts().to_dict(),
               "by_group": x.group.value_counts().to_dict()})(
               labels[labels.midpoint_rms_chain_K2.notna() &
                      (labels.midpoint_rms_chain_K2.map(lambda c: lab_k2.get(int(c)) if pd.notna(c) else None)
                       != labels.group)]),
           "n_without_usable_h12_by_group": labels[labels.n_eligible_helices < 8].group.value_counts().to_dict(),
           "map_layers": allc.map_layer.value_counts().to_dict(),
           "map_layers_by_group": {L: x.group.value_counts().to_dict() for L, x in allc.groupby("map_layer")}}
    (OUT / "numbers_cluster.json").write_text(json.dumps(num, indent=1, default=lambda v: v.item()
                                                         if hasattr(v, "item") else str(v)))
    pd.set_option("display.width", 220)
    print(ari.to_string())
    print(json.dumps({k: num[k] for k in ("n_clustered", "n_shared_helices", "n_reference",
                                         "reference_distance_cl7_cl4_A", "cluster_label_k2", "cluster_label_k7",
                                         "partial", "entry_level", "complete_sv04_minus_sv03_map", "map_layers")},
                     indent=1, default=str))


if __name__ == "__main__":
    main()
