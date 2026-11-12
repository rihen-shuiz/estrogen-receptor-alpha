"""Stage 07 - superposition-free description of every chain: helix ranges from phi/psi, helix lengths, a centre of mass (COM)
per helix and per chain, and the pairwise distances between them. No chain is fitted onto another.

Declared before running:
  * Helical residue: phi in [-100, -30] deg and psi in [-80, -10] deg (backbone dihedrals computed here from the deposited N, CA, C atoms).
    A helix = a run of >= 4 consecutive helical residues. For each consensus helix (data/coordinates/helix_boundaries.csv) the run with
    the largest overlap inside a search window (consensus range +- 6 residues, H12 window 532-550) gives start,
    end and length; 0 when the window is observed but holds no overlapping run.
  * Helix COM: mass-weighted COM of the backbone atoms (N, CA, C, O) of the observed residues of the consensus range,
    in the deposited crystal frame; defined when >= 70 % of the range (and >= 5 residues) is observed.
  * Chain COM: mass-weighted backbone COM over the positions observed in >= 99 % of the 788 chains.
  * Features (all invariant to rotation and translation): 28 helix-helix COM distances, 8 helix-to-chain-COM
    distances, 8 phi/psi helix lengths.
  * Separation cl7 vs cl4 per feature: difference of medians, AUC, and "real" only if |difference| exceeds both
    the pooled within-group SD and the 95th percentile of copy-to-copy differences (chains of one entry, same group).
    Ranking: random forest (500 trees, seed 42), permutation importance on entry-grouped 5-fold splits.
  * Clustering: Ward linkage on z-scored features (chains with every feature), K = 2 and K by silhouette (2..10);
    ARI with the RMSD (1.5 A) groups. Repeated without every H12 feature.
  * Protein parts: each feature belongs to its helices; a part's score = the number of real features that involve it
    and the summed permutation importance (pairwise features split half-half).
Outputs: results/07_geometry_free/geometry_free_per_chain.csv, geometry_free_features.csv, geometry_free_parts.csv,
geometry_free_clusters.csv, geometry_free_kscan.csv, geometry_free_length_check.csv, numbers_geometry_free.json.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import adjusted_rand_score, roc_auc_score, silhouette_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("07_geometry_free")

PHI = (-100.0, -30.0)
PSI = (-80.0, -10.0)
MIN_RUN = 4
MASS = np.array([14.007, 12.011, 12.011, 15.999])  # N, CA, C, O
HNAMES = ["H1", "H3", "H4/5", "H7", "H8", "H9", "H10/11", "H12"]


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    return np.degrees(np.arctan2((np.cross(b1, v) * w).sum(-1), (v * w).sum(-1)))


def phipsi_table() -> pd.DataFrame:
    """phi / psi of every residue from the deposited backbone (NaN where a neighbour atom is missing)."""
    z = np.load(c.p("backbone"), allow_pickle=True)
    x = z["coords"].astype(float)
    n_, ca, c_ = x[:, :, 0], x[:, :, 1], x[:, :, 2]
    phi = np.full(n_.shape[:2], np.nan)
    psi = np.full(n_.shape[:2], np.nan)
    with np.errstate(invalid="ignore"):
        phi[:, 1:] = dihedral(c_[:, :-1], n_[:, 1:], ca[:, 1:], c_[:, 1:])
        psi[:, :-1] = dihedral(n_[:, :-1], ca[:, :-1], c_[:, :-1], n_[:, 1:])
    keys = [str(k) for k in z["chain_keys"]]
    pos = z["positions"].astype(int)
    obs = np.isfinite(ca[..., 0])
    i, j = np.where(obs)
    return pd.DataFrame({"chain_key": np.array(keys)[i], "position": pos[j], "phi": phi[i, j], "psi": psi[i, j]})


def helix_runs(ok: np.ndarray, pos: np.ndarray):
    """Runs of consecutive True (positions must also be consecutive)."""
    runs, start = [], None
    for i in range(len(pos)):
        if ok[i] and (start is None or pos[i] == pos[i - 1] + 1):
            start = pos[i] if start is None else start
        else:
            if start is not None:
                runs.append((start, pos[i - 1]))
            start = pos[i] if ok[i] else None
    if start is not None:
        runs.append((start, pos[-1]))
    return [(a, b) for a, b in runs if b - a + 1 >= MIN_RUN]


def main():
    t = c.chains()
    hb = pd.read_csv(c.p("helix_bounds")).set_index("helix_name")
    rng = {h: (int(hb.loc[h, "global_start"]), int(hb.loc[h, "global_end"])) for h in HNAMES}
    win = {h: (a - 6, b + 6) for h, (a, b) in rng.items()}
    win["H12"] = (532, 550)

    # ---------------- phi/psi helix ranges
    pp = phipsi_table()
    pp["helical"] = pp.phi.between(*PHI) & pp.psi.between(*PSI)
    len_rows = []
    for key, d in pp.groupby("chain_key", sort=False):
        d = d.sort_values("position")
        row = {"chain_key": key}
        for h in HNAMES:
            a, b = win[h]
            w = d[(d.position >= a) & (d.position <= b)]
            if w.phi.notna().sum() + w.psi.notna().sum() == 0:
                row[f"len_{h}"], row[f"start_{h}"], row[f"end_{h}"] = np.nan, np.nan, np.nan
                continue
            runs = helix_runs(w.helical.to_numpy(), w.position.to_numpy())
            ca, cb = rng[h]
            best, ov_best = None, 0
            for r0, r1 in runs:
                ov = max(0, min(r1, cb) - max(r0, ca) + 1)
                if ov > ov_best:
                    best, ov_best = (r0, r1), ov
            if best is None:
                row[f"len_{h}"], row[f"start_{h}"], row[f"end_{h}"] = 0, np.nan, np.nan
            else:
                row[f"len_{h}"], row[f"start_{h}"], row[f"end_{h}"] = best[1] - best[0] + 1, best[0], best[1]
        len_rows.append(row)
    lens = pd.DataFrame(len_rows)

    # ---------------- COMs in the deposited frame
    z = np.load(c.p("backbone"), allow_pickle=True)
    keys = [str(k) for k in z["chain_keys"]]
    assert keys == t.chain_key.tolist()
    xyz, pos = z["coords"].astype(float), z["positions"].astype(int)       # (788, 246, 4, 3)
    obs_res = np.isfinite(xyz[..., 1, 0])                                   # CA observed
    fixed = obs_res.mean(0) >= 0.99

    def com(i, mask):
        a = xyz[i][mask]                                                    # (res, 4, 3)
        m = np.isfinite(a[..., 0])
        w = np.where(m, MASS[None, :], 0.0)
        return (np.nan_to_num(a) * w[..., None]).sum((0, 1)) / w.sum()

    feats = []
    for i, key in enumerate(keys):
        f = {"chain_key": key, "n_fixed_observed": int((fixed & obs_res[i]).sum())}
        cc = com(i, fixed & obs_res[i])
        hc = {}
        for h in HNAMES:
            a, b = rng[h]
            m = (pos >= a) & (pos <= b)
            n_obs = int((m & obs_res[i]).sum())
            f[f"nobs_{h}"] = n_obs
            if n_obs >= 5 and n_obs >= 0.7 * m.sum():
                hc[h] = com(i, m & obs_res[i])
            else:
                hc[h] = None
        for h in HNAMES:
            f[f"d_{h}_COM"] = float(np.linalg.norm(hc[h] - cc)) if hc[h] is not None else np.nan
        for h1, h2 in itertools.combinations(HNAMES, 2):
            f[f"d_{h1}_{h2}"] = (float(np.linalg.norm(hc[h1] - hc[h2])) if hc[h1] is not None and hc[h2] is not None
                                 else np.nan)
        feats.append(f)
    fe = pd.DataFrame(feats).merge(lens, on="chain_key", how="left", validate="one_to_one")
    per = t[["chain_key", "entry_id", "group", "cluster_1p5", "mobility_h12_z", "h12_label"]].merge(
        fe, on="chain_key", validate="one_to_one")
    per.to_csv(OUT / "geometry_free_per_chain.csv", index=False)

    # ---------------- helix lengths vs the hydrogen-bond (pydssp) lengths of the 28 Oct helix geometry
    g = pd.read_csv(c.p("geometry"))
    chk = []
    for h, col in zip(HNAMES, ["H1", "H3", "H45", "H7", "H8", "H9", "H1011", "H12"]):
        m = per[["chain_key", f"len_{h}"]].rename(columns={f"len_{h}": "phipsi"}).merge(
            g[["chain_key", f"len_{col}"]].rename(columns={f"len_{col}": "hbond"}), on="chain_key").dropna()
        chk.append(dict(helix=h, n=len(m), spearman=spearmanr(m.phipsi, m.hbond)[0],
                        median_phipsi=float(m.phipsi.median()), median_hbond=float(m.hbond.median()),
                        median_abs_diff=float((m.phipsi - m.hbond).abs().median())))
    pd.DataFrame(chk).to_csv(OUT / "geometry_free_length_check.csv", index=False)

    # ---------------- per-feature separation cl7 vs cl4 and ordinary variation
    fcols = [x for x in per.columns if x.startswith("d_") or x.startswith("len_")]
    main = per[per.group.isin(["cl7", "cl4"])].copy()
    rows = []
    for fcol in fcols:
        a, b = main[main.group == "cl7"][fcol].dropna(), main[main.group == "cl4"][fcol].dropna()
        if len(a) < 10 or len(b) < 10:
            continue
        cd = c.copy_diffs(main, fcol)
        pooled_sd = float(np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2)))
        copy_p95 = float(np.percentile(cd, 95)) if len(cd) else np.nan
        thr = float(np.nanmax([pooled_sd, copy_p95]))
        diff = float(b.median() - a.median())
        y = np.r_[np.zeros(len(a)), np.ones(len(b))]
        auc = roc_auc_score(y, np.r_[a, b])
        parts = fcol[2:].split("_") if fcol.startswith("d_") else [fcol[4:]]
        rows.append(dict(feature=fcol, kind="distance" if fcol.startswith("d_") else "length",
                         part_1=parts[0], part_2=parts[1] if len(parts) > 1 else "",
                         n_cl7=len(a), n_cl4=len(b), median_cl7=float(a.median()), median_cl4=float(b.median()),
                         diff_cl4_minus_cl7=diff, sd_cl7=float(a.std()), sd_cl4=float(b.std()), pooled_sd=pooled_sd,
                         copy_p95=copy_p95, threshold=thr, auc=auc, separation=abs(auc - 0.5) * 2,
                         p=mannwhitneyu(a, b).pvalue, real=bool(abs(diff) > thr)))
    ft = pd.DataFrame(rows)
    ft["q"] = c.bh(ft.p)

    # random forest with entry-grouped folds (chains with every feature)
    X = main.dropna(subset=fcols)
    yv = (X.group == "cl4").astype(int).to_numpy()
    imp, pimp = np.zeros(len(fcols)), np.zeros(len(fcols))
    accs = []
    for tr, te in GroupKFold(5).split(X, yv, X.entry_id):
        rf = RandomForestClassifier(500, random_state=c.SEED, n_jobs=-1).fit(X.iloc[tr][fcols], yv[tr])
        accs.append(rf.score(X.iloc[te][fcols], yv[te]))
        imp += rf.feature_importances_ / 5          # impurity importance (permutation is ~0: redundant H12 features)
        pi = permutation_importance(rf, X.iloc[te][fcols], yv[te], n_repeats=5, random_state=c.SEED, n_jobs=-1)
        pimp += pi.importances_mean / 5
    ft = ft.merge(pd.DataFrame({"feature": fcols, "rf_importance": imp, "perm_importance": pimp}), on="feature",
                  how="left")
    ft = ft.sort_values(["real", "separation"], ascending=[False, False])
    ft.to_csv(OUT / "geometry_free_features.csv", index=False)

    # ---------------- protein parts
    parts = []
    for h in HNAMES + ["COM"]:
        inv = ft[(ft.part_1 == h) | (ft.part_2 == h)]
        w = np.where(inv.part_2 == "", 1.0, 0.5)
        parts.append(dict(part=h, n_features=len(inv), n_real=int(inv.real.sum()),
                          max_abs_diff_A=float(inv.diff_cl4_minus_cl7.abs().max()) if len(inv) else np.nan,
                          importance=float((inv.rf_importance * w).sum())))
    pt = pd.DataFrame(parts).sort_values(["n_real", "importance"], ascending=False)
    pt.to_csv(OUT / "geometry_free_parts.csv", index=False)

    # ---------------- clustering on the superposition-free features
    def cluster(cols, tag):
        Y = per.dropna(subset=cols)
        Z = (Y[cols] - Y[cols].mean()) / Y[cols].std(ddof=0)
        Z = Z.loc[:, Z.std() > 0]
        L = linkage(Z.to_numpy(), "ward")
        scan = []
        for k in range(2, 11):
            lab = fcluster(L, k, "maxclust")
            scan.append(dict(features=tag, k=k, silhouette=silhouette_score(Z, lab),
                             ari_vs_rmsd_groups=adjusted_rand_score(Y.cluster_1p5, lab),
                             ari_vs_cl7_cl4=adjusted_rand_score(Y.group.where(Y.group != "minor", "minor"), lab),
                             sizes=";".join(map(str, pd.Series(lab).value_counts().tolist()))))
        scan = pd.DataFrame(scan)
        kb = int(scan.loc[scan.silhouette.idxmax(), "k"])
        lab2, labb = fcluster(L, 2, "maxclust"), fcluster(L, kb, "maxclust")
        lab_df = pd.DataFrame({"chain_key": Y.chain_key, f"{tag}_k2": lab2, f"{tag}_best": labb})
        return scan, kb, lab_df, len(Y)

    h12_cols = [x for x in fcols if "H12" in x]
    s_all, kb_all, lab_all, n_all = cluster(fcols, "all_features")
    nolid_cols = [x for x in fcols if x not in h12_cols]
    s_nl, kb_nl, lab_nl, n_nl = cluster(nolid_cols, "without_H12")
    scan = pd.concat([s_all, s_nl], ignore_index=True)
    scan.to_csv(OUT / "geometry_free_kscan.csv", index=False)
    labs = per[["chain_key", "group", "cluster_1p5"]].merge(lab_all, on="chain_key", how="left").merge(
        lab_nl, on="chain_key", how="left")
    labs.to_csv(OUT / "geometry_free_clusters.csv", index=False)

    # what the K = 2 split without H12 features is: composition and the features that separate it
    nl = per.merge(lab_nl, on="chain_key").merge(t[["chain_key", "space_group", "mutation_status", "peptide_bound",
                                                   "ligand_basic_amine"]], on="chain_key")
    small = nl.without_H12_k2.value_counts().index[-1]
    nl["small"] = nl.without_H12_k2 == small
    comp_nl = dict(n_small=int(nl.small.sum()), n_entries_small=int(nl[nl.small].entry_id.nunique()),
                   small_by_group=nl[nl.small].group.value_counts().to_dict(),
                   small_top_space_group=nl[nl.small].space_group.value_counts().head(3).to_dict(),
                   small_top_mutation=nl[nl.small].mutation_status.value_counts().head(3).to_dict(),
                   small_peptide_share=float(nl[nl.small].peptide_bound.mean()),
                   rest_peptide_share=float(nl[~nl.small].peptide_bound.mean()))
    sep_nl = []
    for fcol in nolid_cols:
        a, b = nl[~nl.small][fcol].dropna(), nl[nl.small][fcol].dropna()
        cd = c.copy_diffs(per.assign(group="all"), fcol)
        thr = float(np.nanmax([np.sqrt((a.var() * (len(a) - 1) + b.var() * (len(b) - 1)) / (len(a) + len(b) - 2)),
                               np.percentile(cd, 95)]))
        sep_nl.append(dict(feature=fcol, median_rest=float(a.median()), median_small=float(b.median()),
                           diff=float(b.median() - a.median()), threshold=thr,
                           real=bool(abs(b.median() - a.median()) > thr)))
    sep_nl = pd.DataFrame(sep_nl).assign(absdiff=lambda x: x["diff"].abs()).sort_values("absdiff", ascending=False)
    sep_nl.drop(columns="absdiff").to_csv(OUT / "geometry_free_split_without_h12.csv", index=False)
    comp_nl["n_real_features"] = int(sep_nl.real.sum())
    comp_nl["top_features"] = sep_nl.head(4)[["feature", "diff", "threshold", "real"]].to_dict("records")

    # helix length (stability) vs H12 mobility, within groups, entry level
    lm = []
    for grp in ("cl7", "cl4"):
        e = per[per.group == grp].groupby("entry_id")[["len_H12", "mobility_h12_z"]].mean().dropna()
        r, pv = spearmanr(e.len_H12, e.mobility_h12_z)
        lm.append(dict(group=grp, n_entries=len(e), rho=r, p=pv))

    sk = scan.set_index(["features", "k"])
    c.write_numbers(OUT, {
        "phi_box": PHI, "psi_box": PSI, "n_fixed_positions": int(fixed.sum()),
        "n_chains_all_features": n_all, "n_chains_without_h12_features": n_nl, "n_features": len(fcols),
        "n_real_features": int(ft.real.sum()), "real_features": ft[ft.real].feature.tolist(),
        "real_features_without_h12": ft[ft.real & ~ft.feature.str.contains("H12")].feature.tolist(),
        "top5": ft.head(5)[["feature", "median_cl7", "median_cl4", "diff_cl4_minus_cl7", "threshold", "auc"]]
        .to_dict("records"),
        "rf_cv_accuracy_mean": float(np.mean(accs)),
        "max_perm_importance": float(ft.perm_importance.max()),
        "top_importance": ft.sort_values("rf_importance", ascending=False).head(5)[["feature", "rf_importance"]]
        .to_dict("records"),
        "k_best_all": kb_all, "ari_k2_all": float(sk.loc[("all_features", 2), "ari_vs_rmsd_groups"]),
        "ari_k2_all_vs_cl7_cl4": float(sk.loc[("all_features", 2), "ari_vs_cl7_cl4"]),
        "sizes_k2_all": sk.loc[("all_features", 2), "sizes"],
        "k_best_without_h12": kb_nl, "ari_k2_without_h12": float(sk.loc[("without_H12", 2), "ari_vs_rmsd_groups"]),
        "silhouette_k2_all": float(sk.loc[("all_features", 2), "silhouette"]),
        "silhouette_k2_without_h12": float(sk.loc[("without_H12", 2), "silhouette"]),
        "sizes_k2_without_h12": sk.loc[("without_H12", 2), "sizes"],
        "parts": pt.to_dict("records"),
        "split_without_h12": comp_nl,
        "length_vs_mobility": lm,
    })
    pd.set_option("display.width", 250)
    print(pd.DataFrame(chk).to_string())
    print(ft.head(25)[["feature", "median_cl7", "median_cl4", "diff_cl4_minus_cl7", "pooled_sd", "copy_p95", "auc",
                       "real", "rf_importance"]].round(3).to_string())
    print(ft[~ft.feature.str.contains("H12")].head(10)[["feature", "median_cl7", "median_cl4", "diff_cl4_minus_cl7",
                                                        "pooled_sd", "copy_p95", "auc", "real"]].round(3).to_string())
    print(pt.to_string())
    print(scan.round(3).to_string())
    print(lm, np.mean(accs))
    print(comp_nl)
    print(sep_nl.head(8).to_string())


if __name__ == "__main__":
    main()
