"""Stage 05 - three more 2-D maps of the 788 chains: Sammon map, self-organising map (SOM) and UMAP on PCA scores.

Declared before running:
  * PCA: same recipe as the saved coordinate PCA (superposed C-alpha, positions observed in >= 95 % of chains,
    missing coordinates replaced by the position mean, 10 components); the first 4 scores must equal the saved ones
    (data/maps/october_coord_pca_scores.csv) up to sign. k = smallest number of components reaching 95 % of the variance.
  * UMAP on the k PCA scores (Euclidean, n_neighbors 15, min_dist 0.1, seed 42).
  * Sammon map of the RMSD matrix (Sammon 1969, pseudo-Newton step with step halving; start = classical MDS;
    at most 500 iterations; zero off-diagonal distances set to 0.01 A, count reported).
  * SOM on the k PCA scores: 12 x 12 rectangular grid, Gaussian neighbourhood, linear initialisation on PC1-PC2,
    batch training 100 epochs, radius 6 -> 0.5; quantisation and topographic error reported; U-matrix saved.
  * Quality, every map: trustworthiness (10 neighbours) against the RMSD matrix and Spearman correlation of map
    distance with RMSD (same measures as the five maps made in October, data/maps/october_map_quality.csv).
  * Structure seen: silhouette of cl7 vs cl4 in the map; silhouette of the two cl4 1.0 A sub-groups; HDBSCAN
    (min cluster size 10) on the map, number of groups and ARI with the 1.5 A groups; where the lid outliers fall
    (distance from their own group's map median, as a percentile of that group's chains).
Outputs: results/05_maps/maps_coords.csv, maps_quality.csv, som_units.csv, maps_outlier_positions.csv,
numbers_maps.json.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA
from sklearn.manifold import trustworthiness
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as c  # noqa: E402

OUT = c.out_dir("05_maps")


def classical_mds(d, k=2):
    n = len(d)
    j = np.eye(n) - 1.0 / n
    b = -0.5 * j @ (d ** 2) @ j
    w, v = np.linalg.eigh(b)
    o = np.argsort(w)[::-1][:k]
    return v[:, o] * np.sqrt(np.maximum(w[o], 0))


def sammon(d, y0, maxiter=500, tol=1e-9, maxhalves=20):
    """Sammon mapping (pseudo-Newton update as in Sammon 1969). Returns coordinates and the final stress."""
    n = len(d)
    scale = 0.5 / d.sum()
    dd = d + np.eye(n)
    dinv = 1.0 / dd
    y = y0.copy()
    one = np.ones((n, 2))

    def dist(y):
        return np.sqrt(((y[:, None, :] - y[None, :, :]) ** 2).sum(-1)) + np.eye(n)

    dy = dist(y)
    e = ((d - (dy - np.eye(n))) ** 2 * dinv)[~np.eye(n, dtype=bool)].sum() * scale
    for it in range(maxiter):
        dyinv = 1.0 / dy
        delta = dyinv - dinv
        deltaone = delta @ one
        g = delta @ y - y * deltaone
        dinv3 = dyinv ** 3
        y2 = y ** 2
        h = dinv3 @ y2 - deltaone - 2 * y * (dinv3 @ y) + y2 * (dinv3 @ one)
        s = -g.ravel() / np.abs(h.ravel())
        y_old = y.copy()
        for _ in range(maxhalves):
            y = y_old + s.reshape(n, 2)
            dy = dist(y)
            e_new = ((d - (dy - np.eye(n))) ** 2 * dinv)[~np.eye(n, dtype=bool)].sum() * scale
            if e_new < e:
                break
            s = s * 0.5
        else:
            y = y_old
            break
        if abs((e - e_new) / e) < tol:
            e = e_new
            break
        e = e_new
    return y, float(e), it + 1


def som_train(x, rows=12, cols=12, epochs=100, r0=6.0, r1=0.5, seed=42):
    """Batch SOM with linear initialisation on the first two principal axes of x."""
    rng = np.random.default_rng(seed)
    grid = np.array([(i, j) for i in range(rows) for j in range(cols)], float)
    mu = x.mean(0)
    u, s, vt = np.linalg.svd(x - mu, full_matrices=False)
    sd = s[:2] / np.sqrt(len(x))
    gi = (grid[:, 0] / (rows - 1) - 0.5) * 4 * sd[0]
    gj = (grid[:, 1] / (cols - 1) - 0.5) * 4 * sd[1]
    w = mu + np.outer(gi, vt[0]) + np.outer(gj, vt[1])
    gd2 = ((grid[:, None, :] - grid[None, :, :]) ** 2).sum(-1)
    for ep in range(epochs):
        r = r0 * (r1 / r0) ** (ep / max(epochs - 1, 1))
        bmu = ((x[:, None, :] - w[None, :, :]) ** 2).sum(-1).argmin(1)
        hk = np.exp(-gd2[bmu] / (2 * r * r))           # (n, units)
        num = hk.T @ x
        den = hk.sum(0)[:, None]
        w = np.where(den > 1e-12, num / np.maximum(den, 1e-12), w)
    dist2 = ((x[:, None, :] - w[None, :, :]) ** 2).sum(-1)
    order = np.argsort(dist2, 1)
    bmu, second = order[:, 0], order[:, 1]
    qe = float(np.sqrt(dist2[np.arange(len(x)), bmu]).mean())
    te = float((np.abs(grid[bmu] - grid[second]).max(1) > 1).mean())
    _ = rng  # deterministic (batch); kept for the declared seed
    return w, grid, bmu, qe, te


def umatrix(w, rows, cols):
    u = np.zeros(rows * cols)
    for i in range(rows):
        for j in range(cols):
            k = i * cols + j
            nb = [(i + a) * cols + (j + b) for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1))
                  if 0 <= i + a < rows and 0 <= j + b < cols]
            u[k] = np.mean([np.linalg.norm(w[k] - w[m]) for m in nb])
    return u


def main():
    t = c.chains()
    d, idx = c.rmsd()
    iu = np.triu_indices(len(d), 1)

    # ---------------- PCA of superposed coordinates (saved recipe) and the reuse check
    z = np.load(c.p("aligned"), allow_pickle=True)
    assert [str(k) for k in z["chain_keys"]] == idx.chain_key.tolist()
    x = z["coords"].astype(float)
    keep = np.isfinite(x[..., 0]).mean(0) >= 0.95
    xs = x[:, keep]
    xs = np.where(np.isfinite(xs), xs, np.nanmean(xs, 0))
    pca = PCA(10, random_state=c.SEED).fit(xs.reshape(len(xs), -1))
    sc = pca.transform(xs.reshape(len(xs), -1))
    saved = pd.read_csv(c.p("october_pca_scores"))
    reuse_diff = max(float(np.abs(np.abs(sc[:, k]) - np.abs(saved[f"cpc{k + 1}"].to_numpy())).max()) for k in range(4))
    cum = np.cumsum(pca.explained_variance_ratio_)
    k = int(np.searchsorted(cum, 0.95) + 1)
    pcs = sc[:, :k]

    # ---------------- UMAP on PCA scores
    import umap
    zu = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=c.SEED).fit_transform(pcs)

    # ---------------- Sammon map of the RMSD matrix
    dz = d.copy()
    off = ~np.eye(len(d), dtype=bool)
    n_zero = int(((dz <= 0) & off).sum() // 2)
    dz[(dz <= 0) & off] = 0.01
    zs, sstress, sit = sammon(dz, classical_mds(dz, 2))

    # ---------------- SOM on PCA scores
    rows = cols = 12
    w, grid, bmu, qe, te = som_train(pcs, rows, cols)
    um = umatrix(w, rows, cols)
    zsom = grid[bmu]

    out = pd.DataFrame({"chain_key": t.chain_key, "entry_id": t.entry_id, "group": t.group,
                        "cluster_1p5": t.cluster_1p5, "cut_1.0": t["cut_1.0"]})
    for i in range(k):
        out[f"pc{i + 1}"] = pcs[:, i]
    out["sammon1"], out["sammon2"] = zs[:, 0], zs[:, 1]
    out["som_row"], out["som_col"], out["som_unit"] = zsom[:, 0].astype(int), zsom[:, 1].astype(int), bmu
    out["umap1"], out["umap2"] = zu[:, 0], zu[:, 1]
    out.to_csv(OUT / "maps_coords.csv", index=False)

    units = pd.DataFrame({"unit": np.arange(rows * cols), "row": grid[:, 0].astype(int),
                          "col": grid[:, 1].astype(int), "u_height": um})
    cnt = out.groupby(["som_unit", "group"]).size().unstack(fill_value=0)
    units = units.merge(cnt, left_on="unit", right_index=True, how="left").fillna(0)
    for g in ("cl7", "cl4", "minor"):
        if g not in units:
            units[g] = 0
    units["n_chains"] = units[["cl7", "cl4", "minor"]].astype(int).sum(1)
    units.to_csv(OUT / "som_units.csv", index=False)

    # ---------------- quality and structure
    sub = pd.read_csv(c.p("cl4_subgroups"))[["chain_key", "sub_cluster"]]
    out = out.merge(sub, on="chain_key", how="left")
    main = out.group.isin(["cl7", "cl4"]).to_numpy()
    cl4m = (out.group == "cl4").to_numpy() & out.sub_cluster.notna().to_numpy()
    fl = pd.read_csv(c.RESULTS / "03_lid_outliers/lid_outliers.csv")
    q_rows, pos_rows = [], []
    rng = np.random.default_rng(c.SEED)
    for name, m2 in (("sammon", zs), ("som", zsom + rng.uniform(-0.3, 0.3, zsom.shape)), ("umap_on_pca", zu)):
        e = np.sqrt(((m2[:, None, :] - m2[None, :, :]) ** 2).sum(-1))
        hd = HDBSCAN(min_cluster_size=10).fit_predict(m2)
        q_rows.append(dict(
            map=name, input="RMSD matrix" if name == "sammon" else f"{k} PCA scores",
            trustworthiness=trustworthiness(d, m2, n_neighbors=10, metric="precomputed"),
            spearman_dist=spearmanr(d[iu], e[iu])[0],
            silhouette_cl7_cl4=silhouette_score(m2[main], out.group[main]),
            silhouette_cl4_subgroups=silhouette_score(m2[cl4m], out.sub_cluster[cl4m]),
            hdbscan_groups=int(len(set(hd)) - (1 if -1 in hd else 0)), hdbscan_noise=int((hd == -1).sum()),
            hdbscan_ari_vs_groups=adjusted_rand_score(out.cluster_1p5, hd),
            stress=sstress if name == "sammon" else np.nan,
            som_quantisation_error=qe if name == "som" else np.nan,
            som_topographic_error=te if name == "som" else np.nan))
        for g in ("cl7", "cl4"):
            gm = (out.group == g).to_numpy()
            ctr = np.median(m2[gm], 0)
            dist_g = np.linalg.norm(m2[gm] - ctr, axis=1)
            for _, r in fl[fl.group.isin([g, "minor"]) if g == "cl7" else fl.group.eq(g)].iterrows():
                if r.group == "minor" and r.compared_with not in ("cl7", "cl4"):
                    continue
                if r.group == "minor" and r.compared_with != g:
                    continue
                i = int(np.where(out.chain_key == r.chain_key)[0][0])
                dd = float(np.linalg.norm(m2[i] - ctr))
                pos_rows.append(dict(map=name, chain_key=r.chain_key, group=r.group, compared_with=g,
                                     lid_state=r.lid_state, map_dist_to_group_median=dd,
                                     percentile_in_group=float((dist_g < dd).mean() * 100)))
    # grey chains nearer cl4 (added here so each grey chain is compared with its nearer main group once)
    for name, m2 in (("sammon", zs), ("som", zsom), ("umap_on_pca", zu)):
        gm = (out.group == "cl4").to_numpy()
        ctr = np.median(m2[gm], 0)
        dist_g = np.linalg.norm(m2[gm] - ctr, axis=1)
        for _, r in fl[(fl.group == "minor") & (fl.compared_with == "cl4")].iterrows():
            i = int(np.where(out.chain_key == r.chain_key)[0][0])
            dd = float(np.linalg.norm(m2[i] - ctr))
            pos_rows.append(dict(map=name, chain_key=r.chain_key, group=r.group, compared_with="cl4",
                                 lid_state=r.lid_state, map_dist_to_group_median=dd,
                                 percentile_in_group=float((dist_g < dd).mean() * 100)))
    # ---------------- islands found by HDBSCAN on each map: what they are and whether they exceed ordinary variation
    same_entry = (t.entry_id.to_numpy()[:, None] == t.entry_id.to_numpy()[None, :])
    same_group = (t.group.to_numpy()[:, None] == t.group.to_numpy()[None, :])
    copy_rmsd = d[np.triu(same_entry & same_group, 1)]
    copy_p95 = float(np.percentile(copy_rmsd, 95))
    isl_rows = []
    for name, m2 in (("sammon", zs), ("som", zsom), ("umap_on_pca", zu)):
        lab = HDBSCAN(min_cluster_size=10).fit_predict(m2)
        for a in sorted(set(lab) - {-1}):
            ia = np.where(lab == a)[0]
            sub = t.iloc[ia]
            main_g = sub.group.value_counts().index[0]
            others = np.where((lab != a) & (lab != -1) & (t.group.to_numpy() == main_g))[0]
            w = d[np.ix_(ia, ia)][np.triu_indices(len(ia), 1)]
            btw = d[np.ix_(ia, others)].ravel() if len(others) else np.array([np.nan])
            within_med = float(np.median(w))
            between_med = float(np.median(btw))
            sg = sub.space_group.value_counts(normalize=True)
            mu = sub.mutation_status.value_counts(normalize=True)
            isl_rows.append(dict(
                map=name, island=int(a), n_chains=len(ia), n_entries=sub.entry_id.nunique(),
                cl7=int((sub.group == "cl7").sum()), cl4=int((sub.group == "cl4").sum()),
                minor=int((sub.group == "minor").sum()), main_group=main_g,
                top_space_group=sg.index[0], top_space_group_share=float(sg.iloc[0]),
                top_mutation=mu.index[0], top_mutation_share=float(mu.iloc[0]),
                peptide_share=float(sub.peptide_bound.mean()), h12_complete_share=float((sub.h12_label == "complete").mean()),
                within_rmsd_median=within_med, between_rmsd_median_same_group=between_med,
                gap_A=between_med - within_med, copy_rmsd_p95=copy_p95,
                exceeds_ordinary=bool(np.isfinite(between_med) and between_med - within_med > copy_p95)))
    isl = pd.DataFrame(isl_rows)
    isl.to_csv(OUT / "maps_islands.csv", index=False)
    lab_u = HDBSCAN(min_cluster_size=10).fit_predict(zu)
    assoc = {col: adjusted_rand_score(t[col].map(str), lab_u) for col in
             ("cut_1.0", "space_group", "mutation_status", "peptide_bound", "entry_id")}

    q = pd.DataFrame(q_rows)
    old = pd.read_csv(c.p("october_map_quality"))
    q = pd.concat([q, old.assign(input="RMSD matrix (October maps)")[["embedding", "input", "trustworthiness",
                                                                      "spearman_dist"]].rename(
        columns={"embedding": "map"})], ignore_index=True)
    q.to_csv(OUT / "maps_quality.csv", index=False)
    pos = pd.DataFrame(pos_rows)
    pos.to_csv(OUT / "maps_outlier_positions.csv", index=False)
    edge = pos.groupby("map").apply(lambda g: float((g.percentile_in_group >= 95).mean()), include_groups=False)
    c.write_numbers(OUT, {
        "pca_components_95pct": k, "pca_cum_var_k": float(cum[k - 1]), "pca_reuse_max_abs_diff": reuse_diff,
        "sammon_stress": sstress, "sammon_iterations": sit, "sammon_zero_pairs_set_0p01": n_zero,
        "som_grid": [rows, cols], "som_qe": qe, "som_te": te,
        "som_units_used": int((units.n_chains > 0).sum()),
        "som_units_mixed_cl7_cl4": int(((units.cl7 > 0) & (units.cl4 > 0)).sum()),
        "quality": q.set_index("map")[["trustworthiness", "spearman_dist"]].to_dict("index"),
        "silhouette_cl7_cl4": q.set_index("map").silhouette_cl7_cl4.dropna().to_dict(),
        "silhouette_cl4_sub": q.set_index("map").silhouette_cl4_subgroups.dropna().to_dict(),
        "hdbscan_groups": q.set_index("map").hdbscan_groups.dropna().to_dict(),
        "hdbscan_ari": q.set_index("map").hdbscan_ari_vs_groups.dropna().to_dict(),
        "flagged_at_group_edge_share": edge.to_dict(),
        "copy_rmsd_p95": copy_p95, "copy_rmsd_median": float(np.median(copy_rmsd)),
        "umap_islands_ari": assoc,
        "islands_exceeding_ordinary": isl[isl.exceeds_ordinary][["map", "island", "n_chains", "top_space_group",
                                                                  "top_mutation", "gap_A"]].to_dict("records"),
        "n_islands": isl.groupby("map").size().to_dict(),
    })
    print(q.to_string())
    print(pos.to_string())
    print("k", k, "reuse diff", reuse_diff, "sammon", sstress, sit, n_zero, "som", qe, te)


if __name__ == "__main__":
    main()
