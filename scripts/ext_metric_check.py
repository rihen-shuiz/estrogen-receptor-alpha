"""Tracker Y42: is the S04 RMSD a metric, which chains are 'hubs', and does the cl7/cl4 split survive other position sets?

S04 (results/rmsd_entry.npy) is the C-alpha RMSD over the positions observed in BOTH chains of a pair, after the S03
core fit (no refit per pair). Because the position set changes from pair to pair, the triangle inequality can fail
and a chain with few observed positions (e.g. no H12) can look close to everyone ('hub').
Here, from the same superposed coordinates (results/alignment/aligned_coordinates.npz, no refit):
  full     : all LBD positions, pairwise common (reproduces S04; max |diff| reported)
  fixedpos : positions observed in >= 99% of the 788 chains (pairwise common within that set)
  noh12    : all LBD positions except 528-550 (H11-H12 loop, H12, C-terminal tail)
  h12only  : H12 538-547 only (pairs with < 5 common positions get the matrix median, flagged in the CSV)
  core     : alignment core v3 positions only
For each: triangle-violation fraction on 200,000 random triples (seed 42) and the largest excess, Spearman of a
chain's mean distance against its observed LBD positions (hubness), K=2 average-linkage ARI against the S08 1.5 A
labels and against the cl7/cl4 split, silhouette of the S08 main split.
Outputs: results/ext/metric_check.csv, results/ext/rmsd_fixedpos.npy, results/ext/rmsd_noh12.npy,
results/ext/rmsd_h12only.npy, results/ext/rmsd_coreonly.npy, figures/ext/Y42_metric_check.png,
figures/ext/Y42_region_rmsd.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_metric_check.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK2, ORANGE, RED, VIOLET, YELLOW, GREEN, RES, ROOT, SEED, chain_table, coords, mpl_setup, rmsd,
                        save_fig, scatter_clusters, write_numbers)


def masked_rmsd(x, pos_mask, min_common=1):
    """RMSD over positions observed in both chains, restricted to pos_mask. Returns matrix and common counts."""
    x = x[:, pos_mask]
    obs = np.isfinite(x[..., 0])
    n = len(x)
    s = np.zeros((n, n))
    c = np.zeros((n, n))
    for p in range(x.shape[1]):
        o = obs[:, p]
        xp = np.where(o[:, None], x[:, p], 0.0)
        sq = (xp ** 2).sum(1)
        d2 = sq[:, None] + sq[None, :] - 2 * xp @ xp.T
        m = np.outer(o, o)
        s += np.where(m, np.maximum(d2, 0), 0)
        c += m
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.sqrt(s / c)
    r[c < min_common] = np.nan
    np.fill_diagonal(r, 0)
    return r, c


def tri_violation(d, n=200_000):
    rng = np.random.default_rng(SEED)
    i, j, k = rng.integers(0, len(d), (3, n))
    ok = (i != j) & (j != k) & (i != k)
    i, j, k = i[ok], j[ok], k[ok]
    ex = d[i, k] - d[i, j] - d[j, k]
    return float((ex > 1e-6).mean()), float(ex.max())


def main():
    d0, idx = rmsd()
    x, pos = coords()
    t = chain_table()
    core = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    core_pos = set()
    for k in ("core_positions", "positions"):
        if k in core:
            core_pos = set(int(p) for p in core[k])
    if not core_pos:
        for a, b in core.get("ranges", core.get("core_ranges", [])):
            core_pos |= set(range(a, b + 1))
    obs_frac = np.isfinite(x[..., 0]).mean(0)
    sets = {
        "full": np.ones(len(pos), bool),
        "fixedpos": obs_frac >= 0.99,
        "noh12": pos < 528,
        "h12only": (pos >= 538) & (pos <= 547),
        "core": np.isin(pos, sorted(core_pos)),
    }
    ref = t.cluster.values
    main = np.isin(ref, (7, 4))
    nobs = np.isfinite(x[..., 0]).sum(1)
    rows, mats = [], {}
    for name, m in sets.items():
        r, c = masked_rmsd(x, m)
        n_fill = int(np.isnan(r).sum() // 2)
        low = (c < 5) & ~np.eye(len(c), dtype=bool)
        if name == "h12only":
            med = np.nanmedian(r[~low])
            r = np.where(low, med, r)
            n_fill = int(low.sum() // 2)
        mats[name] = r
        viol, ex = tri_violation(r)
        lab2 = fcluster(linkage(squareform(r, checks=False), "average"), 2, "maxclust")
        lab15 = fcluster(linkage(squareform(r, checks=False), "average"), 1.5, "distance")
        rows.append(dict(matrix=name, n_positions=int(m.sum()), triangle_violation_frac=viol, max_excess_A=ex,
                         hub_spearman=spearmanr(r.mean(1), nobs)[0], ari_vs_s08_k2=adjusted_rand_score(ref, lab2),
                         ari_main_split_k2=adjusted_rand_score(ref[main], lab2[main]),
                         n_clusters_1p5=int(lab15.max()), ari_vs_s08_1p5=adjusted_rand_score(ref, lab15),
                         silhouette_main_split=silhouette_score(r[np.ix_(main, main)], ref[main], metric="precomputed"),
                         median_cl7_cl4_A=float(np.median(r[np.ix_(ref == 7, ref == 4)])),
                         median_within_cl7_A=float(np.median(r[np.ix_(ref == 7, ref == 7)][np.triu_indices((ref == 7).sum(), 1)])),
                         median_within_cl4_A=float(np.median(r[np.ix_(ref == 4, ref == 4)][np.triu_indices((ref == 4).sum(), 1)])),
                         n_pairs_filled=n_fill))
    df = pd.DataFrame(rows)
    df["max_abs_diff_vs_s04"] = np.nan
    df.loc[df.matrix == "full", "max_abs_diff_vs_s04"] = float(np.nanmax(np.abs(mats["full"] - d0)))
    df.to_csv(RES / "metric_check.csv", index=False)
    np.save(RES / "rmsd_fixedpos.npy", mats["fixedpos"])
    np.save(RES / "rmsd_noh12.npy", mats["noh12"])
    np.save(RES / "rmsd_h12only.npy", mats["h12only"])
    np.save(RES / "rmsd_coreonly.npy", mats["core"])

    plt = mpl_setup()
    # ---- figure 1: metric check / hubs
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    rng = np.random.default_rng(SEED)
    i, j, k = rng.integers(0, len(d0), (3, 200_000))
    ex = d0[i, k] - d0[i, j] - d0[j, k]
    ax = axs[0]
    ax.hist(ex[ex > 0], bins=60, color=GREY)
    v = df.set_index("matrix").loc["full"]
    ax.set(xlabel="triangle excess d(i,k) - d(i,j) - d(j,k) (A), violating triples only", ylabel="triples",
           title=f"S04 breaks the triangle inequality in {v.triangle_violation_frac:.1%} of triples")
    ax = axs[1]
    hl = t.h12_label.values
    for lab, col in (("complete", VIOLET), ("partial", YELLOW), ("not_modelled", RED)):
        mm = hl == lab
        ax.scatter(nobs[mm], d0.mean(1)[mm], s=10, c=col, label=f"H12 {lab} (n={mm.sum()})", lw=0, alpha=0.8)
    ax.set(xlabel="observed LBD C-alpha positions in the chain", ylabel="mean RMSD to all 787 other chains (A)",
           title="Chains without H12 sit close to everyone ('hubs')")
    ax.legend(fontsize=8)
    ax = axs[2]
    o = ["full", "fixedpos", "noh12", "core", "h12only"]
    dd = df.set_index("matrix").loc[o]
    yy = np.arange(len(o))
    ax.barh(yy - 0.2, dd.ari_vs_s08_1p5, height=0.38, color=GREEN, label="ARI vs S08 (avg linkage, 1.5 A cut)")
    ax.barh(yy + 0.2, dd.silhouette_main_split, height=0.38, color=VIOLET, label="silhouette of the cl7/cl4 split")
    ax.set_yticks(yy, [f"{n} ({int(dd.loc[n, 'n_positions'])} pos)" for n in o])
    for y, (a, s) in enumerate(zip(dd.ari_vs_s08_1p5, dd.silhouette_main_split)):
        ax.text(max(a, 0) + 0.02, y - 0.2, f"{a:.2f}", va="center", fontsize=8, color=INK2)
        ax.text(max(s, 0) + 0.02, y + 0.2, f"{s:.2f}", va="center", fontsize=8, color=INK2)
    ax.set_xlim(-0.1, 1.2)
    ax.set(xlabel="value", title="The split lives in H12 and the loop before it")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1)
    fig.suptitle("The pairwise-common-residue RMSD is not a metric; the two clusters are an H12 effect",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    nh = df.set_index("matrix").loc["noh12"]
    save_fig(fig, "Y42_metric_check", section="C. Why two clusters", task="Y42",
             shows="triangle-inequality excess, hubness of incomplete chains, split strength per position set",
             takeaway=(f"{v.triangle_violation_frac:.1%} of triples violate the triangle inequality (max excess "
                       f"{v.max_excess_A:.1f} A); without positions 528-550 the cl7/cl4 split nearly disappears "
                       f"(silhouette {nh.silhouette_main_split:.2f}; the 1.5 A cut gives {int(nh.n_clusters_1p5)} cluster)."),
             source="results/ext/metric_check.csv", script="scripts/ext_metric_check.py")

    # ---- figure 2: region-restricted maps
    from sklearn.manifold import MDS
    fig, axs = plt.subplots(1, 4, figsize=(17, 4.3))
    for ax, name in zip(axs, ["fixedpos", "noh12", "core", "h12only"]):
        z = MDS(n_components=2, dissimilarity="precomputed", random_state=SEED, n_init=2, max_iter=300,
                normalized_stress="auto").fit_transform(mats[name])
        scatter_clusters(ax, z[:, 0], z[:, 1], ref, s=7, legend=(name == "fixedpos"))
        r = df.set_index("matrix").loc[name]
        ax.set_title(f"{name}: {int(r.n_positions)} positions\nmedian cl7-cl4 {r.median_cl7_cl4_A:.2f} A, "
                     f"within cl7 {r.median_within_cl7_A:.2f}, cl4 {r.median_within_cl4_A:.2f}", fontsize=8.5)
        ax.set_xlabel("MDS 1 (A)"); ax.set_ylabel("MDS 2 (A)")
    fig.suptitle("Same superposition, different residues: cl7 and cl4 separate only when H12 is included; "
                 "on fixed positions, the core or without H12 they overlap", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    h = df.set_index("matrix").loc["h12only"]
    save_fig(fig, "Y42_region_rmsd", section="C. Why two clusters", task="Y42",
             shows="metric MDS of four region-restricted RMSD matrices coloured by S08 cluster",
             takeaway=(f"Median cl7-cl4 distance is {h.median_cl7_cl4_A:.1f} A on H12 alone but "
                       f"{nh.median_cl7_cl4_A:.2f} A without H12/loop: the split is an H12 position split."),
             source="results/ext/metric_check.csv, results/ext/rmsd_noh12.npy, results/ext/rmsd_h12only.npy",
             script="scripts/ext_metric_check.py")
    write_numbers("Y42", df.set_index("matrix").round(4).to_dict("index"))
    print(df.round(3).to_string())


if __name__ == "__main__":
    main()
