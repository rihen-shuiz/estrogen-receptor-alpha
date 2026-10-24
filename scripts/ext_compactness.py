"""Tracker Y05: compactness of every cluster of every zoo run (Y04) and of the S08 reference partitions.

Input: results/rmsd_entry.npy, results/ext/zoo_labels.csv, results/entry_clusters.csv (S08 cuts 1.0-2.5 A).
Per cluster (noise label -1 skipped; singletons get 0 for intra distances):
  mean_intra_A        mean pairwise RMSD inside the cluster
  medoid_radius_A     mean RMSD to the cluster medoid (chain with the smallest summed RMSD to its mates)
  max_intra_A         largest RMSD to the medoid (radius of the smallest medoid-centred ball holding the cluster)
  diameter_A          largest pairwise RMSD inside the cluster
  separation_A        mean RMSD from the cluster's members to the nearest other cluster (by mean)
  cohesion_separation mean_intra_A / separation_A (lower = tighter and better separated)
Per run: size-weighted mean of each value (results/ext/zoo_compactness_runs.csv).
Outputs: results/ext/zoo_compactness.csv, results/ext/zoo_compactness_runs.csv, figures/ext/Y05_compactness.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_compactness.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CLUSTER_COLOR, GREY, INK2, ORANGE, RES, ROOT, VIOLET, YELLOW, chain_table,
                        cluster_name, mpl_setup, rmsd, save_fig, write_numbers)


def compactness(d, lab, run_id):
    rows = []
    cl = [c for c in sorted(set(lab)) if c >= 0]
    means = {}
    for c in cl:
        m = lab == c
        means[c] = {o: d[np.ix_(m, lab == o)].mean() for o in cl if o != c}
    for c in cl:
        m = np.where(lab == c)[0]
        sub = d[np.ix_(m, m)]
        med = m[np.argmin(sub.sum(1))]
        n = len(m)
        iu = np.triu_indices(n, 1)
        mean_intra = sub[iu].mean() if n > 1 else 0.0
        sep = min(means[c].values()) if means[c] else np.nan
        rows.append(dict(run_id=run_id, cluster=int(c), size=n, medoid=int(med), mean_intra_A=mean_intra,
                         medoid_radius_A=d[med, m].sum() / max(n - 1, 1), max_intra_A=d[med, m].max(),
                         diameter_A=sub.max(), separation_A=sep,
                         cohesion_separation=mean_intra / sep if sep and sep > 0 else np.nan))
    return rows


def main():
    d, idx = rmsd()
    t = chain_table()
    labels = pd.read_csv(RES / "zoo_labels.csv")
    rows = []
    for cut in ("cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5"):
        rows += compactness(d, t[cut].values, f"S08|{cut}")
    for r in labels.columns[1:]:
        rows += compactness(d, labels[r].values, r)
    df = pd.DataFrame(rows)
    df["medoid_chain"] = idx.chain_key.values[df.medoid]
    # matched S08 cluster of each cluster (majority) for colouring
    ref = t.cluster.values
    allab = {**{f"S08|{c}": t[c].values for c in ("cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5")},
             **{r: labels[r].values for r in labels.columns[1:]}}
    df["s08_majority"] = [int(pd.Series(ref[allab[r.run_id] == r.cluster]).mode().iloc[0]) for r in df.itertuples()]
    df.to_csv(RES / "zoo_compactness.csv", index=False)
    w = df.groupby("run_id").apply(lambda g: pd.Series({
        "n_clusters": len(g), "w_mean_intra_A": np.average(g.mean_intra_A, weights=g["size"]),
        "w_medoid_radius_A": np.average(g.medoid_radius_A, weights=g["size"]),
        "w_diameter_A": np.average(g.diameter_A, weights=g["size"]),
        "w_cohesion_separation": np.average(g.cohesion_separation.fillna(0), weights=g["size"])}), include_groups=False)
    runs = pd.read_csv(RES / "zoo_runs.csv")[["run_id", "method", "score", "silhouette", "stability_ari"]]
    w = w.reset_index().merge(runs, on="run_id", how="left")
    w.to_csv(RES / "zoo_compactness_runs.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1.1, 1, 1.1]})
    ax = axs[0]
    s15 = df[df.run_id == "S08|cut_1.5"].sort_values("size", ascending=False)
    x = np.arange(len(s15))
    for off, col, key, lab in ((-0.27, GREY, "mean_intra_A", "mean pairwise RMSD"),
                               (0, VIOLET, "medoid_radius_A", "mean RMSD to medoid"),
                               (0.27, YELLOW, "diameter_A", "diameter (max pairwise)")):
        ax.bar(x + off, s15[key], width=0.26, color=col, label=lab)
    ax.set_xticks(x, [f"cl{c}\nn={n}" for c, n in zip(s15.cluster, s15["size"])])
    ax.set(ylabel="RMSD (A)", title="S08 1.5 A clusters: cl7 is the tightest large cluster")
    ax.legend(fontsize=8)
    for i, r in enumerate(s15.itertuples()):
        ax.text(i - 0.27, r.mean_intra_A + 0.05, f"{r.mean_intra_A:.2f}", ha="center", fontsize=7, color=INK2)
    ax = axs[1]
    good = runs.sort_values("score", ascending=False).groupby("method").head(1).dropna(subset=["score"]).run_id.tolist()
    sel = df[df.run_id.isin(good + ["S08|cut_1.5"])]
    for r in sel.itertuples():
        col = CLUSTER_COLOR.get(r.s08_majority, AQUA)
        ax.scatter(r.size, r.cohesion_separation, s=28, c=col, edgecolor="white", lw=0.5)
    ax.set_xscale("log")
    ax.set(xlabel="cluster size (chains, log scale)", ylabel="cohesion / separation (mean intra / nearest other)",
           title="Every cluster of the best run per family")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=BLUE, label="majority cl7"),
                       Line2D([], [], marker="o", ls="", color=ORANGE, label="majority cl4"),
                       Line2D([], [], marker="o", ls="", color=AQUA, label="majority minor cluster")], fontsize=8)
    ax = axs[2]
    ww = w[w.run_id.isin(good + ["S08|cut_1.5"])].sort_values("w_cohesion_separation")
    yy = np.arange(len(ww))
    ax.barh(yy, ww.w_cohesion_separation, color=[GREY if r.startswith("S08") else VIOLET for r in ww.run_id])
    ax.set_yticks(yy, [r.replace("|", " ") for r in ww.run_id], fontsize=7.5)
    for y, (v, k) in enumerate(zip(ww.w_cohesion_separation, ww.n_clusters)):
        ax.text(v + 0.01, y, f"{v:.2f} (K={int(k)})", va="center", fontsize=7.5, color=INK2)
    ax.invert_yaxis()
    ax.set(xlabel="size-weighted cohesion / separation (lower = more compact)",
           title="Best run of each family + S08 (grey)")
    c7 = s15.set_index("cluster").loc[7]
    c4 = s15.set_index("cluster").loc[4]
    between = d[np.ix_(t.cluster.values == 7, t.cluster.values == 4)].mean()
    fig.suptitle(f"How compact are the clusters? cl7 mean pairwise RMSD {c7.mean_intra_A:.2f} A vs cl4 "
                 f"{c4.mean_intra_A:.2f} A; both far below the mean {between:.1f} A between cl7 and cl4",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y05_compactness", section="A. Structural clustering zoo", task="Y05",
             shows="mean intra RMSD / medoid radius / diameter for S08 clusters; cohesion-separation per cluster and per run",
             takeaway=(f"cl7 (n={int(c7['size'])}) has mean intra RMSD {c7.mean_intra_A:.2f} A and diameter "
                       f"{c7.diameter_A:.1f} A, cl4 (n={int(c4['size'])}) {c4.mean_intra_A:.2f} A / {c4.diameter_A:.1f} A; "
                       f"cohesion/separation {c7.cohesion_separation:.2f} vs {c4.cohesion_separation:.2f}."),
             source="results/ext/zoo_compactness.csv, results/ext/zoo_compactness_runs.csv",
             script="scripts/ext_compactness.py")
    write_numbers("Y05", dict(cl7=c7.drop(["medoid_chain"]).to_dict(), cl4=c4.drop(["medoid_chain"]).to_dict(),
                              cl7_medoid=c7.medoid_chain, cl4_medoid=c4.medoid_chain, mean_cl7_cl4_A=float(between), n_cluster_rows=len(df)))
    print(s15.drop(columns=["run_id"]).round(3).to_string())


if __name__ == "__main__":
    main()
