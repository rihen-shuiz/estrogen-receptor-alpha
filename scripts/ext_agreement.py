"""Tracker Y06: agreement between clustering methods (ARI / AMI), contingency of S08 vs the selected run, alluvial.

Input: results/ext/zoo_labels.csv, results/ext/zoo_runs.csv, results/ext/zoo_ranking.csv (Y04), S08 labels.
Runs compared: the best setting of each of the 14 families (eligible runs by score, otherwise most clusters) plus the
S08 reference at 1.5 A and 2.5 A. Noise (-1) is kept as its own label in ARI/AMI.
Outputs: results/ext/zoo_agreement.csv (all pairs), figures/ext/Y06_agreement.png (ARI and AMI heatmaps + the
S08 x selected-run contingency), figures/ext/Y06_alluvial.png (chains flowing S08 -> selected -> other methods).
Run from the repo root: PYTHONPATH=. python scripts/ext_agreement.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_mutual_info_score as ami
from sklearn.metrics import adjusted_rand_score as ari

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RES, chain_table, mpl_setup, save_fig,
                        write_numbers)


def main():
    t = chain_table()
    lab = pd.read_csv(RES / "zoo_labels.csv")
    runs = pd.read_csv(RES / "zoo_runs.csv")
    rank = pd.read_csv(RES / "zoo_ranking.csv")
    best = rank.groupby("method").head(1)
    rest = runs[~runs.method.isin(best.method)].sort_values("n_clusters", ascending=False).groupby("method").head(1)
    chosen = list(best.run_id) + list(rest.run_id)
    L = {"S08 1.5 A": t["cut_1.5"].values, "S08 2.5 A": t["cut_2.5"].values}
    for r in chosen:
        L[r] = lab[r].values
    names = list(L)
    rows = []
    for a in names:
        for b in names:
            rows.append(dict(run_a=a, run_b=b, ari=ari(L[a], L[b]), ami=ami(L[a], L[b])))
    df = pd.DataFrame(rows)
    df.to_csv(RES / "zoo_agreement.csv", index=False)
    sel = rank.iloc[0].run_id

    plt = mpl_setup()
    fig = plt.figure(figsize=(19, 7.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 0.75])
    short = [n.split("|")[0] + (" " + n.split("|")[1] if "|" in n else "") for n in names]
    for k, key in enumerate(("ari", "ami")):
        ax = fig.add_subplot(gs[k])
        m = df.pivot(index="run_a", columns="run_b", values=key).loc[names, names].values
        im = ax.imshow(m, cmap="esr_seq", vmin=0, vmax=1)
        ax.set_xticks(range(len(names)), short, rotation=70, ha="right", fontsize=7)
        ax.set_yticks(range(len(names)), short, fontsize=7)
        for i in range(len(names)):
            for j in range(len(names)):
                ax.text(j, i, f"{m[i, j]:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if m[i, j] > 0.6 else INK)
        ax.set_title(f"Pairwise {key.upper()} (1 = identical partitions)")
        fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    ax = fig.add_subplot(gs[2])
    ct = pd.crosstab(pd.Series(L["S08 1.5 A"], name="S08"), pd.Series(L[sel], name="sel"))
    ct = ct.loc[sorted(ct.index, key=lambda c: (c not in (7, 4), c))]
    im = ax.imshow(np.log10(ct.values + 1), cmap="esr_seq", aspect="auto")
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            v = ct.values[i, j]
            ax.text(j, i, str(v), ha="center", va="center", fontsize=9, color="white" if v > 30 else INK)
    ax.set_xticks(range(ct.shape[1]), [f"sel. cluster {c}" for c in ct.columns])
    ax.set_yticks(range(ct.shape[0]), [f"S08 cl{c} (n={ct.loc[c].sum()})" for c in ct.index])
    ax.set_title(f"S08 1.5 A x selected run\n({sel}); cell = chains")
    agree_top = df[(df.run_a == "S08 1.5 A") & df.run_b.isin(best.run_id.head(7))].ari
    fig.suptitle(f"The top-7 families agree with S08 (ARI {agree_top.min():.2f}-{agree_top.max():.2f}); density, "
                 "medoid and single-linkage methods disagree with everyone", x=0.01, ha="left", fontweight="bold",
                 fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y06_agreement", section="A. Structural clustering zoo", task="Y06",
             shows="ARI and AMI between the best run of every family and S08; contingency S08 x selected run",
             takeaway=(f"The selected run ({sel}) puts the 445 cl7 and 329 cl4 chains in separate clusters and "
                       f"absorbs the 14 minor-cluster chains; top-family ARI vs S08 {agree_top.min():.2f}-{agree_top.max():.2f}."),
             source="results/ext/zoo_agreement.csv, results/ext/zoo_labels.csv", script="scripts/ext_agreement.py")

    # ---- alluvial
    order = ["S08 1.5 A", sel] + [r for r in ("kmeans|k=2", "hier_average|k=3", "gmm|k=2", "hdbscan|min_cluster_size=5",
                                              "kmedoids_pam|k=2") if r in L and r != sel]
    ref = L["S08 1.5 A"]
    refcol = np.where(ref == 7, BLUE, np.where(ref == 4, ORANGE, AQUA))
    fig, ax = plt.subplots(figsize=(16, 7))
    gap = 12
    xs = np.arange(len(order)) * 3.0
    blocks = {}
    for xi, name in zip(xs, order):
        labs = L[name]
        cats = sorted(set(labs), key=lambda c: -np.sum(labs == c))
        y = 0
        for c in cats:
            n = int(np.sum(labs == c))
            blocks[(name, c)] = [y, y + n, y]
            ax.add_patch(plt.Rectangle((xi - 0.15, y), 0.3, n, color=GREY if c == -1 else "#4d4c48", lw=0))
            if n >= 12:
                ax.text(xi + 0.2, y + n / 2, f"{'noise' if c == -1 else c} ({n})", va="center", fontsize=7, color=INK2)
            y += n + gap
    for k in range(len(order) - 1):
        a, b = order[k], order[k + 1]
        flows = pd.DataFrame({"a": L[a], "b": L[b], "col": refcol}).groupby(["a", "b", "col"]).size().reset_index(name="n")
        flows["ref_order"] = flows.col.map({BLUE: 0, ORANGE: 1, AQUA: 2})
        flows = flows.sort_values(["a", "ref_order", "b"])
        for f in flows.itertuples():
            ya = blocks[(a, f.a)][2]
            blocks[(a, f.a)][2] += f.n
            yb_key = (b, f.b)
            yb = blocks[yb_key][2] if len(blocks[yb_key]) == 3 else blocks[yb_key][0]
            # track incoming offsets separately
            if len(blocks[yb_key]) == 3:
                blocks[yb_key].append(blocks[yb_key][0])
            yb = blocks[yb_key][3]
            blocks[yb_key][3] += f.n
            xx = np.linspace(xs[k] + 0.15, xs[k + 1] - 0.15, 50)
            s = (1 - np.cos(np.linspace(0, np.pi, 50))) / 2
            lo = ya + (yb - ya) * s
            ax.fill_between(xx, lo, lo + f.n, color=f.col, alpha=0.55, lw=0)
        for key in list(blocks):
            if key[0] == b and len(blocks[key]) == 4:
                blocks[key] = blocks[key][:2] + [blocks[key][0]]
    for xi, name in zip(xs, order):
        ax.text(xi, -12, name.replace("|", "\n"), ha="center", va="bottom", fontsize=8.5, color=INK)
    ax.set_xlim(-0.6, xs[-1] + 1.2)
    ax.set_ylim(-80, 788 + gap * 12)
    ax.invert_yaxis()
    ax.axis("off")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=BLUE, label="S08 cl7 chains (445)"), Patch(color=ORANGE, label="S08 cl4 chains (329)"),
                       Patch(color=AQUA, label="S08 minor clusters (14)")], loc="lower right", fontsize=8)
    ax.set_title("Where do the 788 chains go? Ribbons coloured by S08 cluster: blue and orange stay together "
                 "until density / medoid methods mix them", fontsize=11, loc="left")
    save_fig(fig, "Y06_alluvial", section="A. Structural clustering zoo", task="Y06",
             shows="alluvial of chain membership from S08 through the selected run and four other methods",
             takeaway="The cl7/cl4 ribbons stay separate through spectral, k-means, average linkage and GMM; HDBSCAN "
                      "and PAM merge or interleave them.",
             source="results/ext/zoo_labels.csv", script="scripts/ext_agreement.py")
    write_numbers("Y06", dict(selected=sel, ari_top7_min=float(agree_top.min()), ari_top7_max=float(agree_top.max()),
                              contingency=ct.to_dict()))
    print(ct)


if __name__ == "__main__":
    main()
