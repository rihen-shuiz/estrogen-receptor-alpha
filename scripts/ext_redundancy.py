"""Tracker Y26: redundancy and diversity of the ESR1 structure set - rarefaction, Hill numbers, set cover, representatives.

The project's central question is how many structures are needed to capture the conformational diversity and how
much of the PDB set is repetition. S15/R04 answered it with clusters ("90% of clusters"); here it is made explicit
in three more ways. Inputs: results/rmsd_entry.npy (S04), results/rmsd_old_core.npy (K09), S08 labels,
results/ext/zoo_labels.csv (Y04), chain_table() (tiers, H12 label).
1 Rarefaction (100 random orders, seed 42): (a) clusters discovered vs number of chains for several partitions;
  (b) conformational coverage = fraction of all chains within r of at least one sampled chain, vs sample size.
2 Hill numbers of a partition with cluster shares p: q0 = number of clusters, q1 = exp(Shannon), q2 = 1 / sum p^2
  ("effective number of conformations"), for S08 cuts, the best zoo runs, and subsets (tier A only, one chain per
  entry, old-core matrix re-clustered with the S08 method at 1.5 A).
3 Set cover: greedy maximum coverage; a chain is covered at radius r if some representative is within r A. Number of
  representatives needed to cover 90 / 95 / 99 / 100% of the chains at r = 0.5, 0.75, 1.0, 1.5, 2.0 A.
  Representatives are restricted to chains with a complete H12 (S16 label) - chains without H12 have small RMSD to
  everyone (hubs, Y42) and would cover others trivially; the unrestricted result is reported as well.
  Subsets: all chains, tier A only, one chain per entry (first matrix chain), old-core matrix.
4 Representative sets: greedy set cover order at 1.0 A; farthest-first traversal (start = S12 medoid of cl7) among
  complete-H12 chains; the coverage of each prefix at 1.0 / 1.5 A.
5 Redundancy per S08 cluster: 1 - (representatives needed to cover 95% of the cluster at 1.0 A) / cluster size.
Outputs: results/ext/hill_numbers.csv, results/ext/set_cover.csv, results/ext/representatives_ext.csv,
results/ext/rarefaction.csv, results/ext/redundancy_by_cluster.csv, figures/ext/Y26_rarefaction.png,
figures/ext/Y26_hill_numbers.png, figures/ext/Y26_set_cover.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_redundancy.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RED, RES, ROOT, SEED, VIOLET, YELLOW,
                        chain_table, mpl_setup, rmsd, save_fig, write_numbers)

RADII = [round(x, 2) for x in np.arange(0.3, 2.01, 0.1)]
TARGETS = [0.90, 0.95, 0.99, 1.00]


def hill(labels):
    lab = np.asarray(labels)
    lab = lab[lab >= 0]
    p = np.bincount(lab - lab.min()).astype(float)
    p = p[p > 0] / p.sum()
    return len(p), float(np.exp(-(p * np.log(p)).sum())), float(1 / (p ** 2).sum())


def greedy_cover(d, r, cand, targets=TARGETS):
    """Greedy max-coverage. cand: boolean mask of allowed representatives. Returns order and coverage per step."""
    within = d[:, cand] <= r  # rows = chains to cover, cols = candidates
    cand_idx = np.where(cand)[0]
    covered = np.zeros(len(d), bool)
    order, cov = [], []
    while covered.mean() < 1.0:
        gain = within[~covered].sum(0)
        k = int(np.argmax(gain))
        if gain[k] == 0:
            break
        order.append(int(cand_idx[k]))
        covered |= within[:, k]
        cov.append(covered.mean())
    need = {}
    for tg in targets:
        hit = [i + 1 for i, c in enumerate(cov) if c >= tg - 1e-12]
        need[tg] = hit[0] if hit else np.nan
    return order, cov, need


def farthest_first(d, start, cand, k):
    sel = [start]
    dmin = d[start].copy()
    for _ in range(k - 1):
        score = np.where(cand, dmin, -1)
        nxt = int(np.argmax(score))
        sel.append(nxt)
        dmin = np.minimum(dmin, d[nxt])
    return sel


def main():
    rng = np.random.default_rng(SEED)
    d, idx = rmsd()
    dold = np.load(ROOT / "results/rmsd_old_core.npy")
    t = chain_table()
    zl = pd.read_csv(RES / "zoo_labels.csv")
    complete = (t.h12_label == "complete").values
    tierA = (t.tier == "A").values
    first = ~t.entry_id.duplicated().values
    old15 = fcluster(linkage(squareform(dold, checks=False), "average"), 1.5, "distance")

    # ---- Hill numbers
    parts = {f"S08 {c}": t[c].values for c in ("cut_1.0", "cut_1.5", "cut_2.0", "cut_2.5")}
    for r in ("spectral|sigma=2.0,k=2", "hier_average|k=3", "kmeans|k=2", "hier_average|k=8", "gmm|k=2"):
        parts[r] = zl[r].values
    parts["S08 cut_1.5, old core"] = old15
    rows = []
    for name, lab in parts.items():
        q0, q1, q2 = hill(lab)
        rows.append(dict(run_id=name, subset="all chains", n=len(lab), q0=q0, q1=q1, q2=q2))
    for sname, m in (("tier A", tierA), ("one chain per entry", first), ("complete H12", complete)):
        for c in ("cut_1.0", "cut_1.5"):
            sub = d[np.ix_(m, m)]
            lab = fcluster(linkage(squareform(sub, checks=False), "average"), float(c.split("_")[1]), "distance")
            q0, q1, q2 = hill(lab)
            rows.append(dict(run_id=f"S08 {c} (re-clustered)", subset=sname, n=int(m.sum()), q0=q0, q1=q1, q2=q2))
    hn = pd.DataFrame(rows)
    hn.to_csv(RES / "hill_numbers.csv", index=False)

    # ---- set cover
    sc, reps = [], []
    subsets = {"all chains": (d, np.ones(len(d), bool)), "tier A": (d, tierA), "one chain per entry": (d, first),
               "old-core matrix": (dold, np.ones(len(d), bool))}
    for sname, (D, m) in subsets.items():
        Dm = D[np.ix_(m, m)]
        cm = complete[m]
        for r in RADII:
            for restricted in (True, False):
                order, cov, need = greedy_cover(Dm, r, cm if restricted else np.ones(len(Dm), bool))
                for tg in TARGETS:
                    sc.append(dict(subset=sname, radius_A=r, coverage_target=tg, n_structures=need[tg],
                                   n_chains=int(m.sum()), representatives="complete-H12 only" if restricted else "any chain",
                                   max_coverage=cov[-1] if cov else 0.0))
                if sname == "all chains" and r == 1.0 and restricted:
                    for k, (o, c) in enumerate(zip(order, cov)):
                        reps.append(dict(method="greedy set cover r=1.0", k=k + 1, chain_key=idx.chain_key[o],
                                         cluster=int(t.cluster[o]), tier=t.tier[o], coverage_1p0=c))
    sc = pd.DataFrame(sc)
    sc.to_csv(RES / "set_cover.csv", index=False)
    med7 = int(np.where(idx.chain_key == pd.read_csv(ROOT / "results/medoids.csv").set_index("cluster").loc[7, "chain_key"])[0][0])
    ff = farthest_first(d, med7, complete, 25)
    for k in range(1, 26):
        s = ff[:k]
        reps.append(dict(method="farthest-first (complete H12)", k=k, chain_key=idx.chain_key[ff[k - 1]],
                         cluster=int(t.cluster[ff[k - 1]]), tier=t.tier[ff[k - 1]],
                         coverage_1p0=float((d[:, s].min(1) <= 1.0).mean()),
                         coverage_1p5=float((d[:, s].min(1) <= 1.5).mean())))
    rp = pd.DataFrame(reps)
    rp.to_csv(RES / "representatives_ext.csv", index=False)

    # ---- redundancy per cluster
    red = []
    for c, g in t.groupby("cluster"):
        m = (t.cluster == c).values
        Dm = d[np.ix_(m, m)]
        cand = complete[m] if complete[m].any() else np.ones(m.sum(), bool)
        _, _, need = greedy_cover(Dm, 1.0, cand, targets=[0.95])
        if need[0.95] != need[0.95]:  # complete-H12 chains cannot reach 95%: allow any chain of the cluster
            _, _, need = greedy_cover(Dm, 1.0, np.ones(m.sum(), bool), targets=[0.95])
        red.append(dict(cluster=int(c), size=int(m.sum()), entries=int(g.entry_id.nunique()),
                        reps_95_at_1A=need[0.95], redundancy=1 - need[0.95] / m.sum() if need[0.95] == need[0.95] else np.nan))
    red = pd.DataFrame(red)
    red.to_csv(RES / "redundancy_by_cluster.csv", index=False)

    # ---- rarefaction
    rr = []
    ns = np.unique(np.r_[np.arange(1, 50), np.arange(50, 789, 10), 788])
    curves = {"S08 cut_1.0": t["cut_1.0"].values, "S08 cut_1.5": t["cut_1.5"].values,
              "hier_average|k=8": zl["hier_average|k=8"].values, "spectral|sigma=2.0,k=2": zl["spectral|sigma=2.0,k=2"].values}
    orders = [rng.permutation(len(d)) for _ in range(100)]
    for name, lab in curves.items():
        for n in ns:
            vals = [len(set(lab[o[:n]])) for o in orders]
            rr.append(dict(curve=name, kind="clusters", n=int(n), mean=np.mean(vals), sd=np.std(vals)))
    for r in (1.0, 1.5):
        for n in ns:
            vals = [(d[:, o[:n]].min(1) <= r).mean() for o in orders[:40]]
            rr.append(dict(curve=f"coverage r={r} A (random order)", kind="coverage", n=int(n), mean=np.mean(vals), sd=np.std(vals)))
    rr = pd.DataFrame(rr)
    rr.to_csv(RES / "rarefaction.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 2, figsize=(16, 5.6))
    ax = axs[0]
    for (name, g), col in zip(rr[rr.kind == "clusters"].groupby("curve", sort=False), (VIOLET, INK, YELLOW, AQUA)):
        tot = len(set(curves[name]))
        ax.plot(g.n, g["mean"] / tot, color=col, lw=1.8, label=f"{name} ({tot} clusters)")
        ax.fill_between(g.n, (g["mean"] - g.sd) / tot, (g["mean"] + g.sd) / tot, color=col, alpha=0.15, lw=0)
    ax.axhline(0.9, color=GREY, ls="--", lw=1)
    ax.set_xscale("log")
    ax.set(xlabel="chains sampled at random (log)", ylabel="fraction of clusters discovered (mean +/- SD, 100 orders)",
           title="Cluster rarefaction: rare singletons make 'clusters' saturate late")
    ax.legend(fontsize=8, loc="lower right")
    ax = axs[1]
    for (name, g), col in zip(rr[rr.kind == "coverage"].groupby("curve", sort=False), (VIOLET, YELLOW)):
        ax.plot(g.n, g["mean"], color=col, lw=1.8, label=name)
        ax.fill_between(g.n, g["mean"] - g.sd, g["mean"] + g.sd, color=col, alpha=0.15, lw=0)
    for r, col in ((1.0, VIOLET), (1.5, YELLOW)):
        o = rp[rp.method == "greedy set cover r=1.0"] if r == 1.0 else None
        if o is not None:
            ax.plot(o.k, o.coverage_1p0, color=col, ls=":", lw=2, label="greedy set cover r=1.0 A (complete-H12 reps)")
    ax.axhline(0.9, color=GREY, ls="--", lw=1)
    ax.set_xscale("log")
    ax.set(xlabel="structures chosen (log)", ylabel="fraction of the 788 chains within r of a chosen structure",
           title="Conformational coverage: a few chosen structures cover most chains")
    ax.legend(fontsize=8, loc="lower right")
    n90_15 = rr[(rr.curve == "S08 cut_1.5") & (rr["mean"] / 7 >= 0.9)].n.min()
    cov10 = rr[(rr.kind == "coverage") & (rr.curve.str.contains("1.0"))]
    n90_cov = cov10[cov10["mean"] >= 0.9].n.min()
    fig.suptitle(f"Redundancy: random sampling needs ~{n90_15} chains to see 90% of the 1.5 A clusters, but only "
                 f"~{n90_cov} random chains to place 90% of all chains within 1.0 A of a sampled one", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y26_rarefaction", section="K. Redundancy and diversity", task="Y26",
             shows="rarefaction of cluster discovery for four partitions and of conformational coverage at 1.0 / 1.5 A",
             takeaway=(f"Cluster counts are driven by rare singletons (n90 ~{n90_15} chains at 1.5 A), while coverage of the "
                       f"chains saturates fast (90% within 1.0 A after ~{n90_cov} random chains)."),
             source="results/ext/rarefaction.csv, results/ext/representatives_ext.csv", script="scripts/ext_redundancy.py")

    fig, axs = plt.subplots(1, 2, figsize=(16, 5.8), gridspec_kw={"width_ratios": [1.4, 1]})
    ax = axs[0]
    h = hn.copy()
    h["label"] = h.run_id + " | " + h.subset
    yy = np.arange(len(h))[::-1]
    for k, (q, col, mk) in enumerate((("q0", GREY, "o"), ("q1", VIOLET, "D"), ("q2", YELLOW, "s"))):
        ax.scatter(h[q], yy, c=col, marker=mk, s=40, label={"q0": "q0 = number of clusters", "q1": "q1 = exp(Shannon)",
                                                            "q2": "q2 = inverse Simpson"}[q], zorder=3)
    for y_, a, b in zip(yy, h.q2, h.q0):
        ax.plot([a, b], [y_, y_], color=GREY, lw=0.8, zorder=1)
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 5, 10, 20], ["1", "2", "3", "5", "10", "20"])
    ax.set_yticks(yy, h.label, fontsize=7.5)
    ax.set(xlabel="effective number of conformations (log)", title="Hill numbers: many clusters, ~2 effective conformations")
    ax.legend(fontsize=8, loc="lower right")
    ax = axs[1]
    ax.bar(np.arange(len(red)), red.redundancy, color=[BLUE if c == 7 else ORANGE if c == 4 else AQUA for c in red.cluster])
    for i, r in red.iterrows():
        ax.text(i, (r.redundancy if r.redundancy == r.redundancy else 0) + 0.02, f"{int(r.reps_95_at_1A)}/{int(r['size'])}",
                ha="center", fontsize=8, color=INK2)
    ax.set_xticks(np.arange(len(red)), [f"cl{c}" for c in red.cluster])
    ax.set(ylabel="redundancy = 1 - reps(95% at 1.0 A) / size", ylim=(0, 1.08),
           title="Redundancy per S08 cluster (labels: reps / chains)")
    s15 = hn[(hn.run_id == "S08 cut_1.5") & (hn.subset == "all chains")].iloc[0]
    fig.suptitle(f"Diversity: S08 1.5 A has {s15.q0} clusters but only {s15.q1:.2f} effective conformations (q1) and "
                 f"{s15.q2:.2f} (q2); cl7 and cl4 are both > 90% redundant at 1.0 A", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y26_hill_numbers", section="K. Redundancy and diversity", task="Y26",
             shows="Hill numbers q0/q1/q2 across partitions, methods, cores and subsets; redundancy per S08 cluster",
             takeaway=(f"Whatever the method or subset, q1 stays near {s15.q1:.1f}: the set holds about two conformations; "
                       f"q0 depends on how many singletons a method makes."),
             source="results/ext/hill_numbers.csv, results/ext/redundancy_by_cluster.csv", script="scripts/ext_redundancy.py")

    fig, axs = plt.subplots(1, 2, figsize=(16, 5.8))
    ax = axs[0]
    a = sc[(sc.subset == "all chains") & (sc.representatives == "complete-H12 only")]
    b = sc[(sc.subset == "all chains") & (sc.representatives == "any chain")]
    for tg, col in ((0.90, VIOLET), (0.95, YELLOW), (0.99, RED)):
        g = a[a.coverage_target == tg]
        ax.plot(g.radius_A, g.n_structures, "-o", ms=3.5, color=col, lw=1.8, label=f"{tg:.0%} of chains, complete-H12 reps")
        g2 = b[b.coverage_target == tg]
        ax.plot(g2.radius_A, g2.n_structures, "--", color=col, lw=1, label=f"{tg:.0%}, any chain (hubs allowed)")
    for r_ in (0.5, 1.0):
        v = a[(a.coverage_target == 0.9) & (a.radius_A == r_)].n_structures.iloc[0]
        ax.annotate(f"{int(v)} at {r_} A", (r_, v), xytext=(10, 10), textcoords="offset points", fontsize=8, color=INK2)
    ax.set_yscale("log")
    ax.set(xlabel="coverage radius r (A)", ylabel="representatives needed (greedy set cover, log)",
           title="All 788 chains: representatives needed vs radius (gaps = target not reachable)")
    ax.legend(fontsize=7.5, ncol=2)
    ax = axs[1]
    piv = sc[(sc.coverage_target == 0.9) & (sc.radius_A == 1.0)]
    lab_ = piv.subset + " | " + piv.representatives
    yy = np.arange(len(piv))
    ax.barh(yy, piv.n_structures, color=[VIOLET if "complete" in r else GREY for r in piv.representatives])
    for y_, v_, nn in zip(yy, piv.n_structures, piv.n_chains):
        ax.text(v_ + 0.2, y_, f"{v_:.0f} of {nn}", va="center", fontsize=8, color=INK2)
    ax.set_yticks(yy, lab_, fontsize=8)
    ax.invert_yaxis()
    ax.set(xlabel="representatives for 90% coverage at 1.0 A", title="By subset and core; grey = any chain may represent (hubs allowed)")
    n9010 = a[(a.radius_A == 1.0) & (a.coverage_target == 0.9)].n_structures.iloc[0]
    n9015 = a[(a.radius_A == 1.5) & (a.coverage_target == 0.9)].n_structures.iloc[0]
    n9005 = a[(a.radius_A == 0.5) & (a.coverage_target == 0.9)].n_structures.iloc[0]
    fig.suptitle(f"'90% of the conformational diversity' made explicit: {int(n9010)} structures put 90% of the 788 chains within "
                 f"1.0 A of a representative; at 0.5 A it takes {int(n9005)}", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y26_set_cover", section="K. Redundancy and diversity", task="Y26",
             shows="greedy set-cover size by radius and coverage target; 90% at 1.0 A by subset, core and representative rule",
             takeaway=(f"{int(n9010)} complete-H12 representatives cover 90% of all chains at 1.0 A and {int(n9015)} at 1.5 A; "
                       f"the result barely changes for tier A, one-chain-per-entry or the old core."),
             source="results/ext/set_cover.csv, results/ext/representatives_ext.csv", script="scripts/ext_redundancy.py")
    write_numbers("Y26", dict(hill_s08_15=dict(q0=int(s15.q0), q1=float(s15.q1), q2=float(s15.q2)),
                              n90_r10=float(n9010), n90_r15=float(n9015), n90_r05=float(n9005), n90_clusters_random=float(n90_15),
                              n90_coverage_random=float(n90_cov), redundancy=red.round(4).to_dict("records"),
                              set_cover=sc[sc.coverage_target.isin([0.9, 0.95])].round(4).to_dict("records")))
    print(hn.round(3).to_string())
    print(sc[(sc.subset == "all chains")].pivot_table(index=["representatives", "radius_A"], columns="coverage_target", values="n_structures"))
    print(red)


if __name__ == "__main__":
    main()
