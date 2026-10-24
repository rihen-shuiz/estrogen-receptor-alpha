"""Tracker Y27: time - release-year trends, first appearance of each cluster, conformational novelty over time.

Inputs: results/ext/entry_metadata.csv (Y02 release_year; mmCIF first revision date), results/entries.csv
(deposit_date), chain_table() (S08 cluster, resolution, R-free), results/rmsd_entry.npy.
Novelty of a chain = its smallest RMSD to any chain with a complete H12 from an entry released in an EARLIER year
(the first year has none). Chains without a modelled H12 are left out of the reference set because their pairwise RMSD
is computed over few positions and is small to everyone (hubs, Y42); the unrestricted value is kept as
novelty_any_ref_A (it is below 1.5 A for every chain). Per year: chains, chains per cluster, median novelty, fraction novel (> 1.5 A), median resolution and
R-free, cumulative Hill q1 (effective number of conformations among everything released so far).
Outputs: results/ext/time_trends.csv, results/ext/novelty_by_chain.csv, figures/ext/Y27_timeline.png,
figures/ext/Y27_novelty.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_time.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, CLUSTER_COLOR, GREY, INK, INK2, ORANGE, RES, VIOLET, YELLOW, chain_table,
                        mpl_setup, rmsd, save_fig, scatter_clusters, write_numbers)
from ext_redundancy import hill


def main():
    d, idx = rmsd()
    t = chain_table()
    em = pd.read_csv(RES / "entry_metadata.csv")[["entry_id", "release_year"]]
    t = t.merge(em, on="entry_id", how="left")
    yr = t.release_year.values.astype(int)
    nov = np.full(len(t), np.nan)
    nov_any = np.full(len(t), np.nan)
    comp = (t.h12_label == "complete").values
    for i in range(len(t)):
        m = yr < yr[i]
        if m.any():
            nov_any[i] = d[i, m].min()
        if (m & comp).any():
            nov[i] = d[i, m & comp].min()
    t["novelty_A"] = nov  # reference = earlier chains with a complete H12 (hub chains excluded, see Y42)
    t["novelty_any_ref_A"] = nov_any
    t[["chain_key", "entry_id", "release_year", "cluster", "novelty_A", "novelty_any_ref_A"]].to_csv(RES / "novelty_by_chain.csv", index=False)
    rows = []
    for y in range(yr.min(), yr.max() + 1):
        g = t[t.release_year == y]
        upto = t[t.release_year <= y]
        rows.append(dict(release_year=y, n_chains=len(g), n_entries=g.entry_id.nunique(), n_cl7=int((g.cluster == 7).sum()),
                         n_cl4=int((g.cluster == 4).sum()), n_minor=int((~g.cluster.isin((7, 4))).sum()),
                         median_novelty_A=g.novelty_A.median(), frac_novel_gt1p5=(g.novelty_A > 1.5).mean() if len(g) else np.nan,
                         median_resolution_A=g.resolution_A.median(), median_r_free=g.r_free.median(),
                         cum_chains=len(upto), cum_q1=hill(upto["cut_1.5"].values)[1] if len(upto) else np.nan,
                         cum_q0=hill(upto["cut_1.5"].values)[0] if len(upto) else np.nan))
    tt = pd.DataFrame(rows)
    tt.to_csv(RES / "time_trends.csv", index=False)
    first = t.groupby("cluster").release_year.min()

    plt = mpl_setup()
    fig, axs = plt.subplots(3, 1, figsize=(16, 10), sharex=True, gridspec_kw={"height_ratios": [1.7, 0.8, 0.8]})
    ax = axs[0]
    b = np.zeros(len(tt))
    for key, col, lab in (("n_cl7", BLUE, "cl7"), ("n_cl4", ORANGE, "cl4"), ("n_minor", AQUA, "minor clusters")):
        ax.bar(tt.release_year, tt[key], bottom=b, color=col, label=lab, width=0.8)
        b += tt[key].values
    ax.set(ylabel="chains released per year",
           title=f"Releases per year by cluster (cumulative total {int(tt.cum_chains.iloc[-1])} chains); cl4 and cl7 both first seen in {first.get(4)}")
    ax.legend(fontsize=8, loc="upper left")
    ax = axs[1]
    ax.plot(tt.release_year, tt.median_resolution_A, "-o", color=VIOLET, ms=4)
    ax.set(ylabel="median resolution (A)", title="Median resolution per release year")
    ax = axs[2]
    ax.plot(tt.release_year, tt.median_r_free, "-s", color=YELLOW, ms=4)
    ax.set(ylabel="median R-free", xlabel="release year", title="Median R-free per release year")
    peak = tt.loc[tt.n_chains.idxmax()]
    fig.suptitle(f"Timeline: {int(peak.n_chains)} chains in {int(peak.release_year)} alone; the agonist-like cl7 dominates "
                 f"2014-2018, the antagonist-like cl4 grows after 2019", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y27_timeline", section="L. Time", task="Y27",
             shows="chains released per year stacked by cluster with cumulative count; median resolution and R-free per year",
             takeaway=(f"cl4 appears in {first.get(4)} and cl7 in {first.get(7)}; the peak year {int(peak.release_year)} "
                       f"({int(peak.n_chains)} chains) is mostly cl7."),
             source="results/ext/time_trends.csv", script="scripts/ext_time.py")

    fig, axs = plt.subplots(1, 2, figsize=(16, 5.6), gridspec_kw={"width_ratios": [1.5, 1]})
    ax = axs[0]
    rng = np.random.default_rng(42)
    jit = t.release_year + rng.uniform(-0.3, 0.3, len(t))
    scatter_clusters(ax, jit, t.novelty_A, t.cluster, s=8)
    ax.plot(tt.release_year, tt.median_novelty_A, color=INK, lw=1.5, label="median per year")
    ax.axhline(1.5, color=GREY, ls="--", lw=1)
    ax.set_yscale("symlog", linthresh=1)
    ax.set(xlabel="release year", ylabel="novelty = min RMSD to earlier complete-H12 chains (A, symlog)",
           title="Each new structure vs everything released before it")
    ax = axs[1]
    ax.plot(tt.release_year, tt.cum_q1, "-o", color=VIOLET, ms=4, label="cumulative q1 (effective conformations)")
    ax.plot(tt.release_year, tt.cum_q0, "-s", color=GREY, ms=4, label="cumulative q0 (clusters at 1.5 A)")
    ax.set(xlabel="release year", ylabel="number", title="Diversity of everything released so far")
    ax.legend(fontsize=8)
    late = t[t.release_year >= 2010]
    nn = int((late.novelty_A > 1.5).sum())
    fig.suptitle(f"Novelty: since 2010 only {nn} of {len(late)} new chains are > 1.5 A from every earlier complete-H12 chain "
                 f"(all minor-cluster chains); median novelty {late.novelty_A.median():.2f} A", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y27_novelty", section="L. Time", task="Y27",
             shows="novelty of each chain (min RMSD to earlier-released chains) by release year; cumulative Hill q0/q1 over time",
             takeaway=(f"New structures since 2010 almost always repeat a known conformation (median novelty "
                       f"{late.novelty_A.median():.2f} A); the effective number of conformations has been ~2 since both states appeared."),
             source="results/ext/time_trends.csv, results/ext/novelty_by_chain.csv", script="scripts/ext_time.py")
    write_numbers("Y27", dict(first_year=first.to_dict(), peak_year=int(peak.release_year), peak_chains=int(peak.n_chains),
                              frac_novel_since2010=float((late.novelty_A > 1.5).mean()),
                              median_novelty_since2010=float(late.novelty_A.median())))
    print(tt.to_string())


if __name__ == "__main__":
    main()
