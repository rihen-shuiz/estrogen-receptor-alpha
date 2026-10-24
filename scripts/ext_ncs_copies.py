"""Tracker Y43: copies in one crystal - within-entry vs between-entry conformational variation.

Input: results/rmsd_entry.npy (S04), chain_table() (entry, ligand of interest, S08 cluster).
Pairs of chains are put in three classes: same entry (NCS copies in one asymmetric unit), different entries with the
same ligand of interest, different entries with different ligands. Per entry with >= 2 matrix chains: max within-entry
RMSD, whether its chains fall in different S08 clusters (mixed), and the ligand. Mann-Whitney between pair classes.
Outputs: results/ext/ncs_copies.csv (one row per multi-chain entry), results/ext/ncs_pair_classes.csv (medians per class),
figures/ext/Y43_ncs_copies.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_ncs_copies.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, INK, INK2, ORANGE, RES, VIOLET, YELLOW, chain_table, mpl_setup, rmsd, save_fig,
                        write_numbers)


def main():
    d, idx = rmsd()
    t = chain_table()
    e = t.entry_id.values
    lg = t.ligand_id.fillna("none").values
    iu = np.triu_indices(len(d), 1)
    same_e = e[iu[0]] == e[iu[1]]
    same_l = (lg[iu[0]] == lg[iu[1]]) & (lg[iu[0]] != "none")
    same_c = t.cluster.values[iu[0]] == t.cluster.values[iu[1]]
    v = d[iu]
    cls = np.where(same_e, "same entry (NCS copies)", np.where(same_l, "other entry, same ligand", "other entry, other ligand"))
    pc = pd.DataFrame({"cls": cls, "rmsd": v, "same_cluster": same_c})
    summ = pc.groupby("cls").agg(n_pairs=("rmsd", "size"), median_rmsd_A=("rmsd", "median"),
                                 q90_rmsd_A=("rmsd", lambda x: x.quantile(0.9)), frac_same_cluster=("same_cluster", "mean"))
    summ.to_csv(RES / "ncs_pair_classes.csv")
    rows = []
    for ent, g in t.groupby("entry_id"):
        if len(g) < 2:
            continue
        ii = g.index.values
        sub = d[np.ix_(ii, ii)]
        rows.append(dict(entry_id=ent, n_chains=len(g), max_within_rmsd_A=sub.max(),
                         median_within_rmsd_A=np.median(sub[np.triu_indices(len(ii), 1)]),
                         mixed_clusters=g.cluster.nunique() > 1, clusters=";".join(map(str, sorted(g.cluster.unique()))),
                         ligand=";".join(sorted(g.ligand_id.dropna().unique())) or "none"))
    ncs = pd.DataFrame(rows)
    ncs.to_csv(RES / "ncs_copies.csv", index=False)
    a = pc[pc.cls == "same entry (NCS copies)"].rmsd
    b = pc[pc.cls == "other entry, same ligand"].rmsd
    c = pc[pc.cls == "other entry, other ligand"].rmsd
    p_ab = mannwhitneyu(a, b).pvalue

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.4), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    ax = axs[0]
    bins = np.linspace(0, 6, 61)
    for vals, col, lab in ((c, GREY, "other entry, other ligand"), (b, YELLOW, "other entry, same ligand"),
                           (a, VIOLET, "same entry (NCS copies)")):
        ax.hist(vals, bins=bins, density=True, histtype="stepfilled", alpha=0.45, color=col,
                label=f"{lab} (n={len(vals)}, median {vals.median():.2f} A)")
    ax.set(xlabel="C-alpha RMSD between two chains (A)", ylabel="density", title="Pair classes")
    ax.legend(fontsize=8)
    ax = axs[1]
    ax.hist(ncs.max_within_rmsd_A, bins=np.linspace(0, 8, 41), color=VIOLET)
    ax.axvline(1.5, color=INK, ls="--", lw=1)
    ax.set(xlabel="largest RMSD between copies of one entry (A)", ylabel="entries",
           title=f"{len(ncs)} multi-chain entries; {int((ncs.max_within_rmsd_A > 1.5).sum())} have copies > 1.5 A apart")
    ax = axs[2]
    mixed = ncs[ncs.mixed_clusters]
    mm = mixed.sort_values("max_within_rmsd_A")
    ax.barh(np.arange(len(mm)), mm.max_within_rmsd_A, color=GREY)
    ax.set_yticks(np.arange(len(mm)), [f"{r.entry_id} ({r.ligand}): cl {r.clusters.replace(';', '/')}" for r in mm.itertuples()],
                  fontsize=8)
    ax.set(xlabel="largest RMSD between its copies (A)", title=f"The {len(mixed)} entries whose copies fall in different clusters")
    sl = pc[pc.cls == "other entry, same ligand"]
    ax.text(0.98, 0.02, f"same ligand, other entry: {1 - sl.same_cluster.mean():.0%} of pairs\nin different clusters",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=INK2)
    fig.suptitle(f"Copies in one crystal are near-identical (median {a.median():.2f} A) and closer than other entries with the same "
                 f"ligand ({b.median():.2f} A, p={p_ab:.1g}): most PDB redundancy is copies and repeated ligands", x=0.01,
                 ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y43_ncs_copies", section="K. Redundancy and diversity", task="Y43",
             shows="RMSD distributions for chain pairs from the same entry, other entries with the same ligand, other ligands; largest within-entry RMSD; entries with copies in different clusters",
             takeaway=(f"NCS copies median {a.median():.2f} A, same ligand other entry {b.median():.2f} A, other ligand "
                       f"{c.median():.2f} A; {len(mixed)} of {len(ncs)} multi-chain entries are split across clusters."),
             source="results/ext/ncs_copies.csv, results/ext/ncs_pair_classes.csv", script="scripts/ext_ncs_copies.py")
    write_numbers("Y43", dict(n_multi=len(ncs), n_mixed=len(mixed), med_same_entry=float(a.median()),
                              med_same_ligand=float(b.median()), med_other=float(c.median()), p_same_entry_vs_same_ligand=float(p_ab),
                              mixed_entries=mixed.entry_id.tolist()))
    print(summ, len(ncs), len(mixed), mixed[["entry_id", "clusters", "max_within_rmsd_A", "ligand"]].to_string())


if __name__ == "__main__":
    main()
