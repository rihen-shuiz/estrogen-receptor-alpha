"""SV05 (RUN_SV1028): the one SV05 figure, drawn from the saved SV05 tables only.

Inputs: results/sv/SV05/{spread.csv, ari.csv, cluster_sizes.csv, displacement_profile_v4.csv},
results/sv/SV02/consensus_global.csv (helix bands).
Output: figures/sv/SV05/SV05_core_v4.png + .svg, one cropped PNG per graph in figures/panels/sv/SV05/ (rule 1.8).
Colours: core v3 dark grey, core v4 violet (not a cluster colour); cluster bars cl7 blue, cl4 orange, minor green.
Run: .venv/bin/python -B -I scripts/sv/SV05_figures.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import sv_fig as F  # noqa: E402

RES = ROOT / "results/sv/SV05"
V3, V4 = F.INK2, F.VIOLET
SPREAD_MAX_A = 0.5   # SV05 [DEFAULT], drawn as a line


def bands(ax, cons, label=False):
    for r in cons.itertuples():
        ax.axvspan(r.start - 0.5, r.end + 0.5, color=F.GREY, alpha=0.18, lw=0)
        if label:
            ax.text((r.start + r.end) / 2, 1.0, r.helix_name, transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=7, color=F.INK2)


def cluster_colour(s08_label):
    return {7: F.BLUE, 4: F.ORANGE}.get(int(s08_label), F.AQUA)


def main():
    plt = F.mpl_setup()
    sp = pd.read_csv(RES / "spread.csv")
    ari = pd.read_csv(RES / "ari.csv")
    sizes = pd.read_csv(RES / "cluster_sizes.csv")
    prof = pd.read_csv(RES / "displacement_profile_v4.csv")
    cons = pd.read_csv(ROOT / "results/sv/SV02/consensus_global.csv")
    n_cand = int(sp.candidate.sum())
    n_v4 = int(sp.in_core_v4.sum())
    n_v3 = int(sp.in_core_v3.sum())
    jac = len(set(sp.position[sp.in_core_v4]) & set(sp.position[sp.in_core_v3])) / len(
        set(sp.position[sp.in_core_v4]) | set(sp.position[sp.in_core_v3]))
    jac123 = len(set(sp.position[sp.in_lowest_123]) & set(sp.position[sp.in_core_v3])) / len(
        set(sp.position[sp.in_lowest_123]) | set(sp.position[sp.in_core_v3]))

    fig = plt.figure(figsize=(16, 13))
    gs = fig.add_gridspec(4, 2, height_ratios=[2.2, 0.75, 2.0, 1.9], hspace=0.62, wspace=0.18, top=0.93)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :], sharex=ax_a)
    ax_c = fig.add_subplot(gs[2, :], sharex=ax_a)
    ax_d = fig.add_subplot(gs[3, 0])
    ax_e = fig.add_subplot(gs[3, 1])

    # (a) spread per position in the converged candidate frame
    ax = ax_a
    bands(ax, cons, label=True)
    ok = np.isfinite(sp.spread_A)
    ax.plot(sp.position[ok], sp.spread_A[ok], color=F.GREY, lw=1, label="all observed positions")
    c_in = sp[sp.candidate & sp.in_core_v4]
    c_out = sp[sp.candidate & ~sp.in_core_v4]
    ax.scatter(c_in.position, c_in.spread_A, s=14, color=V4, zorder=3, label=f"candidate, spread <= {SPREAD_MAX_A} A (core v4, n={n_v4})")
    ax.scatter(c_out.position, c_out.spread_A, s=16, facecolor="white", edgecolor=V4, lw=1, zorder=3,
               label=f"candidate, spread > {SPREAD_MAX_A} A (n={len(c_out)})")
    ax.axhline(SPREAD_MAX_A, color=F.INK2, lw=0.8, ls="--")
    ax.set_yscale("log")
    ax.set_ylabel("C-alpha spread after iterative fit (A, log)")
    ax.legend(loc="upper left", fontsize=8, ncol=3)
    ax.set_title(f"(a) {n_v4} of {n_cand} candidate positions (consensus helices, observed in >= 98 % of chains, not 528-550) "
                 f"have spread <= {SPREAD_MAX_A} A; candidate spread {sp.spread_A[sp.candidate].min():.2f}-"
                 f"{sp.spread_A[sp.candidate].max():.2f} A", fontsize=10, pad=14)
    plt.setp(ax.get_xticklabels(), visible=True)

    # (b) core membership strip
    ax = ax_b
    bands(ax, cons)
    rows = [("core v3", sp.in_core_v3, V3), ("core v4", sp.in_core_v4, V4), ("lowest-spread 123", sp.in_lowest_123, V4)]
    for k, (name, m, col) in enumerate(rows):
        y = len(rows) - 1 - k
        p = sp.position[m]
        ax.vlines(p, y - 0.35, y + 0.35, color=col, lw=2.2, alpha=1.0 if k < 2 else 0.55)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1], fontsize=8)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_title(f"(b) Core positions: v3 {n_v3}, v4 {n_v4} (Jaccard with v3 {jac:.2f}); "
                 f"lowest-spread 123 (Jaccard with v3 {jac123:.2f})", fontsize=10)
    plt.setp(ax.get_xticklabels(), visible=True)

    # (c) Y10-style displacement profile, v3 (Y10 file) vs v4
    ax = ax_c
    bands(ax, cons)
    ax.plot(prof.position, prof.d_mean_A_v3_y10, color=V3, lw=1.5, label="core v3 (Y10, results/ext)")
    ax.plot(prof.position, prof.d_mean_A, color=V4, lw=1.5, ls="--", label="core v4")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("cl7 vs cl4 mean C-alpha distance (A, symlog)")
    ax.set_xlabel("ESR1 residue number (UniProt P03372)")
    ax.legend(loc="upper left", fontsize=8)
    tail = prof.position >= 528
    s3 = prof.share_sq_v3_y10[tail].sum()
    s4 = prof.share_sq[tail].sum()
    dmax = prof.delta_d_mean_A.abs()
    pmax = int(prof.position[dmax.idxmax()])
    ax.set_title(f"(c) Cluster-mean distance cl7-cl4 per position: 528-550 hold {s3:.1%} (v3) and {s4:.1%} (v4) of the "
                 f"squared difference; largest |v4 - v3| {dmax.max():.2f} A at {pmax}, "
                 f"{dmax[~tail].max():.2f} A before 528", fontsize=10)

    # (d) clusters per cut and ARI against S08
    ax = ax_d
    x = np.arange(len(ari))
    w = 0.38
    ax.bar(x - w / 2, ari.n_clusters_v3_s08, w, color=V3, label="core v3 (S08)")
    ax.bar(x + w / 2, ari.n_clusters_v4, w, color=V4, label="core v4")
    for i, r in ari.iterrows():
        top = max(r.n_clusters_v3_s08, r.n_clusters_v4)
        ax.text(i, top + 0.4, f"ARI {r.ari_v4_vs_s08:.3f}\n{int(r.chains_changed_best_match)} chains moved",
                ha="center", va="bottom", fontsize=8, color=F.INK)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c:.1f} A" for c in ari.cut_A])
    ax.set_xlabel("average-linkage cut")
    ax.set_ylabel("clusters")
    ax.set_ylim(0, max(ari.n_clusters_v3_s08.max(), ari.n_clusters_v4.max()) * 1.35 + 1)
    ax.legend(loc="upper right", fontsize=8)
    ax.set_title(f"(d) Clusters per cut; ARI of v4 against S08 {ari.ari_v4_vs_s08.min():.3f}-{ari.ari_v4_vs_s08.max():.3f}\n"
                 f"pairwise RMSD v4 vs v3: Spearman {ari.spearman_pairwise_v4_vs_v3.iloc[0]:.4f}", fontsize=10)

    # (e) cluster sizes at 1.5 A, coloured by the S08 cluster they match best
    ax = ax_e
    s = sizes[sizes.cut_A == 1.5]
    s3_ = s[s.core == "v3 (S08)"].sort_values("n_chains", ascending=False).reset_index(drop=True)
    s4_ = s[s.core == "v4"].sort_values("n_chains", ascending=False).reset_index(drop=True)
    n = max(len(s3_), len(s4_))
    for k, (tab, off, own) in enumerate([(s3_, -w / 2, True), (s4_, w / 2, False)]):
        for i, r in tab.iterrows():
            lab = r.cluster if own else r.best_match_other_core_cluster
            ax.bar(i + off, r.n_chains, w, color=cluster_colour(lab), edgecolor=V3 if own else V4, lw=1.2,
                   hatch=None if own else "//")
            ax.text(i + off, r.n_chains + 5, str(int(r.n_chains)), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(n))
    ax.set_xticklabels([f"#{i + 1}" for i in range(n)])
    ax.set_xlabel("cluster rank by size at 1.5 A (left bar v3/S08, right hatched bar v4)")
    ax.set_ylabel("chains")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=F.BLUE, label="cl7 / matches cl7"), Patch(color=F.ORANGE, label="cl4 / matches cl4"),
                       Patch(color=F.AQUA, label="minor / matches a minor cluster"),
                       Patch(facecolor="white", edgecolor=V4, hatch="//", label="core v4")], fontsize=8, loc="upper right")
    r15 = ari[ari.cut_A == 1.5].iloc[0]
    ax.set_title(f"(e) Cluster sizes at 1.5 A\nv3 (S08) {r15.sizes_v3_s08.replace(';', '/')}; v4 {r15.sizes_v4.replace(';', '/')}",
                 fontsize=10)

    fig.suptitle("SV05 - core v4 from the stable consensus-helix positions, compared with core v3", x=0.01, ha="left",
                 fontweight="bold", fontsize=12)
    out = F.save_figure(fig, "SV05", "SV05_core_v4", panels={"a_spread": ax_a, "b_membership": ax_b,
                                                            "c_displacement_profile": ax_c, "d_clusters_ari": ax_d,
                                                            "e_cluster_sizes_1_5": ax_e})
    for p in out:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
