"""Tracker Y21: B-factors by helix, B-group, structural cluster and tier; core-selection figure; K03 with core v3.

Inputs: results/helix_bfactor_by_chain.csv (K03: per-chain normalised C-alpha B per helix = mean B / the chain's
median C-alpha B; b_group tertile), results/lbd_extract/residue_mapping.csv (per-residue C-alpha B),
esr_analysis/resources/config/alignment_core_v3.json (core v3, K04), settings esr1_core_ranges (old core).
K03 computed core_mean_b on the OLD core (84 C-alpha, 340-360/400-420/450-470/500-520). Here, per chain:
core_v3_norm_b (core v3, 123 positions) and old_core_norm_b, both = mean C-alpha B over the core positions observed /
chain median C-alpha B; plus h12_to_core_v3_norm_ratio. Nothing is overwritten (new file results/ext/core_b_v3.csv).
helix_b_summary: per helix (helix_class 1 of results/helix_table.csv), per grouping (b_group / cluster / tier / all):
n chains with the helix, median and IQR of normalised B, rank within the group (1 = most stable = lowest B).
Per-residue profile: median normalised C-alpha B per ESR1 position in cl7 and cl4.
Outputs: results/ext/helix_b_summary.csv, results/ext/core_b_v3.csv, results/ext/bfactor_profile.csv,
figures/ext/Y21_helix_b.png, figures/ext/Y21_core_selection.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_bfactors.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, HELIX_NAMES, INK, INK2, ORANGE, RED, RES, ROOT, SETTINGS, VIOLET, YELLOW,
                        chain_table, helices, helix_bands, mpl_setup, save_fig, write_numbers)


def main():
    t = chain_table()
    hb = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv")
    hb = hb.merge(t[["chain_key", "cluster", "tier"]], on="chain_key", validate="one_to_one")
    hb["cluster_name"] = np.where(hb.cluster == 7, "cl7", np.where(hb.cluster == 4, "cl4", "minor"))
    hel = helices()
    rows = []
    for gb, col in (("all", None), ("b_group", "b_group"), ("cluster", "cluster_name"), ("tier", "tier")):
        groups = [("all", hb)] if col is None else list(hb.groupby(col))
        for g, d in groups:
            med = {h: d[f"{h}_norm_mean_b"].median() for h in hel.helix}
            rk = pd.Series(med).rank()
            for h in hel.helix:
                v = d[f"{h}_norm_mean_b"].dropna()
                rows.append(dict(helix=h, helix_name=HELIX_NAMES.get(h, h), group_by=gb, group=g, n_chains=len(v),
                                 median_norm_b=v.median(), q1=v.quantile(.25), q3=v.quantile(.75), rank=int(rk[h]),
                                 in_core_v3=False))
    s = pd.DataFrame(rows)
    core = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    s["in_core_v3"] = s.helix.isin(core["helices"])
    s.to_csv(RES / "helix_b_summary.csv", index=False)

    rm = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_ca", "B_factor_A2"])
    rm = rm[rm.chain_key.isin(t.chain_key) & rm.has_ca & rm.esr1_position.between(305, 550)].copy()
    med = rm.groupby("chain_key").B_factor_A2.median()
    rm["norm_b"] = rm.B_factor_A2 / rm.chain_key.map(med)
    v3 = set()
    for a, b in core["core_ranges"]:
        v3 |= set(range(a, b + 1))
    old = set()
    for a, b in SETTINGS["esr1_core_ranges"]:
        old |= set(range(a, b + 1))
    h12 = set(range(SETTINGS["h12_range"][0], SETTINGS["h12_range"][1] + 1))
    g = rm.groupby("chain_key")
    cb = pd.DataFrame({
        "core_v3_norm_b": g.apply(lambda d: d.loc[d.esr1_position.isin(v3), "norm_b"].mean(), include_groups=False),
        "core_v3_n_ca": g.apply(lambda d: d.esr1_position.isin(v3).sum(), include_groups=False),
        "old_core_norm_b": g.apply(lambda d: d.loc[d.esr1_position.isin(old), "norm_b"].mean(), include_groups=False),
        "old_core_n_ca": g.apply(lambda d: d.esr1_position.isin(old).sum(), include_groups=False),
        "h12_norm_b": g.apply(lambda d: d.loc[d.esr1_position.isin(h12), "norm_b"].mean(), include_groups=False),
    }).reset_index()
    cb["h12_to_core_v3_norm_ratio"] = cb.h12_norm_b / cb.core_v3_norm_b
    cb = t[["chain_key", "entry_id", "cluster", "tier", "b_group"]].merge(cb, on="chain_key", how="left")
    k03 = hb[["chain_key", "core_norm_b"]].rename(columns={"core_norm_b": "k03_old_core_norm_b"})
    cb = cb.merge(k03, on="chain_key", how="left")
    cb.to_csv(RES / "core_b_v3.csv", index=False)
    check = float(np.nanmax(np.abs(cb.old_core_norm_b - cb.k03_old_core_norm_b)))
    prof = rm.merge(t[["chain_key", "cluster"]], on="chain_key")
    pr = prof[prof.cluster.isin((7, 4))].groupby(["esr1_position", "cluster"]).norm_b.median().unstack()
    pr.columns = [f"median_norm_b_cl{c}" for c in pr.columns]
    pr.reset_index().rename(columns={"esr1_position": "position"}).to_csv(RES / "bfactor_profile.csv", index=False)

    plt = mpl_setup()
    fig = plt.figure(figsize=(19, 9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.25, 1])
    ax = fig.add_subplot(gs[0, 0])
    cols = [("b_group", "low"), ("b_group", "mid"), ("b_group", "high"), ("cluster", "cl7"), ("cluster", "cl4"),
            ("cluster", "minor"), ("tier", "A"), ("tier", "B"), ("tier", "C"), ("all", "all")]
    hl = list(hel.helix)
    M = np.array([[s[(s.helix == h) & (s.group_by == a) & (s.group == b)].median_norm_b.iloc[0] for (a, b) in cols] for h in hl])
    R = np.array([[s[(s.helix == h) & (s.group_by == a) & (s.group == b)]["rank"].iloc[0] for (a, b) in cols] for h in hl])
    N = [int(s[(s.group_by == a) & (s.group == b)].n_chains.max()) for (a, b) in cols]
    im = ax.imshow(M, cmap="esr_seq", vmin=0.7, vmax=1.6, aspect="auto")
    for i in range(len(hl)):
        for j in range(len(cols)):
            ax.text(j, i, f"{M[i, j]:.2f}\n#{R[i, j]}", ha="center", va="center", fontsize=6.5,
                    color="white" if M[i, j] > 1.25 else INK)
    ax.set_xticks(range(len(cols)), [f"{b}\n{a}\nn={n}" for (a, b), n in zip(cols, N)], fontsize=7.5)
    ax.set_yticks(range(len(hl)), [f"{HELIX_NAMES.get(h, h)} ({h}{', core v3' if h in core['helices'] else ''})" for h in hl],
                  fontsize=8)
    for k in (2.5, 5.5, 8.5):
        ax.axvline(k, color="white", lw=2)
    fig.colorbar(im, ax=ax, fraction=0.03, label="median normalised C-alpha B (1 = chain median); #rank, 1 = most stable")
    ax.set_title("Normalised B per helix and group: the same helices are stable in every group")
    ax = fig.add_subplot(gs[1, :])
    helix_bands(ax)
    pp = pd.read_csv(RES / "bfactor_profile.csv")
    ax.plot(pp.position, pp.median_norm_b_cl7, color=BLUE, lw=1.5, label="cl7 median")
    ax.plot(pp.position, pp.median_norm_b_cl4, color=ORANGE, lw=1.5, label="cl4 median")
    ax.axhline(1, color=GREY, lw=0.8, ls="--")
    for a, b in core["core_ranges"]:
        ax.plot([a, b], [0.55, 0.55], color=INK, lw=4, solid_capstyle="butt")
    ax.text(305, 0.58, "core v3", fontsize=8, color=INK)
    ax.set(xlabel="ESR1 residue", ylabel="median normalised C-alpha B", xlim=(304, 551), ylim=(0.45, 2.6))
    ax.legend(fontsize=8, loc="upper center")
    ax = fig.add_subplot(gs[0, 1])
    rng = np.random.default_rng(42)
    for k, (c, col, nm) in enumerate(((7, BLUE, "cl7"), (4, ORANGE, "cl4"))):
        d = cb[cb.cluster == c]
        for j, (v, lab) in enumerate((("core_v3_norm_b", "core v3"), ("old_core_norm_b", "old core"),
                                       ("h12_norm_b", "H12"))):
            x0 = j * 2.6 + k
            vals = d[v].dropna().values
            parts = ax.violinplot([vals], positions=[x0], widths=0.8, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(col); pc.set_alpha(0.3); pc.set_edgecolor("none")
            ax.scatter(x0 + rng.uniform(-.15, .15, len(vals)), vals, s=2, c=col, lw=0, alpha=.5)
            ax.hlines(np.median(vals), x0 - .3, x0 + .3, color=INK, lw=1.3)
    ax.set_xticks([0.5, 3.1, 5.7], ["core v3 (123 pos.)", "old core (84 pos.)", "H12 (538-547)"])
    ax.set(ylabel="mean normalised C-alpha B per chain", title="Per chain: core v3 is quieter than the old core; H12 is not")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=BLUE, label="cl7"), Patch(color=ORANGE, label="cl4")], fontsize=8)
    r7 = cb[cb.cluster == 7].h12_to_core_v3_norm_ratio.median()
    r4 = cb[cb.cluster == 4].h12_to_core_v3_norm_ratio.median()
    fig.suptitle(f"B-factors: H3, H4/5, H7, H8 and H10/11 are the most stable helices in every B-group, cluster and tier; "
                 f"H12 / core-v3 B ratio median {r7:.2f} (cl7) vs {r4:.2f} (cl4)", x=0.01, ha="left",
                 fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y21_helix_b", section="H. B-factors and core", task="Y21",
             shows="median normalised B per helix by B-group, cluster and tier (with ranks); per-chain core v3 / old core / H12 B by cluster; per-residue B profile cl7 vs cl4",
             takeaway=(f"The helix ranking is stable across all groups; H12 is relatively more mobile in cl4 "
                       f"(H12/core-v3 ratio {r4:.2f} vs {r7:.2f} in cl7)."),
             source="results/ext/helix_b_summary.csv, results/ext/core_b_v3.csv, results/ext/bfactor_profile.csv",
             script="scripts/ext_bfactors.py")

    # core-selection figure: rank bump chart over b_groups + chosen core
    fig, axs = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axs[0]
    bg = ["low", "mid", "high"]
    nh = len(hl)
    for h in hl:
        rr = [s[(s.helix == h) & (s.group_by == "b_group") & (s.group == b)]["rank"].iloc[0] for b in bg]
        inc = h in core["helices"]
        col = VIOLET if inc else (RED if h == "AC4" else GREY)
        ax.plot(range(3), rr, "-o", color=col, lw=2.5 if inc else 1.2, ms=6)
        ax.text(2.08, rr[-1], f"{HELIX_NAMES.get(h, h)} ({h})", va="center", fontsize=8, color=col)
    ax.axhline(5.5, color=INK, ls="--", lw=0.8)
    ax.text(-0.1, 5.35, "better half (ranks 1-5)", fontsize=8, ha="left", color=INK2)
    ax.set_xticks(range(3), [f"{b} B-group" for b in bg])
    ax.invert_yaxis()
    ax.set_xlim(-0.15, 2.7)
    ax.set(ylabel="rank of median normalised B (1 = most stable)",
           title="K04 rule: helices in the better half in every B-group (violet) form core v3\n"
                 "(ranks among the 9 main helices here; K04 ranked 11 segments)")
    ax = axs[1]
    ax.scatter(cb.old_core_norm_b, cb.core_v3_norm_b, s=6, c=[BLUE if c == 7 else ORANGE if c == 4 else "#1baf7a"
                                                               for c in cb.cluster], lw=0)
    lim = [0.6, 1.4]
    ax.plot(lim, lim, color=INK, lw=0.8)
    ax.set(xlim=lim, ylim=lim, xlabel="old-core normalised B (per chain)", ylabel="core-v3 normalised B (per chain)",
           title=f"Per chain (n={cb.core_v3_norm_b.notna().sum()}): median old core {cb.old_core_norm_b.median():.3f}, "
                 f"core v3 {cb.core_v3_norm_b.median():.3f}")
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=l) for c, l in
                       ((BLUE, "cl7"), (ORANGE, "cl4"), ("#1baf7a", "minor"))], fontsize=8, loc="upper left")
    fig.suptitle("Core selection: core v3 = the five helices that stay in the more stable half of the ranking in all three "
                 "B-groups (never H12, 123 positions)", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y21_core_selection", section="H. B-factors and core", task="Y21",
             shows="rank of each helix by median normalised B in the low/mid/high B-groups (bump chart) and per-chain core B, core v3 vs old core",
             takeaway=(f"Core v3 helices ({', '.join(HELIX_NAMES.get(h, h) for h in core['helices'])}) stay in the better half "
                       f"in every B-group; per chain the core-v3 normalised B (median {cb.core_v3_norm_b.median():.3f}) "
                       f"is lower than on the old core ({cb.old_core_norm_b.median():.3f})."),
             source="results/ext/helix_b_summary.csv, results/ext/core_b_v3.csv", script="scripts/ext_bfactors.py")
    write_numbers("Y21", dict(core_v3_median=float(cb.core_v3_norm_b.median()), old_core_median=float(cb.old_core_norm_b.median()),
                              h12_core_v3_ratio_cl7=float(r7), h12_core_v3_ratio_cl4=float(r4),
                              h12_core_v3_ratio_all=float(cb.h12_to_core_v3_norm_ratio.median()),
                              old_core_vs_k03_maxdiff=check,
                              rank_all=s[s.group_by == "all"].set_index("helix")["rank"].to_dict()))
    print("old core vs K03 max diff", check)
    print(s[s.group_by == "all"][["helix", "helix_name", "median_norm_b", "rank", "in_core_v3"]].to_string())
    print(r7, r4)


if __name__ == "__main__":
    main()
