"""SV14 figures, drawn from saved CSV/JSON only (rule 1.8).

SV14_pocket_vs_ligand (main): one point per chain of the primary analysis set; x = ligand volume (% of the E2
  ligand reference), y = pocket volume (% of the pocket 100 % reference); colour = S08 cluster (cl7 blue, cl4 orange,
  minor aqua); marker shape and edge = SV12 action class (scripts/sv/common/class_style.py, unlabelled grey);
  grey band = IQR of the pocket reference around 100 %; dotted line = median of the apo entries (shown for reference
  whether or not it is the 100 % reference); dashed lines = equal filling 25 / 50 / 75 % if the whole ligand were
  inside the pocket (pocket % = ligand % x L_ref / (f x P_ref)).
SV14_filling: filling % per entry (entry medians), box + points, by cluster (left) and by action class (right).
Sources: results/sv/SV14/pocket_per_chain.csv, references.json, correlations.csv, tests.csv, filling_entry_level.csv.
Outputs: figures/sv/SV14/*.png + .svg; panels figures/panels/sv/SV14/.
Run: .venv/bin/python -B scripts/sv/SV14_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.class_style import CLASS_COLOR, CLASS_LABEL, CLASS_MARKER, CLASS_ORDER  # noqa: E402
from common.sv_fig import CLUSTER_COLOR, GREY, INK2, MUTED, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV14"
GROUP_LABEL = {"cl7": "cl7", "cl4": "cl4", "minor": "minor clusters"}


def rho_text(cor, subset, level="entry"):
    r = cor[(cor.subset == subset) & (cor.level == level) & (cor.kind == "ligand vs pocket volume")].iloc[0]
    return f"{r.rho:.2f} (q {r.q_bh:.1g}, n {int(r.n_units)})"


def main_scatter(plt, d, ref, cor):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    a = d[d.in_primary_analysis].copy()
    L, P = ref["ligand_100pct"]["median_A3"], ref["pocket_100pct"]["median_A3"]
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    q1, q3 = 100 * ref["pocket_100pct"]["q1_A3"] / P, 100 * ref["pocket_100pct"]["q3_A3"] / P
    ax.axhspan(q1, q3, color=GREY, alpha=0.18, lw=0, zorder=0)
    ax.axhline(100, color=MUTED, lw=0.8, zorder=1)
    ax.axvline(100, color=MUTED, lw=0.8, ls=":", zorder=1)
    apo = ref["apo_baseline"]
    apo_pct = 100 * apo["median_of_entry_medians_A3"] / P
    ax.axhline(apo_pct, color=INK2, lw=0.9, ls=(0, (1, 1.5)), zorder=1)
    xmax = float(np.nanmax(a.ligand_pct_of_E2)) * 1.05
    xs = np.linspace(0, xmax, 50)
    ymax = float(np.nanmax(a.pocket_pct_of_ref)) * 1.08
    for f in (0.25, 0.5, 0.75):
        ys = xs * L / (f * P)
        ax.plot(xs, ys, color=INK2, lw=0.8, ls="--", zorder=1)
        xl = min(xmax * 0.97, ymax * 0.95 * f * P / L)
        ax.text(xl, xl * L / (f * P), f" {int(f * 100)} %", fontsize=7, color=INK2, ha="left", va="bottom")
    # unlabelled first (behind), then labelled classes on top
    order = ["unlabelled"] + [c for c in CLASS_ORDER if c != "unlabelled"]
    counts = {}
    for cls in order:
        s = a[a.action_class.where(a.class_primary_set, "unlabelled") == cls]
        if not len(s):
            continue
        counts[cls] = len(s)
        for g in ("cl7", "cl4", "minor"):
            t = s[s.cluster == g]
            if not len(t):
                continue
            if cls == "unlabelled":
                ax.scatter(t.ligand_pct_of_E2, t.pocket_pct_of_ref, s=13, marker=CLASS_MARKER[cls],
                           facecolors=CLUSTER_COLOR[g], edgecolors=CLASS_COLOR[cls], linewidths=0.4, alpha=0.55,
                           zorder=2)
            else:
                ax.scatter(t.ligand_pct_of_E2, t.pocket_pct_of_ref, s=34, marker=CLASS_MARKER[cls],
                           facecolors=CLUSTER_COLOR[g], edgecolors=CLASS_COLOR[cls], linewidths=1.3, alpha=0.95,
                           zorder=3)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(f"ligand volume (% of E2 = {L:.0f} Å³, median of {ref['ligand_100pct']['n_chains']} wild-type "
                  f"EST chains)")
    src = "wild-type EST chains" if ref["pocket_100pct"]["fallback_used"] else "apo entries"
    ax.set_ylabel(f"pocket volume (% of {P:.0f} Å³, {src})")
    ax.set_title(f"Pocket vs ligand volume, {len(a)} chains: Spearman ρ entries {rho_text(cor, 'all')};\n"
                 f"cl7 {rho_text(cor, 'cl7')}, cl4 {rho_text(cor, 'cl4')}", fontsize=9)
    gh = [Patch(color=CLUSTER_COLOR[g], label=f"{GROUP_LABEL[g]} (n={int((a.cluster == g).sum())})")
          for g in ("cl7", "cl4", "minor") if (a.cluster == g).any()]
    leg1 = ax.legend(handles=gh, title="S08 cluster (fill)", loc="upper left", bbox_to_anchor=(1.02, 1.0),
                     fontsize=7, title_fontsize=7, frameon=False, alignment="left")
    ax.add_artist(leg1)
    ch = [Line2D([], [], ls="", marker=CLASS_MARKER[c], markersize=6 if c != "unlabelled" else 4.5,
                 markerfacecolor="white", markeredgecolor=CLASS_COLOR[c], markeredgewidth=1.3,
                 label=f"{CLASS_LABEL[c]} (n={counts[c]})") for c in CLASS_ORDER if c in counts]
    leg2 = ax.legend(handles=ch, title="action class, SV12 (shape, edge)", loc="upper left",
                     bbox_to_anchor=(1.02, 0.74), fontsize=7, title_fontsize=7, frameon=False, alignment="left")
    ax.add_artist(leg2)
    n_apo = apo["entries_passing"]
    rh = [Patch(color=GREY, alpha=0.3, label=f"IQR of the pocket reference\n({q1:.0f}-{q3:.0f} %)"),
          Line2D([], [], color=INK2, ls=(0, (1, 1.5)), lw=0.9,
                 label=f"apo entries median ({apo_pct:.0f} %,\n{n_apo} entries"
                       + (", < 5: not the reference)" if ref["pocket_100pct"]["fallback_used"] else ")")),
          Line2D([], [], color=INK2, ls="--", lw=0.8, label="equal filling 25/50/75 %\n(whole ligand inside)")]
    ax.legend(handles=rh, title="references", loc="upper left", bbox_to_anchor=(1.02, 0.42), fontsize=7,
              title_fontsize=7, frameon=False, alignment="left")
    return fig, ax


def box_panel(ax, e, grouping, groups, colors, labels, title):
    rng = np.random.default_rng(42)
    data = [e[(e.grouping == grouping) & (e.group == g)].filling_pct.to_numpy() for g in groups]
    bp = ax.boxplot(data, positions=range(len(groups)), widths=0.55, showfliers=False, patch_artist=True,
                    medianprops={"color": INK2, "lw": 1.2}, whiskerprops={"color": MUTED},
                    capprops={"color": MUTED}, boxprops={"edgecolor": MUTED})
    for patch, g in zip(bp["boxes"], groups):
        patch.set_facecolor(colors[g])
        patch.set_alpha(0.25)
    for i, (g, v) in enumerate(zip(groups, data)):
        ax.scatter(i + rng.uniform(-0.18, 0.18, len(v)), v, s=9, color=colors[g], alpha=0.75, lw=0, zorder=3)
    ax.set_xticks(range(len(groups)), [f"{labels[g]}\n(n={len(v)})" for g, v in zip(groups, data)], fontsize=7.5)
    ax.set_ylabel("filling % (entry median)")
    ax.set_title(title, fontsize=8.5)


def filling_fig(plt, e, tests):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.2, 4.0), gridspec_kw={"width_ratios": [3, 4]})
    t = tests.set_index(["comparison", "level"])
    r, rc = t.loc[("cluster", "entry")], t.loc[("cluster", "chain")]
    box_panel(a1, e, "cluster", ["cl7", "cl4", "minor"], CLUSTER_COLOR,
              {"cl7": "cl7", "cl4": "cl4", "minor": "minor\n(not tested)"},
              f"By cluster: entry medians cl7 {r.entry_median_cl7:.1f} %, cl4 {r.entry_median_cl4:.1f} %\n"
              f"Mann-Whitney q {r.q_bh:.2f}; chains bootstrap Δ {rc.effect:.1f} pp, q {rc.q_bh:.2f}")
    present = [c for c in CLASS_ORDER if ((e.grouping == "action_class") & (e.group == c)).any()]
    tested = set(e[(e.grouping == "action_class") & e.tested].group)
    labels = {c: CLASS_LABEL[c].replace(" (no class)", "") + ("" if c in tested else "\n(not tested)")
              for c in present}
    if ("action class", "entry") in t.index:
        r2, r2c = t.loc[("action class", "entry")], t.loc[("action class", "chain")]
        g1, g2 = r2.groups.split(" vs ")
        title = (f"By action class (SV12): entry medians {g1} {r2[f'entry_median_{g1}']:.1f} %, "
                 f"{g2} {r2[f'entry_median_{g2}']:.1f} %\nMann-Whitney q {r2.q_bh:.2f}; chains bootstrap q "
                 f"{r2c.q_bh:.2f}; other classes below the 5-entry gate")
    else:
        title = "By action class (SV12): fewer than 2 classes pass the coverage gate"
    box_panel(a2, e, "action_class", present, CLASS_COLOR, labels, title)
    fig.tight_layout()
    return fig, a1, a2


def main():
    d = pd.read_csv(RES / "pocket_per_chain.csv")
    ref = json.loads((RES / "references.json").read_text())
    cor = pd.read_csv(RES / "correlations.csv")
    tests = pd.read_csv(RES / "tests.csv")
    e = pd.read_csv(RES / "filling_entry_level.csv")
    plt = mpl_setup()
    fig, ax = main_scatter(plt, d, ref, cor)
    paths = save_figure(fig, "SV14", "SV14_pocket_vs_ligand", panels={"scatter": ax})
    fig, a1, a2 = filling_fig(plt, e, tests)
    paths += save_figure(fig, "SV14", "SV14_filling", panels={"by_cluster": a1, "by_class": a2})
    for p in paths:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
