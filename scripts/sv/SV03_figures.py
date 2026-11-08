"""SV03 figures, drawn only from saved SV03 CSVs (rule 1.8).

1. SV03_h12_map (step 6): x = distance from the H12 centroid to the H11 axis, y = angle H12-H3; one point per chain
   with an H12 axis (H12 = 538-547, >= 6 observed helical C-alpha); colour = S08 cluster, marker = ligand_chemotype.
   Drawn by scripts/sv/common/h12_map.plot_h12_map (reused by SV04 and SV13).
   Source: results/sv/SV03/geometry_per_chain.csv.
2. SV03_descriptors_by_cluster: the 16 tested descriptors (8 uncensored helix lengths, 4 angles, 2 distances, H11
   kink, H11-H12 hinge) by cluster group (cl7, cl4, minor clusters): box (quartiles, whiskers 1.5 IQR) + every chain
   as a jittered point; panel text = the step-4a cl4 vs cl7 rank-biserial r and BH q at chain and entry level.
   Source: results/sv/SV03/geometry_per_chain.csv, results/sv/SV03/tests.csv.
Each graph is also saved as its own cropped PNG in figures/panels/sv/SV03/.
Run: .venv/bin/python -B -I scripts/sv/SV03_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.h12_map import GROUP_COLOR, GROUP_LABEL, GROUP_ORDER, plot_h12_map  # noqa: E402
from common.sv_fig import INK2, mpl_setup, save_figure  # noqa: E402

RES = ROOT / "results/sv/SV03"
STEP = "SV03"
DESC = ["ulen_H1", "ulen_H3", "ulen_H45", "ulen_H7", "ulen_H8", "ulen_H9", "ulen_H1011", "ulen_H12",
        "angle_H12_H3_deg", "angle_H12_H11_deg", "angle_H12_H45_deg", "angle_H11_H3_deg",
        "dist_H12c_H3axis_A", "dist_H12c_H11axis_A", "kink_H11_deg", "hinge_H11C_H12_deg"]


def fmt_q(q):
    return "n.t." if pd.isna(q) else (f"{q:.1e}" if q < 1e-3 else f"{q:.3f}")


def h12_map(plt, w):
    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    fig.subplots_adjust(right=0.72)
    h = plot_h12_map(ax, w)
    n_all = len(w)
    ax.set_title(f"H12 map: {h['n']} of {n_all} chains with an H12 axis (538-547, ≥ 6 helical Cα)",
                 fontsize=9.5)
    ax.text(0.01, 0.99, "core v3 frame; axes from 4-residue Cα window centroids;\nH11 = SV02 consensus H10/11 "
            "(497-526)", transform=ax.transAxes, ha="left", va="top", fontsize=6.8, color=INK2)
    ax.grid(alpha=0.25, lw=0.5)
    return save_figure(fig, STEP, "SV03_h12_map", panels={"h12_map": ax})


def descriptors(plt, w, t):
    a = t[t.part == "a"].set_index(["descriptor", "level"])
    fig, axes = plt.subplots(4, 4, figsize=(16, 15))
    rng = np.random.default_rng(3)
    panels = {}
    for ax, d in zip(axes.ravel(), DESC):
        lab = t.loc[t.descriptor == d, "descriptor_label"].iloc[0]
        data, ticks = [], []
        for i, g in enumerate(GROUP_ORDER):
            v = w.loc[(w.group == g) & w[d].notna(), d].to_numpy()
            data.append(v)
            ticks.append(f"{GROUP_LABEL[g]}\nn={len(v)}")
            if len(v):
                xj = i + rng.uniform(-0.28, 0.28, len(v))
                ax.scatter(xj, v, s=5, color=GROUP_COLOR[g], alpha=0.45, lw=0, zorder=2)
        bp = ax.boxplot([x if len(x) else [np.nan] for x in data], positions=range(3), widths=0.62,
                        showfliers=False, patch_artist=True, zorder=3)
        for patch, g in zip(bp["boxes"], GROUP_ORDER):
            patch.set(facecolor="none", edgecolor=GROUP_COLOR[g], linewidth=1.3)
        for k in ("whiskers", "caps"):
            for ln, g in zip(bp[k], np.repeat(GROUP_ORDER, 2)):
                ln.set(color=GROUP_COLOR[g], linewidth=1.0)
        for ln in bp["medians"]:
            ln.set(color=INK2, linewidth=1.6)
        ax.set_xticks(range(3), ticks, fontsize=7.5)
        ax.set_ylabel(lab, fontsize=8)
        rc, re_ = a.loc[(d, "chain")], a.loc[(d, "entry")]
        ax.set_title(f"{lab}\ncl4 vs cl7 rank-biserial r: chain {rc.effect:+.2f} (q {fmt_q(rc.q)}),\n"
                     f"entry {re_.effect:+.2f} (q {fmt_q(re_.q)})", fontsize=7.4)
        ax.grid(axis="y", alpha=0.25, lw=0.5)
        panels[d] = ax
    fig.suptitle("SV03 helix descriptors by S08 cluster (lengths: uncensored helices only; angles and distances: "
                 "chains with both axes)", fontsize=10.5, y=0.995)
    fig.tight_layout()
    return save_figure(fig, STEP, "SV03_descriptors_by_cluster", panels=panels)


def main():
    plt = mpl_setup()
    w = pd.read_csv(RES / "geometry_per_chain.csv", keep_default_na=False, na_values=[""])
    t = pd.read_csv(RES / "tests.csv", keep_default_na=False, na_values=[""])
    out = h12_map(plt, w) + descriptors(plt, w, t)
    for p in out:
        print(p.relative_to(ROOT), f"{p.stat().st_size / 1024:.1f} kB")


if __name__ == "__main__":
    main()
