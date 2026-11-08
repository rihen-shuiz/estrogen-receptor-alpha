"""SV11 figure (backup), drawn from a saved CSV only (rule 1.8).

SV11_site_share: share of chains with a P2Rank pocket at each site (ligand pocket, AF-2 groove, dimer interface),
cl7 (blue) vs cl4 (orange), chain level; bar labels = chains with such a pocket / chains in the cluster.
Source: results/sv/SV11/site_tests.csv (level chain, outcome 'share with pocket').
Output: figures/sv/SV11/SV11_site_share.png + .svg; panel figures/panels/sv/SV11/SV11_site_share__share.png.
Run: .venv/bin/python -B -I scripts/sv/SV11_figures.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))
from common.sv_fig import BLUE, INK2, ORANGE, mpl_setup, save_figure  # noqa: E402

SITES = ("ligand pocket", "AF-2 groove", "dimer interface")


def main():
    t = pd.read_csv(ROOT / "results/sv/SV11/site_tests.csv")
    t = t[(t.level == "chain") & (t.outcome == "share with pocket")].set_index("site").loc[list(SITES)]
    plt = mpl_setup()
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    x = np.arange(len(SITES))
    w = 0.36
    for dx, g, col in ((-w / 2, "cl7", BLUE), (w / 2, "cl4", ORANGE)):
        vals = t[f"value_{g}"].to_numpy(float)
        ax.bar(x + dx, vals, width=w, color=col, label=f"{g} (n = {int(t[f'n_{g}'].iloc[0])} chains)")
        for xi, v, k, n in zip(x + dx, vals, t[f"k_{g}"], t[f"n_{g}"]):
            ax.text(xi, v + 0.015, f"{int(k)}/{int(n)}", ha="center", va="bottom", fontsize=7.5, color=INK2)
    ax.set_xticks(x, SITES)
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("share of chains with a pocket at the site")
    ax.set_title("P2Rank pockets by site, cl7 vs cl4 (chain level)")
    ax.legend(loc="upper right", fontsize=8)
    cap = "; ".join(f"{s}: cl7 {t.loc[s, 'value_cl7']:.2f}, cl4 {t.loc[s, 'value_cl4']:.2f}" for s in SITES)
    ax.text(0, -0.2, f"Share of chains with a pocket labelled at the site - {cap}.", transform=ax.transAxes,
            fontsize=7, color=INK2, ha="left", va="top", wrap=True)
    paths = save_figure(fig, "SV11", "SV11_site_share", panels={"share": ax})
    for p in paths:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
