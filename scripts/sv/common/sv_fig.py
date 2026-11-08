"""Shared figure helpers for RUN_SV1028 steps (rule 1.8).

- Palette and matplotlib style come from scripts/ext_common.py (the existing figures): cl7 blue, cl4 orange,
  minor clusters green (AQUA).
- save_figure() writes figures/sv/<step>/<name>.png (200 dpi) + .svg and, for every graph (one axes with its labels,
  legend and colorbar), a PNG in figures/panels/sv/<step>/ cropped to the axes' tight bounding box. A panel may
  also be a list of axes that form one graph (heatmap + side bars + colorbar): the crop is their joint bounding box.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "scripts"))
import ext_common as ec  # noqa: E402  (palette + style of the existing figures)

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (ec.BLUE, ec.ORANGE, ec.AQUA, ec.YELLOW, ec.MAGENTA,
                                                          ec.GREEN, ec.VIOLET, ec.RED)
GREY, INK, INK2, MUTED = ec.GREY, ec.INK, ec.INK2, ec.MUTED
CLUSTER_COLOR = {"cl7": BLUE, "cl4": ORANGE, "minor": AQUA}
mpl_setup = ec.mpl_setup


def save_figure(fig, step: str, name: str, panels: dict | None = None):
    """Save the figure (PNG 200 dpi + SVG) and each named axes as its own cropped PNG. Returns the written paths."""
    from matplotlib.transforms import Bbox
    import matplotlib.pyplot as plt

    fdir = ROOT / "figures/sv" / step
    pdir = ROOT / "figures/panels/sv" / step
    fdir.mkdir(parents=True, exist_ok=True)
    pdir.mkdir(parents=True, exist_ok=True)
    out = []
    png, svg = fdir / f"{name}.png", fdir / f"{name}.svg"
    fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(svg, bbox_inches="tight", facecolor="white")
    out += [png, svg]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for pname, ax in (panels or {}).items():
        # one axes, or a list/tuple of axes that together form one graph (e.g. a heatmap with its side bars)
        axes = ax if isinstance(ax, (list, tuple)) else [ax]
        bb = Bbox.union([a.get_tightbbox(renderer) for a in axes]).transformed(fig.dpi_scale_trans.inverted())
        bb = Bbox.from_extents(bb.x0 - 0.05, bb.y0 - 0.05, bb.x1 + 0.05, bb.y1 + 0.05)
        p = pdir / f"{name}__{pname}.png"
        fig.savefig(p, dpi=200, bbox_inches=bb, facecolor="white")
        out.append(p)
    plt.close(fig)
    return out
