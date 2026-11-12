"""Shared look for the figures: white background, thin marks, recessive axes, plain titles, PNG (200 dpi) + SVG.

Colours: cl7 blue and cl4 orange (the project's validated pair); the minor groups grey with a triangle marker, as in
the meeting slides. Text stays in ink colours, never in a series colour.
"""
from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(os.environ.get("ESR1_ROOT", Path(__file__).resolve().parents[2]))
# ESR1_DATA_DIR: read every input by file name from one folder (evidence/data) instead of results/ and data/
DATA = Path(os.environ["ESR1_DATA_DIR"]) if os.environ.get("ESR1_DATA_DIR") else None
OUT = Path(os.environ.get("ESR1_FIG_DIR", ROOT / "figures"))
OUT.mkdir(parents=True, exist_ok=True)

CL7, CL4, MINOR = "#2a78d6", "#eb6834", "#a3a29c"
INK, INK2, MUTED, GRID, BAND = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#f1f0ec"
VIOLET, MAGENTA, BROWN = "#4a3aa7", "#c2185b", "#5d4037"
GROUP_COLOR = {"cl7": CL7, "cl4": CL4, "minor": MINOR}
GROUP_MARKER = {"cl7": "o", "cl4": "o", "minor": "^"}
GROUP_LABEL = {"cl7": "cl7: lid closed over the pocket", "cl4": "cl4: lid in the coactivator groove",
               "minor": "minor groups (grey)"}

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "font.size": 9.5, "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.labelsize": 9.5, "axes.labelcolor": INK, "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8.5,
    "lines.linewidth": 2, "svg.fonttype": "none", "axes.grid": False,
})


def _path(rel: str) -> Path:
    """rel is relative to the repository root (results/<stage>/<file> or data/<topic>/<file>); with ESR1_DATA_DIR set,
    the file of the same name in that folder is read instead."""
    return DATA / Path(rel).name if DATA else ROOT / rel


def csv(rel: str) -> pd.DataFrame:
    return pd.read_csv(_path(rel), keep_default_na=False, na_values=[""])


def numbers(stage: str) -> dict:
    """results/<stage>/numbers.json (in evidence/data: numbers_<stage>.json)."""
    if DATA:
        return json.loads((DATA / f"numbers_{stage}.json").read_text())
    return json.loads((ROOT / "results" / stage / "numbers.json").read_text())


def grid(ax, axis="y"):
    ax.grid(True, axis=axis, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def title(fig_or_ax, text: str, width: int = 95):
    """Plain title that states the result, wrapped; placed on the figure (suptitle) or an axes."""
    t = "\n".join(textwrap.wrap(text, width))
    if hasattr(fig_or_ax, "suptitle"):
        fig_or_ax.suptitle(t, x=0.01, y=1.0, ha="left", va="bottom", fontsize=11, fontweight="bold", color=INK)
    else:
        fig_or_ax.set_title(t, color=INK)


def save(fig, name: str):
    png, svg = OUT / f"{name}.png", OUT / f"{name}.svg"
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    plt.close(fig)
    print("wrote", png.name, svg.name)


def place_labels(ax, xs, ys, texts, base_y, min_gap=0.9, step=0.13, first=0.30, fontsize=7.5):
    """Label points on a horizontal strip without overlaps: labels go above / below base_y in tiers, each tier
    keeping at least min_gap (data units along x) between neighbouring labels."""
    slots = {}  # (side, tier) -> last x
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    for n, i in enumerate(order):
        placed = False
        for tier in range(6):
            for side in ((1, -1) if n % 2 == 0 else (-1, 1)):
                key = (side, tier)
                if key not in slots or xs[i] - slots[key] >= min_gap:
                    slots[key] = xs[i]
                    ty = base_y + side * (first + tier * step)
                    ax.annotate(texts[i], (xs[i], ys[i]), (xs[i], ty), fontsize=fontsize, color=INK2, ha="center",
                                va="center", arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
                    placed = True
                    break
            if placed:
                break


def spread_labels(values, min_gap):
    """Nudge sorted label positions apart so neighbours are at least min_gap apart; returns positions in input order."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    pos = [values[i] for i in order]
    for k in range(1, len(pos)):
        pos[k] = max(pos[k], pos[k - 1] + min_gap)
    out = [0.0] * len(values)
    for k, i in enumerate(order):
        out[i] = pos[k]
    return out
