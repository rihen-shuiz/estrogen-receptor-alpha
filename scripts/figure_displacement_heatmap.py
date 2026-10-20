"""Tracker E08: heatmap of the ESR1/ESR2 per-residue C-alpha displacement matrix (E07).

Input: results/displacement_matrix.csv (rows = ESR1 positions 305-550, columns = matched chain pairs) and
results/displacement_pairs.csv (column key, in E05 order). Pairs are drawn in that order (grouped by ESR1 entry).
Colour: one-hue sequential blue ramp (docs palette, steps 100-700), clipped at CLIP_A; cells above the clip show the
darkest step. Grey = position unobserved in either chain. A strip above the map marks core v3 helices and H12.
Output: figures/displacement_heatmap.png (static; the table twin is results/displacement_matrix.csv).
Run from the repo root: python scripts/figure_displacement_heatmap.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/displacement_heatmap.png"
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CORE = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
CLIP_A = 5.0
SURFACE, INK, INK_SECONDARY, MUTED, MISSING = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9"
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]  # steps 100-700
CORE_COLOUR, H12_COLOUR = "#2a78d6", "#eb6834"  # categorical slots 1 and 2


def draw(ax_strip, ax, matrix, positions, extra=None):
    """Draw the strip and heatmap; extra(ax, positions) adds overlays (used by E09)."""
    cmap = LinearSegmentedColormap.from_list("blue_seq", BLUE)
    cmap.set_bad(MISSING)
    image = ax.imshow(np.ma.masked_invalid(matrix.T), aspect="auto", cmap=cmap, vmin=0, vmax=CLIP_A,
                      interpolation="nearest", extent=(positions[0] - 0.5, positions[-1] + 0.5, matrix.shape[1], 0))
    ax_strip.set_facecolor(SURFACE)
    for start, end in CORE["core_ranges"]:
        ax_strip.axvspan(start - 0.5, end + 0.5, ymin=0.55, ymax=0.95, color=CORE_COLOUR, lw=0)
    h12 = SETTINGS["h12_range"]
    ax_strip.axvspan(h12[0] - 0.5, h12[1] + 0.5, ymin=0.05, ymax=0.45, color=H12_COLOUR, lw=0)
    ax_strip.text(positions[0] - 2, 0.75, "core v3", ha="right", va="center", fontsize=8, color=INK_SECONDARY)
    ax_strip.text(positions[0] - 2, 0.25, "H12", ha="right", va="center", fontsize=8, color=INK_SECONDARY)
    ax_strip.set_ylim(0, 1)
    ax_strip.axis("off")
    if extra:
        extra(ax, positions)
    ax.set_xlabel("ESR1 position", color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel("Matched ESR1/ESR2 chain pair (grouped by ESR1 entry)", color=INK_SECONDARY, fontsize=9)
    ax.tick_params(colors=MUTED, labelsize=8)
    for side in ax.spines.values():
        side.set_visible(False)
    return image


def main(extra=None, out=OUT, title=None):
    wide = pd.read_csv(ROOT / "results/displacement_matrix.csv", index_col="esr1_pos")
    pairs = pd.read_csv(ROOT / "results/displacement_pairs.csv")
    wide = wide[pairs.column]
    positions = wide.index.to_numpy()
    fig = plt.figure(figsize=(11, 6.2), dpi=150, facecolor=SURFACE)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.05, 1], width_ratios=[1, 0.02], hspace=0.04, wspace=0.03)
    ax_strip, ax = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])
    ax_strip.sharex(ax)
    image = draw(ax_strip, ax, wide.to_numpy(), positions, extra)
    bar = fig.colorbar(image, cax=fig.add_subplot(grid[1, 1]), extend="max")
    bar.set_label(f"C-alpha displacement (Å, clipped at {CLIP_A:g})", color=INK_SECONDARY, fontsize=9)
    bar.ax.tick_params(colors=MUTED, labelsize=8)
    bar.outline.set_visible(False)
    fig.suptitle(title or f"ESR2 vs ESR1 C-alpha displacement after core fit, {wide.shape[1]} matched pairs "
                 f"({pairs.esr1_entry.nunique()} ESR1 / {pairs.esr2_entry.nunique()} ESR2 entries)",
                 x=0.06, ha="left", fontsize=11, color=INK, fontweight="bold")
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    plt.close(fig)
    print(f"{out.relative_to(ROOT)} written; {wide.shape[0]} positions x {wide.shape[1]} pairs")


if __name__ == "__main__":
    main()
