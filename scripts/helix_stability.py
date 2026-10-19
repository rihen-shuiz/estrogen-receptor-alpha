"""Tracker K04: which helices are most stable per B-factor group, and the alignment core chosen from that.

Inputs: results/helix_bfactor_by_chain.csv (K03), results/helix_table.csv (K02), settings (current core, H12 range).
Per helix and b_group (low, mid, high) and over all chains: n chains, median and IQR of normalised mean B.
Rank 1 = lowest median normalised B within the group, over every helix in the table (H12 included, so it stays visible).
Better half = rank <= number of helices / 2. Recommended core = helices in the better half in every tertile group,
never H12, at least 60 residues in total. The result is saved as core_ranges in
esr_analysis/resources/config/alignment_core_v3.json (new file). analysis_settings.json is not changed.
Also: H12 / core ratio (raw and normalised; median and quartiles) per group, and a heatmap of median normalised B.
Outputs: results/helix_stability_by_bgroup.csv, results/h12_core_ratio_by_bgroup.csv,
esr_analysis/resources/config/alignment_core_v3.json, figures/helix_bfactor_stability.png.
docs/helix_stability.md is written by hand from these files.
Run from the repo root: python scripts/helix_stability.py
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
SETTINGS = json.loads((ROOT / "esr_analysis/resources/config/analysis_settings.json").read_text())
CURRENT_CORE = SETTINGS["esr1_core_ranges"]
H12 = tuple(SETTINGS["h12_range"])
MIN_CORE_RESIDUES = 60
GROUPS = ["low", "mid", "high", "all"]
CORE_OUT = ROOT / "esr_analysis/resources/config/alignment_core_v3.json"
HEATMAP_OUT = ROOT / "figures/helix_bfactor_stability.png"

# Sequential ramp (palette.md, blue steps 100 to 700) and neutral ink, light surface.
RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"


def main():
    chains = pd.read_csv(ROOT / "results/helix_bfactor_by_chain.csv", keep_default_na=False, na_values=[""])
    helices = pd.read_csv(ROOT / "results/helix_table.csv", keep_default_na=False)
    helix_ids = helices.helix.tolist()
    n_helices = len(helix_ids)
    is_h12 = helices.name.eq("H12").to_numpy()

    rows = []
    for group in GROUPS:
        members = chains if group == "all" else chains[chains.b_group == group]
        medians = {}
        for helix in helix_ids:
            values = members[f"{helix}_norm_mean_b"].dropna()
            medians[helix] = (len(values), float(values.median()), float(values.quantile(0.25)), float(values.quantile(0.75)))
        ranks = {}
        for helix in helix_ids:
            ranks[helix] = 1 + sum(medians[other][1] < medians[helix][1] for other in helix_ids)
        for helix, meta in zip(helix_ids, helices.itertuples()):
            n, median, q1, q3 = medians[helix]
            rows.append({"helix": helix, "name": meta.name, "start": meta.start, "end": meta.end, "b_group": group,
                         "n_chains": n, "median_norm_b": median, "q1_norm_b": q1, "q3_norm_b": q3,
                         "iqr_norm_b": q3 - q1, "rank": ranks[helix], "n_helices_ranked": n_helices,
                         "better_half": ranks[helix] <= n_helices / 2, "is_h12": meta.name == "H12"})
    stability = pd.DataFrame(rows)
    stability.to_csv(ROOT / "results/helix_stability_by_bgroup.csv", index=False)

    ratio_rows = []
    for group in GROUPS:
        members = chains if group == "all" else chains[chains.b_group == group]
        for kind, column in [("raw", "h12_to_core_ratio"), ("normalised", "h12_to_core_norm_ratio")]:
            values = members[column].dropna()
            ratio_rows.append({"b_group": group, "ratio": kind, "n_chains": len(values), "median": float(values.median()),
                               "q1": float(values.quantile(0.25)), "q3": float(values.quantile(0.75)),
                               "min": float(values.min()), "max": float(values.max())})
    ratios = pd.DataFrame(ratio_rows)
    ratios.to_csv(ROOT / "results/h12_core_ratio_by_bgroup.csv", index=False)

    tertiles = stability[stability.b_group.isin(["low", "mid", "high"])]
    better_everywhere = tertiles.groupby("helix").better_half.all()
    candidates = [h for h in helix_ids if better_everywhere.get(h, False) and not is_h12[helix_ids.index(h)]]
    sizes = helices.set_index("helix")
    residues = int(sum(sizes.loc[h, "end"] - sizes.loc[h, "start"] + 1 for h in candidates))
    if residues < MIN_CORE_RESIDUES:
        raise ValueError(f"Recommended core has {residues} residues, below {MIN_CORE_RESIDUES}; do not save")
    core_ranges = [[int(sizes.loc[h, "start"]), int(sizes.loc[h, "end"])] for h in candidates]
    config = {
        "name": "alignment_core_v3",
        "core_ranges": core_ranges,
        "helices": candidates,
        "residues_total": residues,
        "rule": "helix in the better half (rank <= n_helices/2) in every b_group (low, mid, high); never H12; >= 60 residues",
        "source": "results/helix_stability_by_bgroup.csv (K04) from results/helix_bfactor_by_chain.csv (K03) and results/helix_table.csv (K02)",
        "previous_core_ranges": CURRENT_CORE,
        "h12_settings_range": list(H12),
        "analysis_settings_changed": False,
    }
    CORE_OUT.write_text(json.dumps(config, indent=2) + "\n")

    fig, axis = plt.subplots(figsize=(7.5, 5.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    grid = np.array([[stability[(stability.helix == h) & (stability.b_group == g)].median_norm_b.iloc[0]
                      for g in GROUPS] for h in helix_ids])
    low, high = np.nanmin(grid), np.nanmax(grid)
    cmap = LinearSegmentedColormap.from_list("blue_ramp", RAMP)
    image = axis.imshow(grid, cmap=cmap, vmin=low, vmax=high, aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            value = grid[i, j]
            fraction = (value - low) / (high - low) if high > low else 0.0
            ink = "#ffffff" if fraction > 0.55 else INK
            axis.text(j, i, f"{value:.2f}", ha="center", va="center", color=ink, fontsize=9)
    labels = [f"{h}{' (H12)' if is_h12[k] else ''}  {int(helices.start.iloc[k])}-{int(helices.end.iloc[k])}"
              for k, h in enumerate(helix_ids)]
    axis.set_yticks(range(len(helix_ids)))
    axis.set_yticklabels(labels, color=INK_SECONDARY, fontsize=9)
    axis.set_xticks(range(len(GROUPS)))
    axis.set_xticklabels([g if g == "all" else f"{g} B" for g in GROUPS], color=INK_SECONDARY)
    axis.set_xlabel("B-factor tertile of chain mean B (all = every eligible chain)", color=INK_SECONDARY)
    axis.set_title("Median normalised C-alpha B per helix", loc="left", color=INK, fontsize=12, fontweight="bold")
    for side in ("top", "right", "bottom", "left"):
        axis.spines[side].set_visible(False)
    colorbar = fig.colorbar(image, ax=axis, fraction=0.04, pad=0.03)
    colorbar.set_label("median normalised B (1 = chain median)", color=INK_SECONDARY)
    colorbar.outline.set_visible(False)
    fig.savefig(HEATMAP_OUT, facecolor=SURFACE, bbox_inches="tight")

    print(stability[stability.b_group == "all"][["helix", "name", "start", "end", "n_chains", "median_norm_b", "iqr_norm_b", "rank"]].to_string(index=False))
    print(ratios.to_string(index=False))
    print(f"candidates (better half in every tertile, not H12): {candidates}; residues {residues}")


if __name__ == "__main__":
    main()
