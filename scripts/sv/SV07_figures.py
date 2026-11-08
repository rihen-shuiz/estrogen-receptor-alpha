"""SV07 (RUN_SV1028): spike profile figure, drawn from the saved SV07 tables only.

Inputs: results/sv/SV07/{restricted_profile.csv, verdict.csv, spikes.csv, numbers.json},
results/sv/SV02/consensus_global.csv (helix bands, grey).
(a) Y10 cl7-vs-cl4 C-alpha distance profile, two curves per measure on one axis: all chains (as in Y10) and only the
    chains observing the position in both clusters (SV07 step 2b); cluster means and complete-H12 medoid pairs; the
    1.0 A spike threshold; 528-550 (outside the spike rule) shaded; every spike marked and labelled with its verdict.
(b) share of cl7 / cl4 chains observing each position (SV00 mask).
Colours: cl7 blue, cl4 orange (rule 1.8); curves violet / dark grey; verdict colours red, yellow, magenta, grey
(no cluster colour reused). Output: figures/sv/SV07/SV07_spike_profile.png + .svg, panels in figures/panels/sv/SV07/.
Run: .venv/bin/python -B scripts/sv/SV07_figures.py
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/sv"))

import pandas as pd  # noqa: E402

from common import sv_fig as F  # noqa: E402

RES = ROOT / "results/sv/SV07"
THRESH_A = 1.0
EXCLUDE = (528, 550)
VERDICT_STYLE = {"genuine loop difference": (F.RED, "genuine loop difference"),
                 "driven by <= 3 entries": (F.YELLOW, "driven by <= 3 entries"),
                 "numbering error": (F.MAGENTA, "numbering error"),
                 "sparse-sample effect": (F.GREY, "sparse-sample effect")}
SHORT = {"genuine loop difference": "genuine", "driven by <= 3 entries": "<= 3 entries",
         "numbering error": "numbering", "sparse-sample effect": "sparse sample"}


def bands(ax, cons, label=False):
    for r in cons.itertuples():
        ax.axvspan(r.start - 0.5, r.end + 0.5, color=F.GREY, alpha=0.16, lw=0)
        if label:
            ax.text((r.start + r.end) / 2, 1.0, r.helix_name, transform=ax.get_xaxis_transform(), ha="center",
                    va="bottom", fontsize=7, color=F.INK2)


def groups(verd):
    """Consecutive spike positions with the same verdict and the same 'medoid only' status -> one label."""
    out = []
    for r in verd.sort_values("position").itertuples():
        mo = r.exceeds == "medoid"
        if out and r.position == out[-1]["pos"][-1] + 1 and r.verdict == out[-1]["verdict"] \
                and mo == out[-1]["medoid_only"]:
            out[-1]["pos"].append(r.position)
        else:
            out.append(dict(pos=[r.position], verdict=r.verdict, medoid_only=mo))
    return out


def place(labels, levels, res_per_char=0.95, gap=2.0):
    """Greedy: each label (x, text) goes to the lowest level where it does not overlap an earlier label."""
    used = {lv: [] for lv in levels}
    out = []
    for x, txt in labels:
        half = len(txt) * res_per_char / 2
        for lv in levels:
            if all(x + half + gap < a or x - half - gap > b for a, b in used[lv]):
                used[lv].append((x - half, x + half))
                out.append(lv)
                break
        else:
            out.append(levels[-1])
    return out


def main():
    plt = F.mpl_setup()
    prof = pd.read_csv(RES / "restricted_profile.csv")
    verd = pd.read_csv(RES / "verdict.csv")
    sp = pd.read_csv(RES / "spikes.csv")
    n = json.loads((RES / "numbers.json").read_text())
    cons = pd.read_csv(ROOT / "results/sv/SV02/consensus_global.csv")
    rv = n["restricted_vs_y10"]

    fig = plt.figure(figsize=(16, 9.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[2.4, 1.0], hspace=0.28, top=0.84)
    ax = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1], sharex=ax)

    # (a) profile
    bands(ax, cons, label=True)
    ax.axvspan(EXCLUDE[0] - 0.5, EXCLUDE[1] + 0.5, color=F.INK2, alpha=0.07, lw=0)
    ax.text((EXCLUDE[0] + EXCLUDE[1]) / 2, 0.97, "528-550\nnot in the\nspike rule", transform=ax.get_xaxis_transform(),
            ha="center", va="top", fontsize=7, color=F.INK2)
    ax.plot(prof.position, prof.d_mean_y10_A, color=F.VIOLET, lw=2.2, label="cluster means, all chains (Y10)")
    ax.plot(prof.position, prof.d_mean_restricted_A, color=F.INK, lw=1.0, ls=(0, (4, 3)),
            label="cluster means, chains observing the position in both clusters")
    ax.plot(prof.position, prof.d_medoid_y10_A, color=F.MUTED, lw=1.2,
            label=f"complete-H12 medoids {n['y10_medoids'][0]} vs {n['y10_medoids'][1]} (Y10)")
    ax.plot(prof.position, prof.d_medoid_restricted_A, color=F.INK2, lw=0.9, ls=(0, (1, 1.5)),
            label="complete-H12 medoids re-picked among chains observing the position")
    ax.axhline(THRESH_A, color=F.INK2, lw=0.8, ls="--")
    ax.text(306, THRESH_A * 1.04, "1.0 A spike threshold", fontsize=7, color=F.INK2, va="bottom")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylim(0, max(prof.d_mean_y10_A.max(), prof.d_medoid_restricted_A.max()) * 6)
    ax.set_ylabel("cl7 vs cl4 C-alpha distance (A, symlog)")
    # spike markers
    for r in verd.itertuples():
        s = sp[sp.position == r.position].iloc[0]
        col = VERDICT_STYLE[r.verdict][0]
        y = s.d_mean_y10_A if "mean" in r.exceeds else s.d_medoid_y10_A
        ax.scatter([r.position], [y], s=34, color=col, edgecolor=F.INK, linewidth=0.6, zorder=5)
    # labels, placed without overlap
    gl = groups(verd)
    items = []
    for g in gl:
        p0, p1 = g["pos"][0], g["pos"][-1]
        s = sp[sp.position.between(p0, p1)]
        y = s[["d_mean_y10_A", "d_medoid_y10_A"]].max(axis=1).max()
        rng = f"{p0}-{p1}" if p1 > p0 else f"{p0}"
        txt = f"{rng}: {SHORT[g['verdict']]}" + (" (medoid only)" if g["medoid_only"] else "")
        items.append(((p0 + p1) / 2, y, txt, g["verdict"]))
    lv = place([(x, txt) for x, _, txt, _ in items], [0.62, 0.72, 0.82, 0.92])
    for (x, y, txt, v), h in zip(items, lv):
        ax.annotate(txt, xy=(x, y), xycoords="data", xytext=(x, h), textcoords=ax.get_xaxis_transform(),
                    ha="center", va="bottom", fontsize=7.5, color=F.INK,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=VERDICT_STYLE[v][0], lw=1.0),
                    arrowprops=dict(arrowstyle="-", color=F.INK2, lw=0.5))
    h, lab = ax.get_legend_handles_labels()
    for v in verd.verdict.unique():
        h.append(plt.Line2D([], [], marker="o", ls="", color=VERDICT_STYLE[v][0], markeredgecolor=F.INK,
                            markersize=6))
        lab.append(f"spike verdict: {VERDICT_STYLE[v][1]}")
    ax.legend(h, lab, loc="lower left", fontsize=7.5, ncol=3, bbox_to_anchor=(0.0, 1.035))
    vc = verd.verdict.value_counts()
    ax.set_title(f"(a) {len(verd)} positions outside 528-550 exceed 1.0 A: "
                 f"{vc.get('genuine loop difference', 0)} keep > 1.0 A after removing the 3 most influential entries, "
                 f"{vc.get('driven by <= 3 entries', 0)} rest on <= 3 entries; the restricted cluster-mean curve equals "
                 f"Y10 (max |diff| {rv['max_abs_diff_d_mean_A']:.1g} A)", fontsize=9.5, pad=62)
    ax.tick_params(labelbottom=True)

    # (b) observation fraction
    bands(ax_b, cons)
    ax_b.plot(prof.position, prof.obs_cl7 * 100, color=F.BLUE, lw=1.5, label="cl7")
    ax_b.plot(prof.position, prof.obs_cl4 * 100, color=F.ORANGE, lw=1.5, label="cl4")
    for p in verd.position:
        ax_b.axvline(p, color=F.INK2, lw=0.4, alpha=0.6, zorder=0)
    ax_b.set(ylabel="chains observing (%)", xlabel="ESR1 residue number (UniProt P03372)", ylim=(0, 105),
             xlim=(prof.position.min() - 1, prof.position.max() + 1))
    ax_b.legend(loc="lower left", fontsize=8)
    lo = sp.sort_values("obs_cl7").iloc[0]
    ax_b.set_title(f"(b) Share of chains observing each position (SV00 mask); thin lines = spike positions; lowest at a "
                   f"spike: {int(lo.position)}, cl7 {lo.obs_cl7:.0%} / cl4 {lo.obs_cl4:.0%}", fontsize=9.5)
    fig.suptitle("SV07 spike audit of the Y10 displacement profile (core v3 frame, S08 clusters cl7 vs cl4)",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    out = F.save_figure(fig, "SV07", "SV07_spike_profile", panels={"profile": ax, "observation": ax_b})
    for p in out:
        print(p.relative_to(ROOT))


if __name__ == "__main__":
    main()
