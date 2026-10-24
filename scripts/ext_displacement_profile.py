"""Tracker Y10: where along the sequence do cl7 and cl4 differ? Per-position C-alpha displacement profile.

Input: S03 superposed C-alpha coordinates (results/alignment/aligned_coordinates.npz; core-v3 fit, no refit), S08
1.5 A clusters, H12 labels. Per ESR1 position 305-550:
  d_mean_A       distance between the cl7 and cl4 mean positions (nan-mean over chains observing the position)
  d_medoid_A     distance between the two 'complete' medoids: the most central chain of each cluster among chains
                 with a complete H12 (the plain S12 medoid of cl4, 2QE4|B|1, has no modelled H12; see Y42)
  sd_cl7_A, sd_cl4_A   RMS deviation of the cluster's chains from the cluster mean position (spread)
  obs_cl7, obs_cl4     fraction of the cluster's chains that observe the position
  share_sq       share of sum(d_mean^2) contributed by the position (where the between-cluster RMSD comes from)
Outputs: results/ext/cluster_displacement_profile.csv, figures/ext/Y10_displacement_profile.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_displacement_profile.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (BLUE, GREY, H12, INK2, ORANGE, RES, VIOLET, chain_table, coords, helix_bands, mpl_setup, rmsd,
                        save_fig, write_numbers)


def complete_medoid(d, t, c):
    m = np.where((t.cluster.values == c) & (t.h12_label.values == "complete"))[0]
    return m[np.argmin(d[np.ix_(m, m)].sum(1))]


def main():
    d, idx = rmsd()
    t = chain_table()
    x, pos = coords()
    out = pd.DataFrame({"position": pos})
    mu, sd, ob = {}, {}, {}
    for c in (7, 4):
        xc = x[t.cluster.values == c]
        mu[c] = np.nanmean(xc, 0)
        dev = np.sqrt(np.nansum((xc - mu[c]) ** 2, 2))
        dev[~np.isfinite(xc[..., 0])] = np.nan
        sd[c] = np.sqrt(np.nanmean(dev ** 2, 0))
        ob[c] = np.isfinite(xc[..., 0]).mean(0)
    m7, m4 = complete_medoid(d, t, 7), complete_medoid(d, t, 4)
    out["d_mean_A"] = np.linalg.norm(mu[7] - mu[4], axis=1)
    out["d_medoid_A"] = np.linalg.norm(x[m7] - x[m4], axis=1)
    out["sd_cl7_A"], out["sd_cl4_A"] = sd[7], sd[4]
    out["obs_cl7"], out["obs_cl4"] = ob[7], ob[4]
    sq = out.d_mean_A ** 2
    out["share_sq"] = sq / sq.sum()
    out["medoid_cl7"], out["medoid_cl4"] = idx.chain_key[m7], idx.chain_key[m4]
    out.to_csv(RES / "cluster_displacement_profile.csv", index=False)
    h12 = (out.position >= H12[0]) & (out.position <= H12[1])
    tail = out.position >= 528
    share_h12 = out.share_sq[h12].sum()
    share_tail = out.share_sq[tail].sum()
    top = out.sort_values("d_mean_A", ascending=False).head(5)

    plt = mpl_setup()
    fig, axs = plt.subplots(3, 1, figsize=(16, 9.5), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.3, 1]})
    ax = axs[0]
    helix_bands(ax)
    ax.plot(out.position, out.d_mean_A, color=VIOLET, lw=1.8, label="cluster means (all chains)")
    ax.plot(out.position, out.d_medoid_A, color=INK2, lw=1, ls="--",
            label=f"complete-H12 medoids {idx.chain_key[m7]} vs {idx.chain_key[m4]}")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("cl7 vs cl4 C-alpha distance (A, symlog)")
    loops = [(330, 340), (460, 472)]
    for lo, hi in loops:
        seg = out[(out.position >= lo) & (out.position <= hi)]
        r = seg.loc[seg.d_mean_A.idxmax()]
        ax.annotate(f"{int(r.position)}: {r.d_mean_A:.1f} A", (r.position, r.d_mean_A), xytext=(8, 8),
                    textcoords="offset points", fontsize=8, color=INK2)
    hh = out[h12].loc[out[h12].d_mean_A.idxmax()]
    ax.annotate(f"H12 max {int(hh.position)}: {hh.d_mean_A:.1f} A", (hh.position, hh.d_mean_A), xytext=(-150, 0),
                textcoords="offset points", fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title(f"Positions 528-550 carry {share_tail:.0%} of the squared cl7-cl4 difference (H12 538-547 alone "
                 f"{share_h12:.0%}); the rest of the LBD differs by < {out[~tail].d_mean_A.max():.1f} A", fontsize=10)
    ax = axs[1]
    helix_bands(ax, label=False)
    ax.plot(out.position, out.sd_cl7_A, color=BLUE, lw=1.5, label="cl7 spread (RMS dev. from mean)")
    ax.plot(out.position, out.sd_cl4_A, color=ORANGE, lw=1.5, label="cl4 spread")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("within-cluster spread (A)")
    ax.legend(loc="upper left", fontsize=8)
    ax = axs[2]
    helix_bands(ax, label=False)
    ax.plot(out.position, out.obs_cl7 * 100, color=BLUE, lw=1.5, label="cl7")
    ax.plot(out.position, out.obs_cl4 * 100, color=ORANGE, lw=1.5, label="cl4")
    ax.set(ylabel="chains observing (%)", xlabel="ESR1 residue number (UniProt P03372)", ylim=(0, 105))
    ax.legend(loc="lower left", fontsize=8)
    fig.suptitle("Why two clusters: the difference is concentrated in H12 and the loop before it (blue bands = helices, "
                 "red band = H12)", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y10_displacement_profile", section="C. Why two clusters", task="Y10",
             shows="per-residue distance between cl7 and cl4 mean C-alpha positions and medoids, within-cluster spread, observation rate",
             takeaway=(f"Largest shifts at {', '.join(str(p) for p in top.position.head(3))} "
                       f"({top.d_mean_A.iloc[0]:.1f} A max); positions 528-550 hold {share_tail:.0%} of the squared "
                       f"difference; cl4 chains observe H12 less often ({out[h12].obs_cl4.mean():.0%} vs "
                       f"{out[h12].obs_cl7.mean():.0%})."),
             source="results/ext/cluster_displacement_profile.csv", script="scripts/ext_displacement_profile.py")
    write_numbers("Y10", dict(share_528_550=float(share_tail), share_h12=float(share_h12),
                              max_d_outside_528=float(out[~tail].d_mean_A.max()),
                              top_positions=top.position.tolist(), top_d=top.d_mean_A.round(2).tolist(),
                              obs_h12_cl7=float(out[h12].obs_cl7.mean()), obs_h12_cl4=float(out[h12].obs_cl4.mean()),
                              medoid_cl7=idx.chain_key[m7], medoid_cl4=idx.chain_key[m4]))
    print(top.to_string())
    print(share_tail, share_h12)


if __name__ == "__main__":
    main()
