"""Tracker Y11: per-chain structural features that might explain the cl7 / cl4 split.

Input: S03 superposed C-alpha coordinates (core-v3 fit, common frame), results/helix_table.csv helix ranges,
results/pockets.csv (X09 P2Rank top pocket), results/ext/ligand_descriptors.csv (Y09), chain_table().
Per chain (blank where the needed C-alpha atoms are not modelled; nothing dropped):
  h12_n_ca            C-alpha atoms observed in H12 (538-547)
  h12_centroid_x/y/z  H12 centroid in the common frame (>= 4 C-alpha observed)
  h12_shift_A         distance of the H12 centroid from the cl7 (agonist-like) mean H12 centroid
  h12_to_h3_A, h12_to_h5_A, h12_to_h11_A   H12 centroid to H3 (338-362), H4/5 (371-396), H10/11 (496-531) centroids
  h12_h3_angle_deg, h12_h11_angle_deg      angle between helix axes (first principal axis, oriented N->C)
  k362_e542_A         Lys362 - Glu542 C-alpha distance (AF-2 charge-clamp residues; C-alpha proxy)
  l540_pocket_A       Leu540 C-alpha to the centroid of pocket residues 353/394/524 (settings pocket_positions)
  loop_533_537_n_ca   observed C-alpha atoms in the H11-H12 loop
  h2h3_loop_shift_A   mean distance of 330-340 C-alpha atoms from the cl7 mean (second region of difference, Y10)
  pocket_sas_points, pocket_score (X09), ligand_heavy_atoms, ligand_mw, ligand_basic_amine (Y09), lbd_ca_observed
Contact maps: C-alpha pairs within 8 A, frequency per cluster over chains observing both positions; the difference
map cl7 - cl4 is saved (results/ext/contact_freq_diff.npy).
Outputs: results/ext/structural_features.csv, results/ext/contact_freq_diff.npy, figures/ext/Y11_features_by_cluster.png,
figures/ext/Y11_contact_map_diff.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_structural_features.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ext_common import (AQUA, BLUE, GREY, INK, INK2, ORANGE, RES, ROOT, SETTINGS, bh, chain_table, coords,
                        helix_bands, helices, mpl_setup, save_fig, write_numbers, HELIX_NAMES)

RANGES = {"H3": (338, 362), "H5": (371, 396), "H11": (496, 531), "H12": (538, 547)}


def idx_of(pos, lo, hi):
    return np.where((pos >= lo) & (pos <= hi))[0]


def axis(pts):
    pts = pts[np.isfinite(pts[:, 0])]
    if len(pts) < 4:
        return None
    c = pts - pts.mean(0)
    v = np.linalg.svd(c, full_matrices=False)[2][0]
    if np.dot(pts[-1] - pts[0], v) < 0:
        v = -v
    return v


def centroid(pts, min_n=4):
    pts = pts[np.isfinite(pts[:, 0])]
    return pts.mean(0) if len(pts) >= min_n else np.full(3, np.nan)


def angle(a, b):
    if a is None or b is None:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))))


def main():
    t = chain_table()
    x, pos = coords()
    p = {k: idx_of(pos, *v) for k, v in RANGES.items()}
    ip = {r: int(np.where(pos == r)[0][0]) for r in (362, 542, 540)}
    pocket_idx = [int(np.where(pos == r)[0][0]) for r in SETTINGS["pocket_positions"]]
    cl7 = t.cluster.values == 7
    h12c = np.array([centroid(x[i, p["H12"]]) for i in range(len(x))])
    ref_h12 = np.nanmean(h12c[cl7], 0)
    loop = idx_of(pos, 330, 340)
    loop_ref = np.nanmean(x[cl7][:, loop], 0)
    rows = []
    for i in range(len(x)):
        xi = x[i]
        cen = {k: centroid(xi[v]) for k, v in p.items()}
        ax = {k: axis(xi[v]) for k, v in p.items()}
        pc = centroid(xi[pocket_idx], 2)
        lsh = np.linalg.norm(xi[loop] - loop_ref, axis=1)
        rows.append(dict(
            chain_key=t.chain_key[i], h12_n_ca=int(np.isfinite(xi[p["H12"], 0]).sum()),
            h12_centroid_x=h12c[i, 0], h12_centroid_y=h12c[i, 1], h12_centroid_z=h12c[i, 2],
            h12_shift_A=np.linalg.norm(h12c[i] - ref_h12),
            h12_to_h3_A=np.linalg.norm(cen["H12"] - cen["H3"]), h12_to_h5_A=np.linalg.norm(cen["H12"] - cen["H5"]),
            h12_to_h11_A=np.linalg.norm(cen["H12"] - cen["H11"]),
            h12_h3_angle_deg=angle(ax["H12"], ax["H3"]), h12_h11_angle_deg=angle(ax["H12"], ax["H11"]),
            k362_e542_A=np.linalg.norm(xi[ip[362]] - xi[ip[542]]),
            l540_pocket_A=np.linalg.norm(xi[ip[540]] - pc),
            loop_533_537_n_ca=int(np.isfinite(xi[idx_of(pos, 533, 537), 0]).sum()),
            h2h3_loop_shift_A=np.nanmean(lsh) if np.isfinite(lsh).sum() >= 3 else np.nan,
            lbd_ca_observed=int(np.isfinite(xi[:, 0]).sum())))
    f = pd.DataFrame(rows)
    f = f.merge(t[["chain_key", "entry_id", "cluster", "ligand_id", "sas_points", "pocket_score", "h12_label", "tier",
                   "resolution_A", "r_free", "mean_b_A2"]].rename(columns={"sas_points": "pocket_sas_points"}),
                on="chain_key", how="left", validate="one_to_one")
    lig = pd.read_csv(RES / "ligand_descriptors.csv")[["ligand_id", "heavy_atoms", "mw", "basic_amine", "rotb"]]
    f = f.merge(lig.rename(columns={"heavy_atoms": "ligand_heavy_atoms", "mw": "ligand_mw",
                                    "basic_amine": "ligand_basic_amine", "rotb": "ligand_rotb"}),
                on="ligand_id", how="left")
    assert len(f) == len(t)
    f.to_csv(RES / "structural_features.csv", index=False)

    # contact maps
    obs = np.isfinite(x[..., 0])
    freq = {}
    for c, m in (("cl7", cl7), ("cl4", t.cluster.values == 4)):
        xc, oc = x[m], obs[m]
        cnt = np.zeros((len(pos), len(pos)))
        both = np.zeros((len(pos), len(pos)))
        for k in range(len(xc)):
            xx = np.where(oc[k][:, None], xc[k], 0)
            dd = np.sqrt(((xx[:, None] - xx[None]) ** 2).sum(-1))
            o2 = np.outer(oc[k], oc[k])
            cnt += (dd < 8) & o2
            both += o2
        with np.errstate(invalid="ignore", divide="ignore"):
            freq[c] = np.where(both >= 10, cnt / both, np.nan)
    diff = freq["cl7"] - freq["cl4"]
    np.save(RES / "contact_freq_diff.npy", diff)

    # tests per feature (chain level, descriptive) cl7 vs cl4
    feats = ["h12_shift_A", "h12_to_h3_A", "h12_to_h5_A", "h12_to_h11_A", "h12_h3_angle_deg", "h12_h11_angle_deg",
             "k362_e542_A", "l540_pocket_A", "h12_n_ca", "loop_533_537_n_ca", "h2h3_loop_shift_A", "pocket_sas_points",
             "ligand_heavy_atoms"]
    units = {"h12_shift_A": "A", "h12_to_h3_A": "A", "h12_to_h5_A": "A", "h12_to_h11_A": "A", "h12_h3_angle_deg": "deg",
             "h12_h11_angle_deg": "deg", "k362_e542_A": "A", "l540_pocket_A": "A", "h12_n_ca": "C-alpha",
             "loop_533_537_n_ca": "C-alpha", "h2h3_loop_shift_A": "A", "pocket_sas_points": "points",
             "ligand_heavy_atoms": "atoms"}
    tr = []
    for v in feats:
        a = f.loc[f.cluster == 7, v].dropna()
        b = f.loc[f.cluster == 4, v].dropna()
        u, pv = mannwhitneyu(a, b)
        tr.append(dict(feature=v, n_cl7=len(a), n_cl4=len(b), median_cl7=a.median(), median_cl4=b.median(),
                       rank_biserial=1 - 2 * u / (len(a) * len(b)), p_value=pv))
    tr = pd.DataFrame(tr)
    tr["q_value"] = bh(tr.p_value)
    tr.to_csv(RES / "structural_feature_tests.csv", index=False)

    plt = mpl_setup()
    show = feats[:12]
    fig, axs = plt.subplots(3, 4, figsize=(17, 10.5))
    rng = np.random.default_rng(42)
    for axx, v in zip(axs.ravel(), show):
        groups = [(7, BLUE, "cl7"), (4, ORANGE, "cl4"), (-1, AQUA, "minor")]
        for k, (c, col, name) in enumerate(groups):
            vals = f.loc[(f.cluster == c) if c > 0 else ~f.cluster.isin((7, 4)), v].dropna().values
            if len(vals) > 3:
                parts = axx.violinplot([vals], positions=[k], widths=0.8, showextrema=False)
                for pc in parts["bodies"]:
                    pc.set_facecolor(col); pc.set_alpha(0.3); pc.set_edgecolor("none")
            axx.scatter(k + rng.uniform(-0.18, 0.18, len(vals)), vals, s=3, c=col, lw=0, alpha=0.6)
            if len(vals):
                axx.hlines(np.median(vals), k - 0.3, k + 0.3, color=INK, lw=1.4)
        r = tr.set_index("feature").loc[v]
        axx.set_xticks([0, 1, 2], [f"cl7\nn={int(r.n_cl7)}", f"cl4\nn={int(r.n_cl4)}", "minor"])
        axx.set_ylabel(f"{v} ({units[v]})", fontsize=8)
        axx.set_title(f"{v}: median {r.median_cl7:.1f} vs {r.median_cl4:.1f}, |r|={abs(r.rank_biserial):.2f}",
                      fontsize=8.5)
    s = tr.set_index("feature")
    fig.suptitle(f"Structural features by cluster: H12 sits {s.loc['h12_shift_A', 'median_cl4']:.0f} A from its cl7 "
                 f"position in cl4 and K362-E542 shortens from {s.loc['k362_e542_A', 'median_cl7']:.0f} to "
                 f"{s.loc['k362_e542_A', 'median_cl4']:.0f} A (H12 lies in the coactivator groove)", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y11_features_by_cluster", section="C. Why two clusters", task="Y11",
             shows="12 per-chain structural features (H12 geometry, helix angles, charge clamp, pocket and ligand size) by cluster",
             takeaway=(f"H12-H3 angle {s.loc['h12_h3_angle_deg', 'median_cl7']:.0f} vs "
                       f"{s.loc['h12_h3_angle_deg', 'median_cl4']:.0f} deg, H12-H11 distance "
                       f"{s.loc['h12_to_h11_A', 'median_cl7']:.1f} vs {s.loc['h12_to_h11_A', 'median_cl4']:.1f} A; "
                       f"pocket size and the H2-H3 loop differ less (|r| {abs(s.loc['pocket_sas_points', 'rank_biserial']):.2f}, "
                       f"{abs(s.loc['h2h3_loop_shift_A', 'rank_biserial']):.2f})."),
             source="results/ext/structural_features.csv, results/ext/structural_feature_tests.csv",
             script="scripts/ext_structural_features.py")

    fig, axs = plt.subplots(1, 3, figsize=(19, 6.3), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    ext = [pos[0] - 0.5, pos[-1] + 0.5, pos[-1] + 0.5, pos[0] - 0.5]
    for a_, key in zip(axs[:2], ("cl7", "cl4")):
        im = a_.imshow(freq[key], cmap="esr_seq", vmin=0, vmax=1, extent=ext, interpolation="nearest")
        a_.set_title(f"{key}: C-alpha contact frequency (< 8 A)")
        a_.set_xlabel("ESR1 residue"); a_.set_ylabel("ESR1 residue")
    fig.colorbar(im, ax=axs[1], fraction=0.04, label="fraction of chains in contact")
    a_ = axs[2]
    im = a_.imshow(diff, cmap="esr_div", vmin=-1, vmax=1, extent=ext, interpolation="nearest")
    for h in helices().itertuples():
        if h.helix == "AC4":
            a_.axvspan(h.start - 0.5, h.end + 0.5, color="none", ec=INK, lw=0.6)
    a_.set_title("difference cl7 - cl4 (blue = cl4 only, red = cl7 only)")
    a_.set_xlabel("ESR1 residue"); a_.set_ylabel("ESR1 residue")
    fig.colorbar(im, ax=a_, fraction=0.04, label="contact frequency difference")
    big = np.nan_to_num(np.abs(diff)) > 0.5
    big_h12 = big[np.ix_(p["H12"], np.arange(len(pos)))].sum()
    fig.suptitle(f"Contact maps: {int(big.sum() // 2)} residue pairs change contact frequency by > 0.5 between the "
                 f"clusters; {int(big_h12)} of them involve H12", x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y11_contact_map_diff", section="C. Why two clusters", task="Y11",
             shows="C-alpha contact-frequency maps of cl7 and cl4 and their difference",
             takeaway=f"{int(big.sum() // 2)} residue pairs differ by > 0.5 in contact frequency; {int(big_h12)} involve H12 (H12 packs on H3/H11 in cl7 and on H3/H5 in the AF-2 groove in cl4).",
             source="results/ext/contact_freq_diff.npy, results/ext/structural_features.csv",
             script="scripts/ext_structural_features.py")
    write_numbers("Y11", dict(tests=tr.set_index("feature").round(4).to_dict("index"), n_pairs_diff05=int(big.sum() // 2),
                              n_pairs_diff05_h12=int(big_h12)))
    print(tr.round(3).to_string())


if __name__ == "__main__":
    main()
