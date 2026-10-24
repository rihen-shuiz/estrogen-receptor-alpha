"""Tracker Y23: joint ESR1 + ESR2 C-alpha space and per-residue receptor differences.

ESR1: the 788 S03-superposed chains (aligned_coordinates.npz, core-v3 fit on the S03 reference). ESR2: the 86 X02
chains, read from results/lbd_extract/residue_mapping.csv (deposited coordinates, ESR1 numbering via the E03 map)
and fitted onto the S03 reference chain with esr_analysis.align_structures.fit_to_reference on the shared core-v3
positions (min_core_ca_for_alignment) - the same fit S03 used for ESR1. The S03 reference is the chain whose saved
transform is the identity (found, not assumed).
Joint matrix: 874 x 874 RMSD over the LBD positions observed in both chains (as S04, no refit; new file
results/ext/rmsd_joint.npy). Check: its ESR1 block must equal results/rmsd_entry.npy.
Joint map: metric MDS of the joint matrix (seed 42). Per position: distance between the mean C-alpha positions of
ESR1 cl7 and the ESR2 main basin (X03 cluster 2, peptide-bound), and of ESR1 cl4 and the ESR2 second basin (X03 5),
with the ESR1/ESR2 amino acids (results/numbering_map_esr1_esr2.csv), the draft selectivity residues
(data/selectivity_residues.csv, E01 MANUAL: draft, unconfirmed) and E10 consistently displaced positions marked.
Outputs: results/ext/rmsd_joint.npy, results/ext/joint_index.csv, results/ext/joint_embedding.csv,
results/ext/receptor_residue_diff.csv, figures/ext/Y23_joint_embedding.png, figures/ext/Y23_receptor_residue_diff.png.
Run from the repo root: PYTHONPATH=. python scripts/ext_joint_esr1_esr2.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.manifold import MDS

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from esr_analysis.align_structures import fit_to_reference
from ext_common import (AQUA, BLUE, CATEGORICAL, GREY, INK, INK2, ORANGE, RED, RES, ROOT, SEED, SETTINGS, VIOLET,
                        YELLOW, chain_table, coords, helix_bands, mpl_setup, rmsd, save_fig, write_numbers)
from ext_metric_check import masked_rmsd


def main():
    x1, pos = coords()
    d1, idx = rmsd()
    t = chain_table()
    tr = json.loads((ROOT / "results/alignment/alignment_transforms.json").read_text())
    ref_key = min(tr, key=lambda k: np.abs(np.array(tr[k]["rotation"]) - np.eye(3)).sum()
                  + np.abs(np.array(tr[k]["mobile_centroid"]) - np.array(tr[k]["reference_centroid"])).sum())
    ri = int(np.where(idx.chain_key == ref_key)[0][0])
    reference = {int(p): x1[ri, k] for k, p in enumerate(pos) if np.isfinite(x1[ri, k, 0])}
    core = json.loads((ROOT / "esr_analysis/resources/config/alignment_core_v3.json").read_text())
    core_pos = [p for a, b in core["core_ranges"] for p in range(a, b + 1)]
    keys2 = [l.strip() for l in (ROOT / "results/rmsd_esr2_index.txt").read_text().splitlines() if l.strip()]
    rm = pd.read_csv(ROOT / "results/lbd_extract/residue_mapping.csv", low_memory=False,
                     usecols=["chain_key", "esr1_position", "has_ca", "x_A", "y_A", "z_A"])
    rm = rm[rm.chain_key.isin(keys2) & rm.has_ca & rm.esr1_position.notna()]
    x2 = np.full((len(keys2), len(pos), 3), np.nan)
    fits = []
    pix = {int(p): k for k, p in enumerate(pos)}
    for i, k in enumerate(keys2):
        g = rm[rm.chain_key == k]
        c = {int(p): np.array([a, b, cc]) for p, a, b, cc in zip(g.esr1_position, g.x_A, g.y_A, g.z_A)}
        moved, info = fit_to_reference(c, reference, core_pos, min_core=SETTINGS["min_core_ca_for_alignment"])
        fits.append(dict(chain_key=k, status=info.get("status"), core_rmsd_A=info.get("core_alignment_rmsd_A"),
                         n_core=info.get("n_common_core_ca")))
        for p, xyz in moved.items():
            if p in pix:
                x2[i, pix[p]] = xyz
    fits = pd.DataFrame(fits)
    X = np.concatenate([x1, x2])
    dj, cnt = masked_rmsd(X, np.ones(len(pos), bool))
    blk = float(np.nanmax(np.abs(dj[:788, :788] - d1)))
    np.save(RES / "rmsd_joint.npy", dj)
    e2 = pd.read_csv(ROOT / "results/esr2_clusters.csv").set_index("chain_key").loc[keys2]
    jidx = pd.DataFrame({"row": np.arange(len(X)), "chain_key": list(idx.chain_key) + keys2,
                         "receptor": ["ESR1"] * 788 + ["ESR2"] * len(keys2),
                         "cluster": list(t.cluster.map(lambda c: f"ESR1 cl{c}")) + [f"ESR2 X03-{c}" for c in e2["cut_1.5"]]})
    jidx.to_csv(RES / "joint_index.csv", index=False)
    z = MDS(2, dissimilarity="precomputed", random_state=SEED, n_init=2, max_iter=300, normalized_stress="auto").fit_transform(dj)
    jidx["mds1"], jidx["mds2"] = z[:, 0], z[:, 1]
    jidx = jidx.merge(fits[["chain_key", "core_rmsd_A"]], on="chain_key", how="left")
    jidx.to_csv(RES / "joint_embedding.csv", index=False)
    # cross-receptor distances
    is2 = jidx.receptor.values == "ESR2"
    c2 = e2["cut_1.5"].values
    m7, m4 = t.cluster.values == 7, t.cluster.values == 4
    b2m, b2s = c2 == 2, c2 == 5
    blk_med = {
        "ESR1cl7-ESR2main": float(np.median(dj[:788][m7][:, 788:][:, b2m])),
        "ESR1cl4-ESR2main": float(np.median(dj[:788][m4][:, 788:][:, b2m])),
        "ESR1cl7-ESR2second": float(np.median(dj[:788][m7][:, 788:][:, b2s])),
        "ESR1cl4-ESR2second": float(np.median(dj[:788][m4][:, 788:][:, b2s])),
        "ESR1cl7-ESR1cl4": float(np.median(d1[np.ix_(m7, m4)])),
    }
    nm = pd.read_csv(ROOT / "results/numbering_map_esr1_esr2.csv").set_index("esr1_pos")
    sel = pd.read_csv(ROOT / "data/selectivity_residues.csv")
    cons = pd.read_csv(ROOT / "results/consistent_positions.csv").set_index("position")
    out = pd.DataFrame({"position": pos})
    out["esr1_aa"] = out.position.map(nm.esr1_aa)
    out["esr2_aa"] = out.position.map(nm.esr2_aa)
    out["esr2_pos"] = out.position.map(nm.esr2_pos)
    out["same_residue"] = out.esr1_aa == out.esr2_aa
    mean = lambda a: np.nanmean(a, 0)
    out["d_mean_A"] = np.linalg.norm(mean(x1[m7]) - mean(x2[b2m]), axis=1)
    out["d_mean_antag_A"] = np.linalg.norm(mean(x1[m4]) - mean(x2[b2s]), axis=1)
    out["obs_esr2_main"] = np.isfinite(x2[b2m][..., 0]).mean(0)
    out["draft_selectivity"] = out.position.isin(sel.esr1_pos)
    out["e10_consistent"] = out.position.map(cons.consistent).fillna(False).astype(bool)
    out.to_csv(RES / "receptor_residue_diff.csv", index=False)

    plt = mpl_setup()
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axs[0]
    for lab, m, col, mk, s in (("ESR1 cl7", ~is2 & np.r_[m7, np.zeros(len(keys2), bool)], BLUE, "o", 8),
                               ("ESR1 cl4", ~is2 & np.r_[m4, np.zeros(len(keys2), bool)], ORANGE, "o", 8),
                               ("ESR1 minor", ~is2 & ~np.r_[m7 | m4, np.ones(len(keys2), bool)], AQUA, "^", 14),
                               ("ESR2 main basin (X03-2)", np.r_[np.zeros(788, bool), b2m], VIOLET, "D", 26),
                               ("ESR2 second basin (X03-5)", np.r_[np.zeros(788, bool), b2s], YELLOW, "D", 26),
                               ("ESR2 other", np.r_[np.zeros(788, bool), ~(b2m | b2s)], GREY, "D", 26)):
        ax.scatter(z[m, 0], z[m, 1], s=s, c=col, marker=mk, lw=0.4, edgecolor="white" if mk == "o" else INK,
                   label=f"{lab} (n={m.sum()})", zorder=3 if "ESR2" in lab else 2)
    ax.legend(fontsize=8)
    ax.set(xlabel="joint MDS 1 (A)", ylabel="joint MDS 2 (A)",
           title=f"Joint map of 874 chains (ESR2 fitted on the ESR1 reference {ref_key} by core v3)")
    ax = axs[1]
    keys_ = list(blk_med)
    ax.barh(keys_[::-1], [blk_med[k] for k in keys_[::-1]], color=[INK2 if "ESR1cl7-ESR1cl4" == k else VIOLET for k in keys_[::-1]])
    for y, k in enumerate(keys_[::-1]):
        ax.text(blk_med[k] + 0.05, y, f"{blk_med[k]:.2f} A", va="center", fontsize=8, color=INK2)
    ax.set(xlabel="median C-alpha RMSD between the groups (A)", title="Receptor vs state: which distance is larger?")
    fig.suptitle(f"ESR2 agonist-state chains sit next to ESR1 cl7 ({blk_med['ESR1cl7-ESR2main']:.2f} A) and the ESR2 second basin "
                 f"next to ESR1 cl4 ({blk_med['ESR1cl4-ESR2second']:.2f} A): the H12 state matters more than the receptor",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y23_joint_embedding", section="I. ESR2 and ESR1 vs ESR2", task="Y23",
             shows="metric MDS of the joint 874 x 874 ESR1+ESR2 RMSD matrix; median RMSD between receptor/state groups",
             takeaway=(f"Same-state cross-receptor distances ({blk_med['ESR1cl7-ESR2main']:.2f} / {blk_med['ESR1cl4-ESR2second']:.2f} A) "
                       f"are far below the ESR1 cl7-cl4 distance ({blk_med['ESR1cl7-ESR1cl4']:.2f} A)."),
             source="results/ext/rmsd_joint.npy, results/ext/joint_embedding.csv", script="scripts/ext_joint_esr1_esr2.py")
    fig, axs = plt.subplots(2, 1, figsize=(16, 7.5), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    ax = axs[0]
    helix_bands(ax)
    ax.plot(out.position, out.d_mean_A, color=VIOLET, lw=1.6, label="ESR1 cl7 vs ESR2 main basin (agonist-like)")
    ax.plot(out.position, out.d_mean_antag_A, color=YELLOW, lw=1.2, label="ESR1 cl4 vs ESR2 second basin")
    diff = out[~out.same_residue & out.esr2_aa.notna()]
    ax.scatter(diff.position, diff.d_mean_A, s=14, c=INK, zorder=4, label=f"sequence differs (n={len(diff)})")
    for r in out[out.draft_selectivity].itertuples():
        ax.annotate(f"{r.esr1_aa}{r.position}{r.esr2_aa} (draft E01)", (r.position, r.d_mean_A), xytext=(0, 28),
                    textcoords="offset points", ha="center", fontsize=8, color=RED,
                    arrowprops=dict(arrowstyle="-", color=RED, lw=0.7))
    ax.set_yscale("symlog", linthresh=1)
    ax.set(ylabel="distance between mean C-alpha (A, symlog)")
    ax.legend(fontsize=8, loc="upper left")
    ax = axs[1]
    helix_bands(ax, label=False)
    ax.bar(out.position, (~out.same_residue & out.esr2_aa.notna()).astype(int), color=INK2, width=0.9, label="ESR1 != ESR2 residue")
    ax.bar(out.position, out.e10_consistent.astype(int) * -1, color=CATEGORICAL[4], width=0.9, label="E10 consistently displaced")
    ax.set(ylabel="flag", xlabel="ESR1 residue", yticks=[-1, 0, 1], yticklabels=["E10", "", "seq."], xlim=(304, 551))
    ax.legend(fontsize=8, loc="lower left", ncol=2)
    lbd = out[out.position < 528]
    fig.suptitle(f"Per-residue ESR1 vs ESR2 (same H12 state): outside 528-550 the median shift is {lbd.d_mean_A.median():.2f} A; "
                 f"the two draft pocket selectivity positions (384, 421) move < {out[out.draft_selectivity].d_mean_A.max():.1f} A",
                 x=0.01, ha="left", fontweight="bold", fontsize=11)
    fig.tight_layout()
    save_fig(fig, "Y23_receptor_residue_diff", section="I. ESR2 and ESR1 vs ESR2", task="Y23",
             shows="per-position distance between ESR1 and ESR2 mean C-alpha traces in the same H12 state, sequence differences, draft selectivity residues, E10 positions",
             takeaway=("Backbone differences between the receptors are small in the core and concentrate in loops; the pocket "
                       "selectivity residues differ in side chain (L384M, M421I), not in backbone position."),
             source="results/ext/receptor_residue_diff.csv", script="scripts/ext_joint_esr1_esr2.py")
    write_numbers("Y23", dict(reference=ref_key, esr1_block_maxdiff=blk, n_esr2_fit_ok=int((fits.status == "ok").sum()),
                              esr2_core_rmsd_median=float(fits.core_rmsd_A.median()), block_medians=blk_med,
                              median_shift_lt528=float(lbd.d_mean_A.median()),
                              selectivity_shift=out[out.draft_selectivity].set_index("position").d_mean_A.round(3).to_dict(),
                              n_seq_diff_lbd=int(len(diff))))
    print(ref_key, blk, fits.status.value_counts().to_dict(), blk_med, lbd.d_mean_A.median())


if __name__ == "__main__":
    main()
