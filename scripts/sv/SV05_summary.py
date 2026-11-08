"""SV05 (RUN_SV1028): write results/sv/SV05/env.txt and SUMMARY.md (rule 1.7) from the saved SV05 outputs only.
Run after SV05_core_v4.py and SV05_figures.py:
.venv/bin/python -B -I scripts/sv/SV05_summary.py
"""
import json
import platform
import sys
from datetime import datetime
from importlib import metadata
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV05"
CUTS = ["1.0", "1.5", "2.0", "2.5"]


def size(path):
    n = path.stat().st_size
    return f"{n / 1e6:.2f} MB" if n >= 1e6 else f"{n / 1e3:.1f} kB"


def env():
    pk = ["numpy", "pandas", "scipy", "scikit-learn", "matplotlib"]
    lines = [f"# SV05 environment, {datetime.now().astimezone().strftime('%Y-%m-%dT%H:%M%z')}",
             f"python {platform.python_version()} at {sys.executable}", f"platform {platform.platform()}"]
    lines += [f"{p} {metadata.version(p)}" for p in pk]
    lines += ["code path reused unchanged: scripts/align_selected_chains.fit_all, esr_analysis.calculate_rmsd.compute_distances,",
              "  esr_analysis._task_io.aligned_caches, esr_analysis.cluster_structures.structural_groups (as scripts/core_sensitivity.py)",
              "nothing installed, upgraded or downgraded by SV05", "network: none"]
    (OUT / "env.txt").write_text("\n".join(lines) + "\n")


def rng(s):
    return s if isinstance(s, str) and s else "-"


def main():
    env()
    n = json.loads((OUT / "numbers_core_v4.json").read_text())
    core = json.loads((OUT / "core_v4.json").read_text())
    sp = pd.read_csv(OUT / "spread.csv")
    ari = pd.read_csv(OUT / "ari.csv")
    comp = pd.read_csv(OUT / "core_comparison.csv")
    cl = pd.read_csv(OUT / "clusters_v4.csv")
    it = pd.read_csv(OUT / "iterations.csv")
    prof = pd.read_csv(OUT / "displacement_profile_v4.csv")
    pv4, pv3 = n["profile_v4"], n["profile_v3_y10"]

    moved = {}
    for c in CUTS:
        t = pd.crosstab(cl[f"cut_{c}"], cl[f"s08_cut_{c}"])
        best = t.idxmax(axis=1)
        m = cl[f"cut_{c}"].map(best) != cl[f"s08_cut_{c}"]
        moved[c] = cl.loc[m, ["chain_key", f"s08_cut_{c}", f"cut_{c}"]].assign(to=cl.loc[m, f"cut_{c}"].map(best))
    s15 = cl["s08_cut_1.5"]

    def name(c, lab):
        m = cl[f"s08_cut_{c}"] == lab
        return f"{lab} ({int(m.sum())} chains, 1.5 A label of most: cl{int(s15[m].mode().iloc[0])})"
    # 1.0 A regrouping: where do the large (>= 10 chains) v4 clusters come from (S08 1.0 A and 1.5 A labels)
    t10 = pd.crosstab(cl["cut_1.0"], cl["s08_cut_1.0"])
    big = t10[t10.sum(axis=1) >= 10]
    within = {int(k): sorted({int(x) for x in s15[cl["cut_1.0"] == k]}) for k in big.index}

    cand = sp[sp.candidate]
    cand_helices = [h for h in dict.fromkeys(cand.consensus_helix)]
    new_helices = comp[(comp.v3_n == 0) & (comp.v4_n > 0)]
    new_names = " and ".join(new_helices.helix)
    out_cand = cand[~cand.in_core_v4]
    v3n = int(sp.in_core_v3.sum())

    L = []
    L.append(
        "Not done / substituted: (1) interpretations fixed before the run (the step names no formula): 'mean shift' = mean "
        "|x_new - x_old| over all observed (chain, candidate) C-alpha between two iterations; each iteration Kabsch-fits every "
        "chain on its observed candidates onto the nan-mean candidate structure and then re-anchors the whole ensemble on the "
        "starting core-v3-frame mean (one common rigid motion, which leaves every spread unchanged); 'spread' = RMS deviation of "
        "the observing chains from the mean position (the Y10 sd formula); ties in the 123 lowest-spread set broken by position. "
        "(2) The Y10-style profile re-implements scripts/ext_displacement_profile.py (same formulas) instead of running it, "
        "because that script writes into results/ext and figures/ext; on the v3 inputs it reproduces "
        f"results/ext/cluster_displacement_profile.csv to {n['y10_reproduction_max_abs_diff']:.1e} with the same medoids. Clusters "
        "there are the S08 1.5 A cl7/cl4 labels; the complete-H12 medoid rule is applied on the v4 matrix. (3) The 123 "
        "lowest-spread set is reported as a set (core_v4.json, spread.csv, core_comparison.csv, Jaccard); no matrix or "
        "clustering was computed for it (step 5 names core v4 only). (4) 'Which helices differ' groups positions by the "
        "helix-table ranges in results/sv/SV02/helix_boundaries.csv (table_start-table_end), because v3 includes helix ends "
        "outside the SV02 consensus helices. (5) Supporting files not named in the run file: rmsd_v4_index.csv, clusters_v4.csv, "
        "cluster_sizes.csv, core_comparison.csv, core_sensitivity_v4.csv, displacement_profile_v4.csv, iterations.csv, "
        "numbers_core_v4.json, data/sv/SV05/. The superposition core of SV03/SV04/SV07 stays v3; nothing existing was "
        "overwritten.")
    L.append("")
    L.append("# SV05 - Stable core v4 from complete helices, compared with core v3")
    L.append("")
    r15 = ari[ari.cut_A == 1.5].iloc[0]
    L.append(f"**Verdict ([DEFAULT] step 3):** core v4 = **{n['n_core_v4']} positions** with spread <= 0.5 A "
             f"(>= 80, not BLOCKED), from {n['n_candidates']} candidates. Jaccard with v3 ({v3n} positions) = "
             f"**{n['jaccard_v4_v3']:.2f}** (lowest-spread 123: {n['jaccard_lowest123_v3']:.2f}). v4 adds {new_names} and drops "
             f"{int((sp.in_core_v3 & ~sp.in_core_v4).sum())} v3 positions at helix ends. Recomputed with v4 through the K09 code path: pairwise RMSD Spearman vs the primary matrix "
             f"**{n['spearman_pairwise_v4_vs_v3']:.4f}**, median |dRMSD| {n['median_abs_diff_A']:.3f} A (max "
             f"{n['max_abs_diff_A']:.2f} A); ARI against S08 at 1.0/1.5/2.0/2.5 A = "
             + "/".join(f"{v:.3f}" for v in ari.ari_v4_vs_s08)
             + f"; clusters {'/'.join(str(v) for v in ari.n_clusters_v4)} (v4) vs "
             f"{'/'.join(str(v) for v in ari.n_clusters_v3_s08)} (S08). At 1.5 A {int(r15.chains_changed_best_match)} chains "
             f"change cluster; the two main clusters keep {r15.sizes_v4.split(';')[0]}/{r15.sizes_v4.split(';')[1]} chains "
             f"(S08 {r15.sizes_v3_s08.split(';')[0]}/{r15.sizes_v3_s08.split(';')[1]}). The Y10-style profile on v4 puts "
             f"{pv4['share_528_550']:.1%} of the squared cl7-cl4 difference in 528-550 (v3 {pv3['share_528_550']:.1%}). "
             "Adopting v4 or keeping v3 is a team decision (run file FINAL); this step only compares.")
    L.append("")
    L.append("## Step 1 - candidates (spread.csv)")
    L.append("")
    L.append(f"- SV02 global consensus helices: {n['n_consensus_positions']} positions (results/sv/SV02/consensus_global.csv); "
             f"in 528-550: {n['n_consensus_in_528_550']} (all of H12 538-544); observed in < 98 % of the 788 chains: "
             f"{n['n_consensus_below_98pct']} ({', '.join(str(p) for p in n['consensus_below_98pct_positions'])}, "
             f"{sp.obs_frac[sp.position.isin(n['consensus_below_98pct_positions'])].iloc[0]:.1%}).")
    L.append(f"- Candidates: **{n['n_candidates']}** positions in {len(cand_helices)} helices ({', '.join(cand_helices)}).")
    L.append("")
    L.append("## Step 2 - iterative superposition (iterations.csv)")
    L.append("")
    L.append(f"- Start: core v3 frame of the SV00 gate (data/sv/gate/aligned_coordinates.npz, reference 7NFB|B|1). "
             f"Converged after **{n['iterations']} iterations**: mean shift "
             + ", ".join(f"{v:.2g}" for v in it.mean_shift_A) + " A (stop < 0.01 A); mean per-chain RMS to the mean "
             f"structure on the candidates {it.mean_rms_to_mean_A.iloc[-1]:.3f} A.")
    L.append(f"- Candidate spread: min {n['candidate_spread_min_A']:.2f}, median {n['candidate_spread_median_A']:.2f}, max "
             f"{n['candidate_spread_max_A']:.2f} A (in the v3 frame before iterating: median "
             f"{cand.spread_v3_frame_A.median():.2f}, max {cand.spread_v3_frame_A.max():.2f} A).")
    L.append("")
    L.append("## Step 3 - core v4 (core_v4.json)")
    L.append("")
    L.append(f"- Core v4 = {n['n_core_v4']} positions: " + ", ".join(f"{a}-{b}" for a, b in core["core_ranges"]) + ".")
    L.append(f"- Candidates above 0.5 A (left out): " + ", ".join(
        f"{int(r.position)} ({r.consensus_helix}, {r.spread_A:.2f} A)" for r in out_cand.itertuples()) + ".")
    L.append(f"- Lowest-spread 123 (size of v3, max spread {n['spread_max_in_lowest_123_A']:.2f} A): "
             + ", ".join(f"{a}-{b}" for a, b in core["lowest_spread_123"]["core_ranges"]) + ".")
    L.append("")
    L.append("## Step 4 - comparison with v3 (core_comparison.csv)")
    L.append("")
    L.append(f"Jaccard v4 vs v3 = {n['jaccard_v4_v3']:.3f}; lowest-spread 123 vs v3 = {n['jaccard_lowest123_v3']:.3f}. "
             "Helices grouped by the helix-table ranges:")
    L.append("")
    L.append("| helix | v3 | v4 | both | v3 only | v4 only | lowest-123 | lowest-123 not in v3 | v3 not in lowest-123 |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for r in comp.itertuples():
        L.append(f"| {r.helix} | {r.v3_n} | {r.v4_n} | {r.both_v3_v4} | {rng(r.v3_only_positions)} | {rng(r.v4_only_positions)} | "
                 f"{r.lowest123_n} | {rng(r.lowest123_only_positions)} | {rng(r.v3_not_lowest123_positions)} |")
    L.append("")
    v3_out = sp[sp.in_core_v3 & ~sp.in_core_v4]
    L.append(f"- v3 positions not in v4: {len(v3_out)}; their observation fraction {v3_out.obs_frac.min():.1%}-"
             f"{v3_out.obs_frac.max():.1%} and spread {v3_out.spread_A.min():.2f}-{v3_out.spread_A.max():.2f} A "
             f"(v3 positions in v4: spread max {sp.spread_A[sp.in_core_v3 & sp.in_core_v4].max():.2f} A). "
             f"{int((v3_out.obs_frac < 0.98).sum())} of them are observed in < 98 % of chains, "
             f"{int((~v3_out.in_consensus_helix).sum())} lie outside the SV02 consensus helices.")
    L.append(f"- Helices in v4 but not in v3: " + ", ".join(f"{r.helix} {r.v4_ranges}" for r in new_helices.itertuples())
             + f" (v4 has {len(sp[sp.in_core_v4 & ~sp.in_core_v3])} positions that v3 lacks; "
             f"{int(new_helices.v4_n.sum())} of them in these helices).")
    L.append("")
    L.append("## Step 5 - recomputed with v4 (rmsd_v4.npy, ari.csv, cluster_sizes.csv)")
    L.append("")
    pc = n["s03_path_check_v3"]
    L.append(f"- Code path of results/core_sensitivity.csv (K09): S03 fit_all -> S04 compute_distances -> S08 structural_groups. "
             f"Check: the same fit_all call with core v3 reproduces results/alignment/aligned_coordinates.npz "
             f"(max |d| {pc['max_abs_diff_A']:.1e} A, same gaps {pc['same_gaps']}).")
    L.append(f"- v4 fit: status {n['v4_fit_status']}, shared core C-alpha min {n['v4_n_common_core_min']}, core RMSD to 7NFB|B|1 "
             f"median {n['v4_core_rmsd_median_A']:.3f} A, max {n['v4_core_rmsd_max_A']:.3f} A. Matrix 788 x 788, "
             f"{n['n_distance_exclusions']} distance exclusions, no NaN.")
    L.append(f"- Pairwise RMSD v4 vs primary (v3, results/rmsd_entry.npy), {n['n_pairs']} pairs: Spearman "
             f"{n['spearman_pairwise_v4_vs_v3']:.4f}; |dRMSD| median {n['median_abs_diff_A']:.4f} A, 95th percentile "
             f"{n['p95_abs_diff_A']:.3f} A, max {n['max_abs_diff_A']:.3f} A; pairs with |dRMSD| > 0.5 A: "
             f"{n['n_pairs_abs_diff_gt_0_5A']}.")
    L.append("")
    L.append("| cut (A) | clusters S08 (v3) | clusters v4 | ARI v4 vs S08 | chains changing cluster | sizes S08 | sizes v4 |")
    L.append("|---|---|---|---|---|---|---|")
    for r in ari.itertuples():
        L.append(f"| {r.cut_A:.1f} | {r.n_clusters_v3_s08} | {r.n_clusters_v4} | {r.ari_v4_vs_s08:.3f} | "
                 f"{r.chains_changed_best_match} | {r.sizes_v3_s08.replace(';', ', ')} | {r.sizes_v4.replace(';', ', ')} |")
    L.append("")
    L.append("'Chains changing cluster' = chains whose v4 cluster's best-matching S08 cluster (most shared chains) is not "
             "their own S08 cluster (the SV00 gate rule).")
    for c in ["1.5", "2.0", "2.5"]:
        m = moved[c]
        L.append(f"- {c} A: " + "; ".join(f"{r.chain_key} from S08 {c} A cluster {name(c, int(getattr(r, '_2')))} to the v4 cluster "
                                         f"that best matches S08 {c} A cluster {name(c, int(r.to))}" for r in m.itertuples()) + ".")
    parts = []
    for k, row in big.iterrows():
        src = ", ".join(f"{int(v)} from S08 1.0 A cl{int(s)}" for s, v in row.items() if v > 0)
        parts.append(f"v4 cluster {int(k)} ({int(row.sum())} chains: {src}; S08 1.5 A cluster(s) {within[int(k)]})")
    L.append(f"- 1.0 A ({int(ari.chains_changed_best_match.iloc[0])} chains change): " + "; ".join(parts) + ".")
    L.append("")
    L.append("## Step 5 - Y10-style displacement profile on v4 (displacement_profile_v4.csv)")
    L.append("")
    L.append(f"- Medoids (complete H12, v4 matrix): {pv4['medoid_cl7']} (cl7), {pv4['medoid_cl4']} (cl4); Y10 (v3): "
             f"{pv3['medoid_cl7']}, {pv3['medoid_cl4']}.")
    L.append(f"- Share of the squared cl7-cl4 mean difference in 528-550: v4 {pv4['share_528_550']:.1%}, v3 "
             f"{pv3['share_528_550']:.1%}; in H12 538-547: v4 {pv4['share_h12_538_547']:.1%}, v3 {pv3['share_h12_538_547']:.1%}.")
    L.append(f"- Largest cl7-cl4 distance before 528: v4 {pv4['max_d_mean_outside_528_A']:.2f} A at "
             f"{pv4['position_max_d_mean_outside_528']}, v3 {pv3['max_d_mean_outside_528_A']:.2f} A at "
             f"{pv3['position_max_d_mean_outside_528']}. Top 5 positions v4 {pv4['top_positions']} "
             f"({pv4['top_d_mean_A'][0]:.1f} A max), v3 {pv3['top_positions']} ({pv3['top_d_mean_A'][0]:.1f} A max).")
    L.append(f"- |d_mean v4 - v3| max {n['profile_max_abs_delta_d_mean_A']:.2f} A (position {n['profile_position_max_abs_delta']}); "
             f"before 528 max {n['profile_max_abs_delta_outside_528_A']:.2f} A.")
    L.append("")
    L.append("## Figure")
    L.append("")
    L.append("figures/sv/SV05/SV05_core_v4.png (+ .svg; panels in figures/panels/sv/SV05/): (a) spread per position after the "
             "iterative fit with candidates and the 0.5 A line, (b) core v3 / v4 / lowest-123 positions, (c) Y10-style "
             "cl7-cl4 distance per position v3 vs v4, (d) clusters per cut and ARI, (e) cluster sizes at 1.5 A.")
    L.append("")
    L.append("## Files")
    L.append("")
    files = sorted(OUT.glob("*")) + sorted((ROOT / "data/sv/SV05").glob("*")) + sorted((ROOT / "figures/sv/SV05").glob("*")) \
        + sorted((ROOT / "figures/panels/sv/SV05").glob("*")) + sorted((ROOT / "scripts/sv").glob("SV05_*.py"))
    for f in files:
        if f.name == "SUMMARY.md":
            continue
        L.append(f"- {f.relative_to(ROOT)} ({size(f)})")
    L.append("- results/sv/SV05/SUMMARY.md (this file)")
    L.append("- docs/EXPLORATION.md: section 'Suggested, not done (SV05, 2026-11-08)' appended (rule 1.1); existing text unchanged")
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
