"""SV06: write results/sv/SV06/SUMMARY.md (rule 1.7) from the saved SV06 files only.
Run after SV06_outliers.py and SV06_figures.py: .venv/bin/python -B -I scripts/sv/SV06_summary.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/sv/SV06"


def size(p):
    s = p.stat().st_size
    return f"{s / 1e3:.1f} kB" if s < 1e6 else f"{s / 1e6:.2f} MB"


def main():
    n = json.loads((OUT / "numbers.json").read_text())
    mc = n["medoid_check"]
    o = pd.read_csv(OUT / "outliers_on_map.csv")
    ht = pd.read_csv(OUT / "hypergeom_tests.csv")
    verdict = (OUT / "verdict.md").read_text().strip().splitlines()
    relabel = o[o.chain_key.str.split("|").str[1] != o.auth_asym_id.astype(str)]
    rl = ", ".join(f"{r.chain_key} = {r.label}" for r in relabel.itertuples()) or "none"

    L = [("Not done / substituted: (1) panel B medoids: the S12 medoids of results/medoids.csv were applied first; "
          f"their cl4 medoid {mc['s12_cl4']} has H12 {mc['s12_cl4_h12_label']} and lies {mc['rmsd_s12_cl7_cl4']:.2f} A "
          f"from the cl7 medoid (median cl7-cl4 chain RMSD {mc['median_rmsd_cl7_cl4_chains']:.2f} A; the S04 RMSD uses "
          "the C-alpha positions two chains share), a degenerate reference, so panel B uses the Y10 complete-H12 "
          f"medoids ({mc['y10_cl7']}, {mc['y10_cl4']}; the pair the Y15 medoid detector used). This substitution was "
          "made after the first run; the S12 distances stay in outliers_on_map.csv (*_s12 columns). "
          "(2) Y16 reasons come from results/ext/outlier_verdicts.csv, which is not an `out` of the Y16 card but is "
          "listed as evidence in the Y16 DONE line of tracker/esr1.log; it equals the verdicts in numbers_Y16.json "
          "(checked). (3) The run file names no groups, units or verdict rule for the hypergeometric test. Declared "
          "before running: groups cl7, cl4 and minor (5 clusters pooled), plus cl7 vs cl4 within the main clusters; "
          "entry level = (entry, cluster group) pairs; BH over the 8 tests; 'pattern' = q < 0.05 at chain and entry "
          "level. (4) No threshold for 'near the diagonal' is given, so none is applied; the own/other ratio is "
          f"reported. (5) Labels are PDB|author chain; chain_key holds the label_asym_id ({rl}). (6) figures/ext is "
          "read-only, so the Y07 heatmap was redrawn with Y07's own bar code (scripts/ext_dendrogram.py bars()) and "
          "one extra outlier bar. scripts/sv/common/sv_fig.py save_figure() now also crops a panel from a list of "
          "axes (backward compatible)."),
         "", "# SV06 - Consensus outliers on the clustering map", "",
         f"**Verdict: {n['calls']['cl7 vs cl4 (main only)']['call']} between the main clusters** "
         f"(outliers are fewer in cl7 and more in cl4 than expected). The minor clusters are enriched, as they must "
         f"be by construction of the Y15 detectors.", ""]
    L += ["## Inputs (rule 1.10)", "",
          "- Y15 card out results/ext/outlier_flags.csv (consensus_outlier, 22 chains); Y16 reasons from "
          "results/ext/outlier_verdicts.csv (Y16 DONE evidence).",
          "- Y03 card out results/ext/embedding_rmsd.csv, mds1/mds2 (shared map = metric MDS, Y03 DONE note).",
          "- Y07 card out results/ext/leaf_order.csv; S08 results/entry_clusters.csv (cut 1.5 = SV00 s08_cluster, "
          "checked); S12 results/medoids.csv; Y10 results/ext/cluster_displacement_profile.csv (complete-H12 "
          "medoids); S04/S05 results/rmsd_entry.npy + rmsd_entry_index.csv (loaded, not recomputed); SV00 "
          "results/sv/SV00/chain_table.csv.", ""]
    rb = n["reason_by_group"]
    L += ["## Key numbers", "",
          f"- {n['n_outliers']} consensus outliers in {n['n_outlier_entries']} entries: cl7 {n['outliers_by_group']['cl7']} "
          f"of {n['chains_by_group']['cl7']} chains, cl4 {n['outliers_by_group']['cl4']} of "
          f"{n['chains_by_group']['cl4']}, minor {n['outliers_by_group']['minor']} of {n['chains_by_group']['minor']}.",
          f"- Y16 reason by group (quality / conformation): cl7 {rb['cl7']['quality']} / {rb['cl7']['conformation']}, "
          f"cl4 {rb['cl4']['quality']} / {rb['cl4']['conformation']}, minor {rb['minor']['quality']} / "
          f"{rb['minor']['conformation']}.", "",
          "Hypergeometric tests (two-sided = Fisher exact; hypergeom_tests.csv):", "",
          "| level | test | units | in group | outliers | observed | expected | obs/exp | p | q (BH) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in ht[ht.tested].itertuples():
        L.append(f"| {r.level} | {r.test} | {r.n_units} | {r.n_group} | {r.n_outliers} | {r.observed} | "
                 f"{r.expected:.2f} | {r.obs_over_exp:.2f} | {r.p_two_sided:.2g} | {r.q_bh:.2g} |")
    L += ["", "Calls (q < 0.05 at both levels = pattern): " + "; ".join(
        f"{k}: {v['call']}" for k, v in n["calls"].items()) + ".", ""]
    L += ["Panel B (complete-H12 medoids, own main medoid vs the other):", "",
          f"- medoids {n['medoid_cl7']} (cl7) and {n['medoid_cl4']} (cl4), {n['rmsd_medoid_cl7_cl4']:.2f} A apart.",
          f"- {n['n_main_closer_to_other']} of {n['n_main']} main-cluster outliers are closer to the other medoid; "
          f"own/other ratio {n['ratio_main_min']:.2f}-{n['ratio_main_max']:.2f} (median {n['ratio_main_median']:.2f}); "
          f"minor-cluster outliers {n['ratio_minor_min']:.2f}-{n['ratio_minor_max']:.2f}.",
          f"- With the S12 medoids (degenerate, not used) {n['n_main_closer_to_other_s12']} of {n['n_main']} would "
          "be 'closer to the other medoid'.", ""]
    top = o.sort_values("own_other_ratio", ascending=False).head(6)
    L += ["Outliers nearest the diagonal (highest own/other ratio):", "",
          "| label | cluster | Y16 reason | own medoid RMSD | other medoid RMSD | ratio |", "|---|---|---|---|---|---|"]
    for r in top.itertuples():
        lab = r.label.replace("|", "\\|")          # a bare | would split the markdown cell
        L.append(f"| {lab} | {r.cluster_group} | {r.reason}: {r.reasons} | {r.rmsd_own_medoid:.2f} | "
                 f"{r.rmsd_other_medoid:.2f} | {r.own_other_ratio:.2f} |")
    L += ["", "## verdict.md", ""] + [f"> {v}" for v in verdict] + [""]
    files = ["results/sv/SV06/outliers_on_map.csv", "results/sv/SV06/hypergeom_tests.csv",
             "results/sv/SV06/numbers.json", "results/sv/SV06/verdict.md", "results/sv/SV06/env.txt",
             "figures/sv/SV06/SV06_outliers_map.png", "figures/sv/SV06/SV06_outliers_map.svg",
             "figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.png", "figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.svg",
             "figures/panels/sv/SV06/SV06_outliers_map__A_mds_map.png",
             "figures/panels/sv/SV06/SV06_outliers_map__B_own_vs_other_medoid.png",
             "figures/panels/sv/SV06/SV06_rmsd_heatmap_outlier_bar__heatmap_with_bars.png",
             "scripts/sv/SV06_outliers.py", "scripts/sv/SV06_figures.py", "scripts/sv/SV06_summary.py",
             "scripts/sv/common/sv_fig.py (edited: multi-axes panel crop)"]
    L += ["## Files", ""]
    for f in files:
        path = ROOT / f.split(" (")[0]
        L.append(f"- {f} ({size(path)})")
    L += ["- results/sv/SV06/SUMMARY.md (this file)", "",
          "Run order: .venv/bin/python -B -I scripts/sv/SV06_outliers.py; scripts/sv/SV06_figures.py; "
          "scripts/sv/SV06_summary.py"]
    (OUT / "SUMMARY.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
