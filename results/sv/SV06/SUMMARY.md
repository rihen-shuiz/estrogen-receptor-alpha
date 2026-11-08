Not done / substituted: (1) panel B medoids: the S12 medoids of results/medoids.csv were applied first; their cl4 medoid 2QE4|B|1 has H12 not_modelled and lies 0.36 A from the cl7 medoid (median cl7-cl4 chain RMSD 3.49 A; the S04 RMSD uses the C-alpha positions two chains share), a degenerate reference, so panel B uses the Y10 complete-H12 medoids (5TM1|A|1, 7UJ7|C|1; the pair the Y15 medoid detector used). This substitution was made after the first run; the S12 distances stay in outliers_on_map.csv (*_s12 columns). (2) Y16 reasons come from results/ext/outlier_verdicts.csv, which is not an `out` of the Y16 card but is listed as evidence in the Y16 DONE line of tracker/esr1.log; it equals the verdicts in numbers_Y16.json (checked). (3) The run file names no groups, units or verdict rule for the hypergeometric test. Declared before running: groups cl7, cl4 and minor (5 clusters pooled), plus cl7 vs cl4 within the main clusters; entry level = (entry, cluster group) pairs; BH over the 8 tests; 'pattern' = q < 0.05 at chain and entry level. (4) No threshold for 'near the diagonal' is given, so none is applied; the own/other ratio is reported. (5) Labels are PDB|author chain; chain_key holds the label_asym_id (6PFM|A|1 = 6PFM|D). (6) figures/ext is read-only, so the Y07 heatmap was redrawn with Y07's own bar code (scripts/ext_dendrogram.py bars()) and one extra outlier bar. scripts/sv/common/sv_fig.py save_figure() now also crops a panel from a list of axes (backward compatible).

# SV06 - Consensus outliers on the clustering map

**Verdict: pattern between the main clusters** (outliers are fewer in cl7 and more in cl4 than expected). The minor clusters are enriched, as they must be by construction of the Y15 detectors.

## Inputs (rule 1.10)

- Y15 card out results/ext/outlier_flags.csv (consensus_outlier, 22 chains); Y16 reasons from results/ext/outlier_verdicts.csv (Y16 DONE evidence).
- Y03 card out results/ext/embedding_rmsd.csv, mds1/mds2 (shared map = metric MDS, Y03 DONE note).
- Y07 card out results/ext/leaf_order.csv; S08 results/entry_clusters.csv (cut 1.5 = SV00 s08_cluster, checked); S12 results/medoids.csv; Y10 results/ext/cluster_displacement_profile.csv (complete-H12 medoids); S04/S05 results/rmsd_entry.npy + rmsd_entry_index.csv (loaded, not recomputed); SV00 results/sv/SV00/chain_table.csv.

## Key numbers

- 22 consensus outliers in 17 entries: cl7 2 of 445 chains, cl4 12 of 329, minor 8 of 14.
- Y16 reason by group (quality / conformation): cl7 0 / 2, cl4 3 / 9, minor 3 / 5.

Hypergeometric tests (two-sided = Fisher exact; hypergeom_tests.csv):

| level | test | units | in group | outliers | observed | expected | obs/exp | p | q (BH) |
|---|---|---|---|---|---|---|---|---|---|
| chain | cl7 vs all | 788 | 445 | 22 | 2 | 12.42 | 0.16 | 6.3e-06 | 1.7e-05 |
| chain | cl4 vs all | 788 | 329 | 22 | 12 | 9.19 | 1.31 | 0.27 | 0.27 |
| chain | minor vs all | 788 | 14 | 22 | 8 | 0.39 | 20.47 | 2.5e-10 | 2e-09 |
| chain | cl7 vs cl4 (main only) | 774 | 445 | 14 | 2 | 8.05 | 0.25 | 0.0015 | 0.0023 |
| entry | cl7 vs all | 356 | 221 | 17 | 2 | 10.55 | 0.19 | 1.8e-05 | 3.5e-05 |
| entry | cl4 vs all | 356 | 128 | 17 | 9 | 6.11 | 1.47 | 0.19 | 0.22 |
| entry | minor vs all | 356 | 7 | 17 | 6 | 0.33 | 17.95 | 3.1e-08 | 1.2e-07 |
| entry | cl7 vs cl4 (main only) | 349 | 221 | 11 | 2 | 6.97 | 0.29 | 0.0026 | 0.0034 |

Calls (q < 0.05 at both levels = pattern): cl7 vs cl4 (main only): pattern; cl7 vs all: pattern; cl4 vs all: no pattern; minor vs all: pattern.

Panel B (complete-H12 medoids, own main medoid vs the other):

- medoids 5TM1|A|1 (cl7) and 7UJ7|C|1 (cl4), 3.68 A apart.
- 0 of 14 main-cluster outliers are closer to the other medoid; own/other ratio 0.21-0.63 (median 0.25); minor-cluster outliers 0.33-0.79.
- With the S12 medoids (degenerate, not used) 2 of 14 would be 'closer to the other medoid'.

Outliers nearest the diagonal (highest own/other ratio):

| label | cluster | Y16 reason | own medoid RMSD | other medoid RMSD | ratio |
|---|---|---|---|---|---|
| 1A52\|A | minor | quality: resolution_A 2.8 > p95 2.7 | 5.53 | 7.02 | 0.79 |
| 1A52\|B | minor | quality: resolution_A 2.8 > p95 2.7 | 5.45 | 7.14 | 0.76 |
| 6PFM\|D | minor | quality: resolution_A 2.84 > p95 2.7; mean_b_chain 70.7 > p95 66.4 | 5.82 | 7.70 | 0.76 |
| 1R5K\|C | minor | conformation: none | 1.38 | 1.95 | 0.71 |
| 9B25\|C | minor | conformation: none | 2.37 | 3.44 | 0.69 |
| 2BJ4\|B | cl4 | quality: H12 not modelled | 0.94 | 1.49 | 0.63 |

## verdict.md

> cl7 vs cl4: pattern - of the 14 consensus outliers inside the main clusters, 2 are in cl7 and 12 in cl4 (chain 2 vs 8.0 expected, q 0.0023; entry 2 vs 7.0, q 0.0034; hypergeometric, BH over 8 tests).
> Minor clusters hold 8 of 22 outliers from 14 chains (pattern: chain 8 vs 0.4 expected, q 2e-09; entry 6 vs 0.3, q 1.2e-07); expected by construction, since the Y15 detectors measure distance from the main clusters.
> Panel B (complete-H12 medoids 5TM1|A|1 / 7UJ7|C|1, 3.68 A apart): 0 of 14 main-cluster outliers are closer to the other medoid than to their own (own/other RMSD ratio 0.21-0.63, median 0.25; highest 2BJ4|B 0.63, Y16 reason: H12 not modelled); the S12 cl4 medoid 2QE4|B|1 has no modelled H12 and is not used.

## Files

- results/sv/SV06/outliers_on_map.csv (7.6 kB)
- results/sv/SV06/hypergeom_tests.csv (3.0 kB)
- results/sv/SV06/numbers.json (3.1 kB)
- results/sv/SV06/verdict.md (0.7 kB)
- results/sv/SV06/env.txt (0.3 kB)
- figures/sv/SV06/SV06_outliers_map.png (390.1 kB)
- figures/sv/SV06/SV06_outliers_map.svg (181.2 kB)
- figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.png (1.11 MB)
- figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.svg (1.30 MB)
- figures/panels/sv/SV06/SV06_outliers_map__A_mds_map.png (195.9 kB)
- figures/panels/sv/SV06/SV06_outliers_map__B_own_vs_other_medoid.png (160.2 kB)
- figures/panels/sv/SV06/SV06_rmsd_heatmap_outlier_bar__heatmap_with_bars.png (1.04 MB)
- scripts/sv/SV06_outliers.py (18.3 kB)
- scripts/sv/SV06_figures.py (19.6 kB)
- scripts/sv/SV06_summary.py (7.6 kB)
- scripts/sv/common/sv_fig.py (edited: multi-axes panel crop) (2.4 kB)
- results/sv/SV06/SUMMARY.md (this file)

Run order: .venv/bin/python -B -I scripts/sv/SV06_outliers.py; scripts/sv/SV06_figures.py; scripts/sv/SV06_summary.py
