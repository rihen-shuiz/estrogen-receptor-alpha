Not done / substituted: (1) 'Always include the five main SV figures (SV09)': the SV09 table lists 10 files (with the SV13 per-class H12 panels, the SV13 forest plot and the SV14 main figure of the 2026-11-04 addition); all 10 are in visuals/, multi-graph ones as their VIS1 panels, except 3 empty SV13 class panels (antagonist, other, partial_agonist: 'no chain carries this class', text only), which are not copied. (2) Two claims use a keep-all (788-chain) image from figures/ that RESULTS_SUMMARY does not cite, because it prints the stated number: 03_dendrogram_1p5_cut.png from figures/dendrogram.png; 40_cluster_discovery_curve.png from figures/diversity_curve.png (the 7 clusters at 1.5 A, with the cited Y07 heatmap kept for the 3.49 A claim; 90% of the 1.5 A clusters at N=488, in place of the cited Y26_rarefaction panel 01, which shows the curve without the number). (3) The Y30 B-putty pair (2 panels, cited in section 4 'PyMOL putty') is excluded: the same claim (H12 / core-v3 B ratio 1.53 vs 1.33) is shown for all chains by Y21_helix_b; so only the Y30 medoids legend is needed (suffix b). (4) The claim check reads the documents with line breaks as single spaces and \| as | (RESULTS_SUMMARY wraps lines; GFM tables need the escape) and checks each quote inside its cited section or SV paragraph; INDEX.md itself is not counted as a visual. (5) Supporting files not named in the run file: results/sv/VIS2/selection.csv, results/sv/VIS2/check.csv, scripts/sv/VIS2_build.py. Nothing installed.

# VIS2 - Collect the valuable visuals into one folder

**Verdict**

- Pool (VIS1 manifest): 314 items - not_cut 9, png_gutter 219, single 58, sv_script 28.
- Copied into visuals/: 77 files = 76 graphs + 1 legend (suffix b); by pool method not_cut 8, png_gutter 44, single 14, sv_script 11. Excluded: 237.
  - duplicate (another graph kept for the claim): 119
  - tied to no stated claim (QC, diagnostic, context): 65
  - superseded or September (351-entry) version: 49
  - empty panel (no data): 3
  - legend-only file: 1
- Story order (RESULTS_SUMMARY.md): The data 01-02 (2); section 1 03-11 (9); section 2 12-24 (13); section 3 25-29 (5); section 4 30-32 (3); section 5 33-37 (5); section 6 38-43 (6); section 7 44-47 (4); Supervisor requests (SV) 48-76 (29).
- INDEX claims: 91 verbatim quotes in 77 rows, from docs/RESULTS_SUMMARY.md.
- Check (scripts/sv/VIS2_check.py, results/sv/VIS2/check.csv): claim 91 pass / 0 fail; listed 77 pass / 0 fail; numbering 77 pass / 0 fail; opens 77 pass / 0 fail; sha1 1 pass / 0 fail; source 77 pass / 0 fail; sv09_main 11 pass / 0 fail. All pass.
- Distinct SHA-1 among the 77 PNG files: 77 distinct SHA-1; every copy equals its source SHA-1 (check 'source').

## Main SV figures (SV09 table) in visuals/

  - figures/sv/SV01b/SV01b_forest_completeness.png: original -> 49_sv01b_forest_completeness.png
  - figures/sv/SV01b/SV01b_forest_mobility.png: original -> 50_sv01b_forest_mobility.png
  - figures/sv/SV02/SV02_helix_track.png: original -> 51_sv02_helix_track.png
  - figures/sv/SV02/rama/SV02_rama_537.png: original -> 52_sv02_rama_537.png
  - figures/sv/SV04/h12_map_partial.png: 2 of 2 VIS1 panels -> 55_sv04_h12_map_partial_a.png, 56_sv04_h12_map_partial_b.png
  - figures/sv/SV06/SV06_outliers_map.png: original -> 61_sv06_outliers_map.png
  - figures/sv/SV07/SV07_spike_profile.png: 2 of 2 VIS1 panels -> 62_sv07_spike_profile.png, 63_sv07_observation_share.png
  - figures/sv/SV13/SV13_forest.png: 3 of 3 VIS1 panels -> 65_sv13_forest_c1_complete.png, 66_sv13_forest_c1_mobility_c2a.png, 67_sv13_forest_c2b.png
  - figures/sv/SV13/SV13_h12_class_panels.png: 4 of 7 VIS1 panels -> 68_sv13_h12_panel_agonist.png, 69_sv13_h12_panel_serm.png, 70_sv13_h12_panel_serd.png, 71_sv13_h12_class_overlay.png
  - figures/sv/SV14/SV14_pocket_vs_ligand.png: original -> 73_sv14_pocket_vs_ligand.png

## Selection rule as applied (scripts/sv/VIS2_build.py)

- One graph per claim stated in docs/RESULTS_SUMMARY.md or an SV verdict file: for a claim cited with a multi-graph figure, the panel that shows the stated number; where two graphs carry the same claim, the one that states it (reason per pool item in selection.csv, e.g. 'same claim as NN_...').
- Keep-all (788-chain) graphs over the September (351-entry) core-task and provisional images (figures/0*-11*, figures/provisional_*, figures/qc_outliers, results/figures/, results/archive_v2_strict/); graphs at the declared 1.5 A cut (multi-cut graphs only where the claim is the 1.5 A value).
- not_cut originals are copied as they are (one file, several graphs): 26_outlier_quality.png, 30_helix_b_ranking.png, 34_ligand_physchem_by_cluster.png, 39_set_cover.png, 49_sv01b_forest_completeness.png, 50_sv01b_forest_mobility.png, 61_sv06_outliers_map.png, 74_sv14_filling.png.

## Files

- visuals/INDEX.md (32.4 KB)
- results/sv/VIS2/selection.csv (74.0 KB)
- results/sv/VIS2/check.csv (30.9 KB)
- results/sv/VIS2/SUMMARY.md (5.1 KB)
- results/sv/VIS2/env.txt (0.3 KB)
- scripts/sv/VIS2_build.py (43.7 KB)
- scripts/sv/VIS2_check.py (7.9 KB)
- scripts/sv/VIS2_summary.py (9.9 KB)
- visuals/*.png: 77 files (12.4 MB)
- docs/EXPLORATION.md: appended 'Suggested, not done (VIS2, 2026-11-08)' (rule 1.1)

Status proposed: DONE
