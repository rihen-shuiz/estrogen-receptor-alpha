Not done / substituted: (1) interpretations fixed before the run (the step names no formula): 'mean shift' = mean |x_new - x_old| over all observed (chain, candidate) C-alpha between two iterations; each iteration Kabsch-fits every chain on its observed candidates onto the nan-mean candidate structure and then re-anchors the whole ensemble on the starting core-v3-frame mean (one common rigid motion, which leaves every spread unchanged); 'spread' = RMS deviation of the observing chains from the mean position (the Y10 sd formula); ties in the 123 lowest-spread set broken by position. (2) The Y10-style profile re-implements scripts/ext_displacement_profile.py (same formulas) instead of running it, because that script writes into results/ext and figures/ext; on the v3 inputs it reproduces results/ext/cluster_displacement_profile.csv to 3.6e-15 with the same medoids. Clusters there are the S08 1.5 A cl7/cl4 labels; the complete-H12 medoid rule is applied on the v4 matrix. (3) The 123 lowest-spread set is reported as a set (core_v4.json, spread.csv, core_comparison.csv, Jaccard); no matrix or clustering was computed for it (step 5 names core v4 only). (4) 'Which helices differ' groups positions by the helix-table ranges in results/sv/SV02/helix_boundaries.csv (table_start-table_end), because v3 includes helix ends outside the SV02 consensus helices. (5) Supporting files not named in the run file: rmsd_v4_index.csv, clusters_v4.csv, cluster_sizes.csv, core_comparison.csv, core_sensitivity_v4.csv, displacement_profile_v4.csv, iterations.csv, numbers_core_v4.json, data/sv/SV05/. The superposition core of SV03/SV04/SV07 stays v3; nothing existing was overwritten.

# SV05 - Stable core v4 from complete helices, compared with core v3

**Verdict ([DEFAULT] step 3):** core v4 = **128 positions** with spread <= 0.5 A (>= 80, not BLOCKED), from 132 candidates. Jaccard with v3 (123 positions) = **0.64** (lowest-spread 123: 0.61). v4 adds H1 and H9 and drops 25 v3 positions at helix ends. Recomputed with v4 through the K09 code path: pairwise RMSD Spearman vs the primary matrix **0.9975**, median |dRMSD| 0.015 A (max 0.52 A); ARI against S08 at 1.0/1.5/2.0/2.5 A = 0.817/0.990/0.992/0.992; clusters 15/6/4/3 (v4) vs 16/7/4/3 (S08). At 1.5 A 3 chains change cluster; the two main clusters keep 447/329 chains (S08 445/329). The Y10-style profile on v4 puts 98.5% of the squared cl7-cl4 difference in 528-550 (v3 98.6%). Adopting v4 or keeping v3 is a team decision (run file FINAL); this step only compares.

## Step 1 - candidates (spread.csv)

- SV02 global consensus helices: 140 positions (results/sv/SV02/consensus_global.csv); in 528-550: 7 (all of H12 538-544); observed in < 98 % of the 788 chains: 1 (526, 96.8%).
- Candidates: **132** positions in 7 helices (H1, H3, H4/5, H7, H8, H9, H10/11).

## Step 2 - iterative superposition (iterations.csv)

- Start: core v3 frame of the SV00 gate (data/sv/gate/aligned_coordinates.npz, reference 7NFB|B|1). Converged after **2 iterations**: mean shift 0.11, 3e-06 A (stop < 0.01 A); mean per-chain RMS to the mean structure on the candidates 0.242 A.
- Candidate spread: min 0.10, median 0.19, max 0.78 A (in the v3 frame before iterating: median 0.26, max 0.74 A).

## Step 3 - core v4 (core_v4.json)

- Core v4 = 128 positions: 312-321, 343-361, 372-394, 423-437, 442-455, 473-492, 497-523.
- Candidates above 0.5 A (left out): 342 (H3, 0.54 A), 422 (H7, 0.51 A), 524 (H10/11, 0.69 A), 525 (H10/11, 0.78 A).
- Lowest-spread 123 (size of v3, max spread 0.38 A): 312-321, 346-361, 372-394, 424-437, 442-455, 473-492, 497-522.

## Step 4 - comparison with v3 (core_comparison.csv)

Jaccard v4 vs v3 = 0.641; lowest-spread 123 vs v3 = 0.608. Helices grouped by the helix-table ranges:

| helix | v3 | v4 | both | v3 only | v4 only | lowest-123 | lowest-123 not in v3 | v3 not in lowest-123 |
|---|---|---|---|---|---|---|---|---|
| H1 | 0 | 10 | 0 | - | 312-321 | 10 | 312-321 | - |
| H3 | 25 | 19 | 19 | 338-342 362 | - | 16 | - | 338-345 362 |
| H4/5 | 26 | 23 | 23 | 371 395-396 | - | 23 | - | 371 395-396 |
| H7 | 20 | 15 | 15 | 420-422 438-439 | - | 14 | - | 420-423 438-439 |
| H8 | 16 | 14 | 14 | 441 456 | - | 14 | - | 441 456 |
| H9 | 0 | 20 | 0 | - | 473-492 | 20 | 473-492 | - |
| H10/11 | 36 | 27 | 27 | 496 524-531 | - | 26 | - | 496 523-531 |

- v3 positions not in v4: 25; their observation fraction 64.5%-100.0% and spread 0.17-4.34 A (v3 positions in v4: spread max 0.50 A). 10 of them are observed in < 98 % of chains, 20 lie outside the SV02 consensus helices.
- Helices in v4 but not in v3: H1 312-321, H9 473-492 (v4 has 30 positions that v3 lacks; 30 of them in these helices).

## Step 5 - recomputed with v4 (rmsd_v4.npy, ari.csv, cluster_sizes.csv)

- Code path of results/core_sensitivity.csv (K09): S03 fit_all -> S04 compute_distances -> S08 structural_groups. Check: the same fit_all call with core v3 reproduces results/alignment/aligned_coordinates.npz (max |d| 3.6e-14 A, same gaps True).
- v4 fit: status {'ok': 788}, shared core C-alpha min 126, core RMSD to 7NFB|B|1 median 0.398 A, max 0.564 A. Matrix 788 x 788, 0 distance exclusions, no NaN.
- Pairwise RMSD v4 vs primary (v3, results/rmsd_entry.npy), 310078 pairs: Spearman 0.9975; |dRMSD| median 0.0149 A, 95th percentile 0.105 A, max 0.525 A; pairs with |dRMSD| > 0.5 A: 5.

| cut (A) | clusters S08 (v3) | clusters v4 | ARI v4 vs S08 | chains changing cluster | sizes S08 | sizes v4 |
|---|---|---|---|---|---|---|
| 1.0 | 16 | 15 | 0.817 | 66 | 442, 262, 48, 11, 4, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1 | 442, 211, 114, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1 |
| 1.5 | 7 | 6 | 0.990 | 3 | 445, 329, 5, 3, 3, 2, 1 | 447, 329, 4, 3, 3, 2 |
| 2.0 | 4 | 4 | 0.992 | 2 | 445, 336, 5, 2 | 447, 335, 4, 2 |
| 2.5 | 3 | 3 | 0.992 | 2 | 447, 336, 5 | 449, 335, 4 |

'Chains changing cluster' = chains whose v4 cluster's best-matching S08 cluster (most shared chains) is not their own S08 cluster (the SV00 gate rule).
- 1.5 A: 2JF9|A|1 from S08 1.5 A cluster 1 (5 chains, 1.5 A label of most: cl1) to the v4 cluster that best matches S08 1.5 A cluster 7 (445 chains, 1.5 A label of most: cl7); 2JFA|B|1 from S08 1.5 A cluster 4 (329 chains, 1.5 A label of most: cl4) to the v4 cluster that best matches S08 1.5 A cluster 7 (445 chains, 1.5 A label of most: cl7); 3DT3|B|1 from S08 1.5 A cluster 5 (1 chains, 1.5 A label of most: cl5) to the v4 cluster that best matches S08 1.5 A cluster 4 (329 chains, 1.5 A label of most: cl4).
- 2.0 A: 2JF9|A|1 from S08 2.0 A cluster 1 (5 chains, 1.5 A label of most: cl1) to the v4 cluster that best matches S08 2.0 A cluster 4 (445 chains, 1.5 A label of most: cl7); 2JFA|B|1 from S08 2.0 A cluster 2 (336 chains, 1.5 A label of most: cl4) to the v4 cluster that best matches S08 2.0 A cluster 4 (445 chains, 1.5 A label of most: cl7).
- 2.5 A: 2JF9|A|1 from S08 2.5 A cluster 1 (5 chains, 1.5 A label of most: cl1) to the v4 cluster that best matches S08 2.5 A cluster 3 (447 chains, 1.5 A label of most: cl7); 2JFA|B|1 from S08 2.5 A cluster 2 (336 chains, 1.5 A label of most: cl4) to the v4 cluster that best matches S08 2.5 A cluster 3 (447 chains, 1.5 A label of most: cl7).
- 1.0 A (66 chains change): v4 cluster 7 (114 chains: 4 from S08 1.0 A cl8, 110 from S08 1.0 A cl11; S08 1.5 A cluster(s) [4]); v4 cluster 8 (211 chains: 2 from S08 1.0 A cl6, 10 from S08 1.0 A cl7, 47 from S08 1.0 A cl10, 152 from S08 1.0 A cl11; S08 1.5 A cluster(s) [4]); v4 cluster 12 (442 chains: 1 from S08 1.0 A cl2, 441 from S08 1.0 A cl15; S08 1.5 A cluster(s) [1, 7]).

## Step 5 - Y10-style displacement profile on v4 (displacement_profile_v4.csv)

- Medoids (complete H12, v4 matrix): 5TM1|A|1 (cl7), 7UJ7|C|1 (cl4); Y10 (v3): 5TM1|A|1, 7UJ7|C|1.
- Share of the squared cl7-cl4 mean difference in 528-550: v4 98.5%, v3 98.6%; in H12 538-547: v4 34.5%, v3 34.6%.
- Largest cl7-cl4 distance before 528: v4 3.71 A at 333, v3 3.58 A at 333. Top 5 positions v4 [550, 549, 548, 547, 546] (37.0 A max), v3 [550, 549, 548, 547, 546] (36.9 A max).
- |d_mean v4 - v3| max 0.34 A (position 531); before 528 max 0.15 A.

## Figure

figures/sv/SV05/SV05_core_v4.png (+ .svg; panels in figures/panels/sv/SV05/): (a) spread per position after the iterative fit with candidates and the 0.5 A line, (b) core v3 / v4 / lowest-123 positions, (c) Y10-style cl7-cl4 distance per position v3 vs v4, (d) clusters per cut and ARI, (e) cluster sizes at 1.5 A.

## Files

- results/sv/SV05/ari.csv (0.8 kB)
- results/sv/SV05/cluster_sizes.csv (1.4 kB)
- results/sv/SV05/clusters_v4.csv (24.9 kB)
- results/sv/SV05/core_comparison.csv (0.5 kB)
- results/sv/SV05/core_sensitivity_v4.csv (0.5 kB)
- results/sv/SV05/core_v4.json (4.5 kB)
- results/sv/SV05/displacement_profile_v4.csv (69.3 kB)
- results/sv/SV05/env.txt (0.5 kB)
- results/sv/SV05/iterations.csv (0.2 kB)
- results/sv/SV05/numbers_core_v4.json (2.1 kB)
- results/sv/SV05/rmsd_v4.npy (4.97 MB)
- results/sv/SV05/rmsd_v4_index.csv (14.1 kB)
- results/sv/SV05/spread.csv (24.3 kB)
- data/sv/SV05/aligned_coordinates_v4.npz (4.20 MB)
- data/sv/SV05/alignment_transforms_v4.json (1.69 MB)
- data/sv/SV05/global_alignment_v4.csv (39.8 kB)
- data/sv/SV05/iterative_fit_coords.npz (4.20 MB)
- data/sv/SV05/run_stdout.txt (2.2 kB)
- figures/sv/SV05/SV05_core_v4.png (466.4 kB)
- figures/sv/SV05/SV05_core_v4.svg (164.2 kB)
- figures/panels/sv/SV05/SV05_core_v4__a_spread.png (127.5 kB)
- figures/panels/sv/SV05/SV05_core_v4__b_membership.png (36.2 kB)
- figures/panels/sv/SV05/SV05_core_v4__c_displacement_profile.png (144.7 kB)
- figures/panels/sv/SV05/SV05_core_v4__d_clusters_ari.png (63.5 kB)
- figures/panels/sv/SV05/SV05_core_v4__e_cluster_sizes_1_5.png (62.6 kB)
- scripts/sv/SV05_core_v4.py (24.4 kB)
- scripts/sv/SV05_figures.py (8.6 kB)
- scripts/sv/SV05_summary.py (15.1 kB)
- results/sv/SV05/SUMMARY.md (this file)
- docs/EXPLORATION.md: section 'Suggested, not done (SV05, 2026-11-08)' appended (rule 1.1); existing text unchanged
