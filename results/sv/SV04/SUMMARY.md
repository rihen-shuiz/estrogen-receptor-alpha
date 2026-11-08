Not done / substituted: (1) 'helical' comes from SV02's pydssp fallback (mkdssp not installed; 'H' = any helix). (2) Step 1 'observed C-alpha' read as SV03's 'observed helical C-alpha' (step 1: 'fit the axis as in SV03'): a residue observed but not pydssp-H counts as missing, so a helix is 'complete' (validation) only when every position of its consensus range is observed and helical, and 'observed length' = helical C-alpha used. (3) H12 = 538-547 (primary window, run file section 2), not the SV02 proposal 538-544: only 75 chains have all of 538-547 helical (cl7 58, cl4 17), so the H12 error curve rests on them, and most chains' H12 'missing' residues are 546-547 observed but not helical, not unobserved. (4) [DEFAULT] error per observed length = median over N- and C-end deletions pooled (per-end curves are in validation_error.csv). (5) 'Report the axis angle separately' = RMS axis angle over the shared helices, clustered the same way with its own ARI rows, plus per-helix distance_summary.csv. (6) K = 2: linkage recomputed on the cl7 + cl4 sub-matrix; K = 7: all chains cut into 7 clusters vs the S08 1.5 A labels; entry level = the K07 chain per entry (results/entry_level_clusters.csv), sub-matrix re-linked; ARI is descriptive (no test, so no BH). (7) Step 5 'where they land' made operational three ways, fixed before the run: landing_h12 (primary; nearest of the cl7 / cl4 reference H12 midpoints = coordinate-wise median over complete-H12 chains), landing_k2 and landing_k7 (SV04 cluster at K = 2 / 7, labelled by the majority S08 group of its complete-H12 members); minor-cluster chains cannot agree with a cl7 / cl4 label, so their agreement is n/a. (8) Figure panel (b) is not in the step text: the SV04 H12 midpoint (axis point of 542.5) and the SV03 H12 centroid differ along the H12 axis, which moves x by a cluster-dependent amount (complete cl4 median +1.51 A, cl7 -0.68 A; y identical), so panel (a) as specified (complete chains at SV03 positions, partial at SV04 estimates) carries that method offset between filled and hollow markers; (b) draws every chain with the SV04 estimator. (9) Supporting files not named in the run file: numbers_midpoints.json, numbers_cluster.json, cluster_labels.csv, cluster_crosstab.csv, distance_summary.csv, h12_map_coords.csv, partial_h12_entry_level.csv, data/sv/SV04/, scripts/sv/common/helix_midpoint.py.

# SV04 - Helix midpoints for incomplete helices (Ye)

**Verdict:** clustering chains on helix midpoints reproduces S08: ARI 0.954 at K = 2 (cl7 + cl4, 774 chains) and 0.931 at K = 7 (all 788 chains); entry level 0.977 / 0.959. All 9 chain-level K = 2 mismatches (cl4 5, cl7 4) are chains without a usable H12 midpoint (9 of 36 cl7 / cl4 chains without one, clustered on the other 7 helices). The midpoint distance is dominated by H12 (cl7 vs cl4 H12 midpoints median 12.6 A apart; every other helix <= 0.30 A), so the agreement mostly restates the H12 split. Partial-H12 chains: 199 of 212 have an H12 midpoint usable by the [DEFAULT] rule; every cl7 / cl4 one lands on its own S08 side by H12 position (191 of 191). Partial chains without a usable H12 midpoint are placed by the other helices only and agree less (K = 2: cl4 7 of 8, cl7 1 of 3; disagreeing: 7UJW|C|1 (cl4, 6 of 10 H12 C-alpha observed -> cl7-like), 9VT5|B|1 (cl7, 5 of 10 H12 C-alpha observed -> cl4-like), 9VTU|B|1 (cl7, 2 of 10 H12 C-alpha observed -> cl4-like)). Agreement with S08 is not independent evidence: the S08 RMSD already used the same observed H12 C-alpha.

## Steps 1-3 - midpoints and truncation validation (validation_error.csv, midpoints.csv)

- Frame / axes: core v3 (SV03 frozen frame, load_ca_core); SV04 fit sets equal SV03's for 6304 of 6304 chain x helix rows, axis status equal 6304/6304, max axis angle difference 1.5e-06 deg, every midpoint lies on the SV03 axis (max offset 1.1e-14 A). H12 midpoint minus SV03 H12 centroid along the axis: median +1.48 A (max |3.43| A); other helices median within +-0.03 A.
- Rise per residue on complete helices (median): H1 1.50, H3 1.50, H4/5 1.45, H7 1.50, H8 1.52, H9 1.50, H10/11 1.53, H12 1.55 A.
- Complete helices used for validation (all positions observed and helical): H1 767 (cl7 435, cl4 318, minor 14), H3 773 (cl7 440, cl4 319, minor 14), H4/5 756 (cl7 427, cl4 315, minor 14), H7 756 (cl7 436, cl4 306, minor 14), H8 778 (cl7 442, cl4 323, minor 13), H9 697 (cl7 368, cl4 315, minor 14), H10/11 645 (cl7 402, cl4 229, minor 14), H12 75 (cl7 58, cl4 17, minor 0).

Median midpoint error (A) / median axis angle error (deg) after deleting k residues from either end (N and C pooled; n = truncations estimated); * = fails the [DEFAULT] (median > 1.0 A); - = fewer than 6 C-alpha left:

| helix | range | L | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | eligible k |
|---|---|---|---|---|---|---|---|---|---|
| H1 | 312-321 | 10 | 0.04 / 1.1 (n=1534) | 0.07 / 1.5 (n=1534) | 0.17 / 2.9 (n=1534) | 0.43 / 6.7 (n=1534) | - | - | 0,1,2,3,4 |
| H3 | 342-361 | 20 | 0.10 / 1.2 (n=1546) | 0.22 / 2.5 (n=1546) | 0.32 / 3.6 (n=1546) | 0.43 / 4.7 (n=1546) | 0.61 / 6.1 (n=1546) | 0.80 / 7.6 (n=1546) | 0,1,2,3,4,5,6 |
| H4/5 | 372-394 | 23 | 0.18 / 2.0 (n=1512) | 0.39 / 4.1 (n=1512) | 0.62 / 6.3 (n=1512) | 0.91 / 8.7 (n=1512) | 1.26 / 11.4 * (n=1512) | 1.61 / 14.1 * (n=1512) | 0,1,2,3,4 |
| H7 | 422-437 | 16 | 0.04 / 0.6 (n=1512) | 0.07 / 1.1 (n=1512) | 0.07 / 1.1 (n=1512) | 0.06 / 0.9 (n=1512) | 0.16 / 2.0 (n=1512) | 0.26 / 2.8 (n=1512) | 0,1,2,3,4,5,6 |
| H8 | 442-455 | 14 | 0.05 / 0.7 (n=1556) | 0.07 / 1.0 (n=1556) | 0.09 / 1.3 (n=1556) | 0.13 / 1.5 (n=1556) | 0.17 / 1.9 (n=1556) | 0.24 / 2.7 (n=1556) | 0,1,2,3,4,5,6 |
| H9 | 473-492 | 20 | 0.02 / 0.3 (n=1394) | 0.04 / 0.5 (n=1394) | 0.04 / 0.4 (n=1394) | 0.03 / 0.3 (n=1394) | 0.06 / 0.6 (n=1394) | 0.09 / 0.8 (n=1394) | 0,1,2,3,4,5,6 |
| H10/11 | 497-526 | 30 | 0.04 / 0.3 (n=1290) | 0.08 / 0.6 (n=1290) | 0.12 / 0.9 (n=1290) | 0.18 / 1.3 (n=1290) | 0.23 / 1.6 (n=1290) | 0.28 / 1.9 (n=1290) | 0,1,2,3,4,5,6 |
| H12 | 538-547 | 10 | 0.08 / 1.9 (n=150) | 0.12 / 1.7 (n=150) | 0.19 / 3.0 (n=150) | 0.54 / 8.0 (n=150) | - | - | 0,1,2,3,4 |

- [DEFAULT] applied: eligible chain x helix H1 788/788, H3 788/788, H4/5 788/788, H7 788/788, H8 788/788, H9 788/788, H10/11 788/788, H12 748/788. Every chain x helix with an axis falls on a validated observed length that passes (H4/5 fails at k >= 5, but no chain misses more than 4 H4/5 residues; H12 k distribution k=0: 75, k=1: 89, k=2: 347, k=3: 196, k=4: 41, k=5: 9, k=6: 4, k=10: 27); the only exclusions are the 40 H12 helices with < 6 observed helical C-alpha (no axis, as in SV03). Eligible helices per chain: 7 helices 40 chains, 8 helices 748 chains.
- H12 by SV01a status: complete 553 chains, eligible 549 (cl7 427/431, cl4 120/120, minor 2/2); not_modelled 23 chains, eligible 0 (cl7 0/7, cl4 0/14, minor 0/2); partial 212 chains, eligible 199 (cl7 4/7, cl4 187/195, minor 8/10).

## Step 4 - clustering on midpoints (ari.csv, cluster_crosstab.csv, distance_summary.csv)

- Distance = RMS of midpoint distances over shared eligible helices; pairs share 7 helices: 30700, 8 helices: 279378 pairs; no chain dropped (0). Average linkage.

| distance | level | K | chains | n | SV04 cluster sizes | ARI vs S08 |
|---|---|---|---|---|---|---|
| midpoint_rms | chain | 2 | cl7 + cl4 chains | 774 | 446;328 | 0.954 |
| midpoint_rms | chain | 7 | all chains | 788 | 448;333;2;2;1;1;1 | 0.931 |
| midpoint_rms | entry | 2 | cl7 + cl4 chains | 344 | 219;125 | 0.977 |
| midpoint_rms | entry | 7 | all chains | 351 | 220;124;2;2;1;1;1 | 0.959 |
| axis_angle_rms | chain | 2 | cl7 + cl4 chains | 774 | 453;321 | 0.949 |
| axis_angle_rms | chain | 7 | all chains | 788 | 454;306;18;5;3;1;1 | 0.910 |
| axis_angle_rms | entry | 2 | cl7 + cl4 chains | 344 | 220;124 | 0.965 |
| axis_angle_rms | entry | 7 | all chains | 351 | 221;123;2;2;1;1;1 | 0.955 |

Per-helix pair distances (chain level, eligible helices; median [90th percentile]):

| helix | midpoint within cl7 (A) | within cl4 | cl7 vs cl4 | axis angle within cl7 (deg) | within cl4 | cl7 vs cl4 |
|---|---|---|---|---|---|---|
| H1 | 0.21 [0.40] | 0.35 [0.63] | 0.30 [0.54] | 1.4 [2.8] | 1.9 [3.5] | 1.8 [3.3] |
| H3 | 0.20 [0.34] | 0.33 [0.68] | 0.30 [0.53] | 1.0 [2.1] | 1.7 [3.0] | 1.6 [2.7] |
| H4/5 | 0.12 [0.23] | 0.22 [0.41] | 0.18 [0.34] | 0.6 [1.4] | 1.0 [2.2] | 1.0 [1.9] |
| H7 | 0.17 [0.31] | 0.28 [0.54] | 0.23 [0.45] | 1.2 [2.4] | 1.4 [2.6] | 1.3 [2.5] |
| H8 | 0.09 [0.16] | 0.12 [0.22] | 0.12 [0.20] | 0.8 [1.5] | 1.0 [2.2] | 0.9 [1.9] |
| H9 | 0.21 [0.36] | 0.28 [0.58] | 0.25 [0.48] | 0.8 [1.7] | 1.1 [2.2] | 1.2 [2.2] |
| H10/11 | 0.13 [0.25] | 0.19 [0.33] | 0.19 [0.32] | 0.6 [1.2] | 1.0 [1.9] | 1.2 [2.0] |
| H12 | 0.36 [0.78] | 0.95 [1.93] | 12.58 [13.44] | 3.3 [6.3] | 6.4 [15.0] | 118.6 [124.7] |
| RMS over shared helices | 0.23 [0.36] | 0.44 [0.76] | 4.43 [4.74] | 1.6 [2.6] | 2.6 [5.3] | 41.8 [44.0] |

## Step 5 - partial-H12 chains (partial_h12_assignment.csv, partial_h12_entry_level.csv)

- Partial H12 (SV01a: 1-9 of the 10 C-alpha 538-547 observed): 212 chains (cl4 195, minor 10, cl7 7); usable H12 midpoint 199; placed on the map (H12, H3, H10/11 usable) 199 (cl4 187, minor 8, cl7 4).
- Reference H12 midpoints (complete-H12 chains with a usable H12): cl7 n = 427, cl4 n = 120; 12.3 A apart. Margin d(cl4 ref) - d(cl7 ref) for placed partial chains: cl4 median -12.1 A (range -12.3 to -9.3); cl7 median +12.2 A (range +10.9 to +12.3); minor median +4.8 A (range -12.3 to +7.7).
- SV04 cluster labels (majority S08 group of complete-H12 members): K = 2 {'1': 'cl7', '2': 'cl4'}; K = 7 {'2': 'none', '1': 'none', '5': 'cl7', '4': 'cl4', '3': 'none', '6': 'none', '7': 'none'}.

| S08 group | partial chains | landing_h12 (H12 midpoint) | agree | landing_k2 | agree | landing_k7 | agree |
|---|---|---|---|---|---|---|---|
| cl7 | 7 | cl7 4, not placed 3 | 4/4 | cl7 5, cl4 2 | 5/7 | cl7 5, cl4 2 | 5/7 |
| cl4 | 195 | cl4 187, not placed 8 | 187/187 | cl4 194, cl7 1 | 194/195 | cl4 192, none 2, cl7 1 | 192/195 |
| minor | 10 | cl7 5, cl4 3, not placed 2 | n/a | n/a (minor cluster) 10 | n/a | none 5, cl4 5 | 0/10 |

- Agreement split by whether the chain has a usable H12 midpoint (agree / chains): agree_h12|cl4, H12 usable True 187/187; agree_h12|cl7, H12 usable True 4/4; agree_k2|cl4, H12 usable False 7/8; agree_k2|cl4, H12 usable True 187/187; agree_k2|cl7, H12 usable False 1/3; agree_k2|cl7, H12 usable True 4/4; agree_k7|cl4, H12 usable False 7/8; agree_k7|cl4, H12 usable True 185/187; agree_k7|cl7, H12 usable False 1/3; agree_k7|cl7, H12 usable True 4/4.
- Minor-cluster partial chains with an H12 midpoint (distance to the cl7 / cl4 reference, A): 1A52|A|1 (S08 cl1) 31.0 / 37.4 -> cl7-like; 1A52|B|1 (S08 cl1) 31.0 / 38.3 -> cl7-like; 1R5K|A|1 (S08 cl2) 9.3 / 4.0 -> cl4-like; 6OWC|B|1 (S08 cl3) 14.5 / 2.2 -> cl4-like; 6PFM|A|1 (S08 cl1) 29.8 / 37.5 -> cl7-like; 9B25|A|1 (S08 cl6) 7.3 / 11.6 -> cl7-like; 9B25|C|1 (S08 cl6) 10.2 / 15.4 -> cl7-like; 9B25|D|1 (S08 cl3) 14.2 / 2.4 -> cl4-like.
- Entry level (majority landing vs majority S08 group of an entry's landed partial chains): landing_h12: 83 agree, 0 disagree, 1 tie, 4 n/a (minor) of 88 entries; landing_k2: 88 agree, 2 disagree, 0 tie, 0 n/a (minor) of 90 entries; landing_k7: 85 agree, 7 disagree, 2 tie, 0 n/a (minor) of 94 entries.
- Map positions of the 199 placed partial chains, SV04 estimate minus SV03 position: x median +1.46 A (max |2.72|), y max |2.5e-14| deg; every placed partial chain already had an SV03 point (0 new). Same difference for complete-H12 chains: cl4 +1.51 A (n=120), cl7 -0.68 A (n=427), minor +0.98 A (n=2); the x offset is a method difference (midpoint at residue 542.5 vs centroid of the observed helical windows), not a position of the partial chains.

## Figure

- figures/sv/SV04/h12_map_partial.png/.svg: SV03 H12 map axes (x = H12 point to H11 axis, y = angle H12-H3; colour = S08 cluster, marker = ligand_chemotype). (a) as specified: 549 complete-H12 chains filled at SV03 positions, 199 partial-H12 chains hollow at SV04 estimates; (b) the same chains, all at SV04 estimates. Not drawn: 40 chains (not_modelled H12, or no usable H12 / H3 / H10/11).

## Files

- results/sv/SV04/ari.csv (0.9 kB)
- results/sv/SV04/cluster_crosstab.csv (1.9 kB)
- results/sv/SV04/cluster_labels.csv (39.8 kB)
- results/sv/SV04/distance_summary.csv (3.0 kB)
- results/sv/SV04/env.txt (0.5 kB)
- results/sv/SV04/h12_map_coords.csv (205.9 kB)
- results/sv/SV04/midpoints.csv (2.32 MB)
- results/sv/SV04/numbers_cluster.json (5.9 kB)
- results/sv/SV04/numbers_midpoints.json (5.6 kB)
- results/sv/SV04/partial_h12_assignment.csv (50.2 kB)
- results/sv/SV04/partial_h12_entry_level.csv (8.6 kB)
- results/sv/SV04/validation_error.csv (33.0 kB)
- data/sv/SV04/axis_angle_rms.npy (4.97 MB)
- data/sv/SV04/matrix_index.csv (14.1 kB)
- data/sv/SV04/midpoint_rms.npy (4.97 MB)
- data/sv/SV04/n_shared.npy (4.97 MB)
- data/sv/SV04/validation_truncations.csv.gz (1.95 MB)
- figures/sv/SV04/h12_map_partial.png (363.6 kB)
- figures/sv/SV04/h12_map_partial.svg (323.4 kB)
- figures/panels/sv/SV04/h12_map_partial__a_sv03_positions.png (193.2 kB)
- figures/panels/sv/SV04/h12_map_partial__b_sv04_all.png (166.5 kB)
- scripts/sv/SV04_midpoints.py (15.6 kB)
- scripts/sv/SV04_cluster.py (22.0 kB)
- scripts/sv/SV04_figures.py (4.9 kB)
- scripts/sv/SV04_summary.py (17.9 kB)
- scripts/sv/common/helix_midpoint.py (2.2 kB)
- results/sv/SV04/SUMMARY.md (this file)
