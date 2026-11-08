Not done / substituted: (1) 'helical' and helix ends come from SV02's pydssp fallback (mkdssp not installed; pydssp 'H' = any helix, 3-, 4- and 5-turn together), not from DSSP alpha-H. (2) H12 range = 538-547, the primary H12 window of run file section 2, not the SV02 data-derived consensus H12 538-544 (a proposal); the H12 descriptors with 538-544 are kept as a descriptive sensitivity copy (geometry_per_chain.csv *_w544, h12_window_sensitivity.csv), not tested. (3) 'H11' = SV02 consensus helix H10/11 497-526 (SV02 keeps H10/H11 combined); its halves split that range at the midpoint (H11N 497-511, H11C 512-526); 'H11-H12 hinge angle' = angle H11C-H12 (the half next to the hinge), distinct from angle H12-H11 (whole H11) - an interpretation, fixed before the run. (4) A third terminus state the step does not name: 'gap' = next residue not helical but the second one helical (1-residue break before another H run): neither censored nor a 'real end', so not an uncensored length (H7 2 start / 0 end, H8 0 start / 14 end, H9 1 start / 0 end, H10/11 0 start / 1 end). (5) Step-5 [DEFAULT] for the one continuous covariate (ligand heavy atoms): Spearman |rho| >= 0.3 used in place of |rank-biserial r|. (6) 'mutations' coded as in SV01a: Y537S, L536S, D538G, other construct mutation (each present vs absent). (7) Chain-level rank test = GEE on average ranks with entry clusters (the SV11 convention); entry level = Mann-Whitney U / Spearman on entry medians; BH within four families ((a) chain, (a) entry, (b) chain, (b) entry). Tests with < 2 entries in a group were not run; a group of 2-5 entries passes this rule, and cluster-robust GEE errors are optimistic with so few clusters: 4 of the 37 [DEFAULT] hits rest on such a group (flagged below). (8) Frame = SV00 gate transforms (S03 code path, equal to results/alignment/ within 3.6e-14 A per SV00); frozen copy data/sv/SV03/core_v3_transforms.csv. (9) Supporting files not named in the run file: helix_per_chain.csv, helix_length_summary.csv, h12_window_sensitivity.csv, axis_method_check.csv, sensitivity_check.csv, cl7_vs_cl4_wide.csv, frame_check.json, numbers_geometry.json, numbers_tests.json, data/sv/SV03/, scripts/sv/common/{core_frame,helix_geometry,h12_map}.py.

# SV03 - Helix geometry

**Verdict ([DEFAULT] step 5):** 37 within-cluster comparisons meet |r| >= 0.3 and q < 0.05 at both chain and entry level (19 in cl7, 18 in cl4; 21 use H12, 16 only other helices), so by the rule the geometric descriptors are **more sensitive than RMSD** for those subgroups; 4 of them rest on a group of <= 5 entries. cl7 vs cl4 (step 4a): 12 of 16 descriptors differ at q < 0.05 at both levels, all 7 H12 descriptors among them; not separated: length H1, length H4/5, length H7, length H8.

Size of the within-cluster differences (binary hits, |median yes - median no| at chain level): angles incl. kink and hinge 0.5-5.4 deg, distances 0.52-0.70 A, lengths 1-4 residues; for comparison the cl7 - cl4 chain medians differ by 32.4 deg (angle H12-H3) and 9.8 A (H12 centroid-H11 axis) (step 4a table).

## Frame (step inputs)

- Core v3 frame of the SV00 gate (reference 7NFB|B|1, 788 chains, status {'ok': 788}); convention `x_core = (x_deposited - mobile_centroid) @ R + reference_centroid (row vectors)`; deposited backbone C-alpha transformed with the frozen transforms vs data/sv/gate/aligned_coordinates.npz: max |d| = 2.1e-14 A, same observed mask True.
- Reuse (SV04, SV13): data/sv/SV03/core_v3_transforms.csv + scripts/sv/common/core_frame.py (load_transforms, to_core, load_ca_core); axis fit scripts/sv/common/helix_geometry.py (fit_axis); H12 map scripts/sv/common/h12_map.py (plot_h12_map, add_legends; marker column and hollow markers are parameters).
- Run ends and censoring reproduce SV02 sse_per_chain.csv exactly for every helix with the same range (status, start, end, start/end censored: H1 788/788; H3 788/788; H4/5 788/788; H7 788/788; H8 788/788; H9 788/788; H10/11 788/788; H12_w544 761/761).

## Step 1 - helix lengths and censoring (helix_length_summary.csv)

Length = first to last pydssp-H residue of the chain's runs overlapping the range; censored = a terminus within 2 residues of an unobserved residue or chain break; uncensored = both termini real ends.

| helix | range | chains with a run | start cens. | end cens. | gap ends | censoring fraction (all / cl7 / cl4) | uncensored | median [IQR] uncensored length cl7 | cl4 |
|---|---|---|---|---|---|---|---|---|---|
| H1 | 312-321 | 788 | 4 | 0 | 0 | 0.005 / 0.000 / 0.012 | 784 | 10 [10-10] (n=445) | 10 [10-10] (n=325) |
| H3 | 342-361 | 788 | 189 | 0 | 0 | 0.240 / 0.094 / 0.447 | 599 | 23 [23-25] (n=403) | 23 [22-25] (n=182) |
| H4/5 | 372-394 | 788 | 2 | 1 | 0 | 0.004 / 0.004 / 0.003 | 785 | 23 [23-23] (n=443) | 23 [23-23] (n=328) |
| H7 | 422-437 | 788 | 116 | 0 | 2 | 0.147 / 0.097 / 0.222 | 670 | 17 [16-17] (n=402) | 17 [17-18] (n=254) |
| H8 | 442-455 | 788 | 0 | 2 | 14 | 0.003 / 0.002 / 0.003 | 772 | 14 [14-14] (n=439) | 14 [14-14] (n=319) |
| H9 | 473-492 | 788 | 434 | 0 | 1 | 0.551 / 0.598 / 0.502 | 353 | 24 [21-26] (n=178) | 27 [27-28] (n=164) |
| H10/11 | 497-526 | 788 | 0 | 172 | 1 | 0.218 / 0.079 / 0.413 | 615 | 34 [33-35] (n=409) | 31 [30-31] (n=193) |
| H12 | 538-547 | 761 | 55 | 223 | 0 | 0.330 / 0.140 / 0.589 | 510 | 8 [8-8] (n=374) | 8 [8-9] (n=129) |

- H12 538-547: 23 chains observe none of it, 4 have no H residue in it; 761 have a run.

## Step 2-3 - axes and descriptors (geometry_per_chain.csv, helix_per_chain.csv)

- Axis [DEFAULT]: >= 6 observed helical C-alpha in the range, centroids of consecutive 4-residue windows, principal line, oriented N -> C; helix centroid = mean window centroid (core v3 frame, centroid_<helix>_x/y/z).
- Axes per helix: H1 788/788; H10/11 788/788; H11C 788/788; H11N 788/788; H12 748/788; H3 788/788; H4/5 788/788; H7 788/788; H8 788/788; H9 788/788. H12 without an axis: 40 chains (< 6 observed helical C-alpha in 538-547).
- Axis rule on ideal alpha-helices (axis_method_check.csv): angle error 6 CA 9.6 deg, 7 CA 4.5 deg, 8 CA 0.6 deg, 9 CA 1.4 deg, 10 CA 1.6 deg, 11 CA 0.9 deg, 12 CA 0.1 deg; H12 axes use cl4: 6 CA 34, 7 CA 166, 8 CA 84, 9 CA 6, 10 CA 17; cl7: 6 CA 1, 7 CA 28, 8 CA 262, 9 CA 82, 10 CA 58; minor: 6 CA 6, 7 CA 2, 8 CA 1, 9 CA 1 chains - cl4 H12 axes rest mostly on 7 C-alpha, cl7 on 8, so H12 angle differences of a few degrees inside a cluster are within the method's range for short H12 runs.
- H12 map points (chains with an H12 axis): 748 (cl7 431, cl4 307, minor 10).

## Step 4a - cl7 vs cl4 (tests.csv part a, cl7_vs_cl4_wide.csv)

Chain level 445 cl7 / 329 cl4 chains at most; entry level entry groups by the SV10 rule ({'cl7': 219, 'cl4': 124, 'mixed': 8}). r = rank-biserial, positive = higher in cl4.

| descriptor | uses H12 | median cl7 | median cl4 | r chain | q chain | r entry | q entry | q < 0.05 both |
|---|---|---|---|---|---|---|---|---|
| length H1 (res) | no | 10.00 | 10.00 | -0.02 | 0.667 | -0.01 | 0.781 | no |
| length H3 (res) | no | 23.00 | 23.00 | -0.26 | 1.5e-04 | -0.31 | 1.0e-05 | yes |
| length H4/5 (res) | no | 23.00 | 23.00 | +0.05 | 0.176 | +0.03 | 0.624 | no |
| length H7 (res) | no | 17.00 | 17.00 | +0.10 | 0.065 | +0.18 | 0.005 | no |
| length H8 (res) | no | 14.00 | 14.00 | +0.02 | 0.271 | +0.01 | 0.618 | no |
| length H9 (res) | no | 24.00 | 27.00 | +0.67 | 1.2e-24 | +0.69 | 2.1e-18 | yes |
| length H10/11 (res) | no | 34.00 | 31.00 | -0.82 | 2.2e-115 | -0.84 | 8.7e-33 | yes |
| length H12 (res) | yes | 8.00 | 8.00 | +0.13 | 0.028 | +0.19 | 0.006 | yes |
| angle H12-H3 (deg) | yes | 79.05 | 46.63 | -1.00 | 7.0e-299 | -1.00 | 1.2e-51 | yes |
| angle H12-H11 (deg) | yes | 108.72 | 126.67 | +1.00 | 2.9e-270 | +1.00 | 1.2e-51 | yes |
| angle H12-H4/5 (deg) | yes | 104.70 | 90.50 | -0.98 | 3.8e-237 | -1.00 | 1.7e-51 | yes |
| angle H11-H3 (deg) | no | 128.73 | 130.12 | +0.66 | 1.9e-37 | +0.67 | 1.9e-24 | yes |
| H12 centroid-H3 axis (A) | yes | 11.39 | 10.82 | -0.78 | 6.6e-92 | -0.84 | 3.5e-37 | yes |
| H12 centroid-H11 axis (A) | yes | 12.48 | 22.26 | +1.00 | 1.2e-259 | +1.00 | 1.2e-51 | yes |
| H11 kink (deg) | no | 10.16 | 8.84 | -0.58 | 2.4e-29 | -0.66 | 3.2e-24 | yes |
| H11-H12 hinge (deg) | yes | 115.87 | 120.56 | +0.73 | 1.0e-57 | +0.81 | 3.0e-34 | yes |

## Step 4b-5 - within cl7 and within cl4 (tests.csv part b, sensitivity_check.csv)

- Comparisons: 16 descriptors x 8 covariates x 2 clusters = 256; tested at both levels 208 (chain-level tests 208/256, entry-level 208/256); not tested (< 2 chains or entries in a group, per descriptor): cl4 D538G 16 descriptors; cl4 peptide 8 descriptors; cl7 basic amine 7 descriptors; cl7 D538G 1 descriptors; cl7 L536S 16 descriptors.
- q < 0.05 at both levels: 43; of these, |effect| >= 0.3 at both levels ([DEFAULT] met): 37; sign equal at both levels for 37/37.
- Hits by covariate: Y537S 9, carboxylic acid 5, other mutation 5, basic amine 5, L536S 5, heavy atoms 4, peptide 2, D538G 2.

| cluster | descriptor | covariate | chains (no / yes) | entries (no / yes) | median no / yes (chain) | r chain | q chain | r entry | q entry | small group |
|---|---|---|---|---|---|---|---|---|---|---|
| cl7 | H11 kink | heavy atoms | 438 (rho) | 219 | - | -0.54 | 1.6e-33 | -0.63 | 9.4e-24 |  |
| cl7 | H11 kink | carboxylic acid | 420 / 18 | 209 / 10 | 10.22 / 8.64 deg | -0.60 | 1.8e-06 | -0.57 | 0.014 |  |
| cl7 | H11 kink | Y537S | 40 / 405 | 21 / 200 | 11.02 / 10.12 deg | -0.40 | 0.005 | -0.35 | 0.040 |  |
| cl7 | length H1 | Y537S | 40 / 405 | 21 / 200 | 11.00 / 10.00 res | -0.45 | 4.2e-05 | -0.34 | 0.014 |  |
| cl7 | length H3 | carboxylic acid | 382 / 15 | 202 / 8 | 23.00 / 25.00 res | +0.56 | 9.6e-13 | +0.70 | 0.006 |  |
| cl7 | length H3 | Y537S | 33 / 370 | 16 / 197 | 24.00 / 23.00 res | -0.31 | 8.9e-05 | -0.43 | 0.017 |  |
| cl7 | length H9 | Y537S | 19 / 159 | 14 / 134 | 27.00 / 23.00 res | -0.50 | 5.0e-04 | -0.56 | 0.006 |  |
| cl7 | length H9 | other mutation | 171 / 7 | 143 / 5 | 23.00 / 27.00 res | +0.77 | 5.2e-45 | +0.79 | 0.015 | yes |
| cl7 | angle H12-H11 | carboxylic acid | 409 / 15 | 203 / 8 | 108.78 / 107.03 deg | -0.55 | 1.3e-05 | -0.66 | 0.007 |  |
| cl7 | angle H12-H11 | peptide | 25 / 406 | 10 / 206 | 106.69 / 108.79 deg | +0.53 | 0.001 | +0.52 | 0.030 |  |
| cl7 | angle H12-H11 | Y537S | 34 / 397 | 17 / 199 | 104.79 / 108.91 deg | +0.91 | 4.5e-82 | +0.94 | 8.1e-09 |  |
| cl7 | angle H12-H11 | D538G | 427 / 4 | 214 / 2 | 108.72 / 103.29 deg | -0.97 | 1.7e-183 | -0.98 | 0.007 | yes |
| cl7 | angle H12-H3 | Y537S | 34 / 397 | 17 / 199 | 80.94 / 78.91 deg | -0.56 | 1.4e-06 | -0.57 | 0.002 |  |
| cl7 | angle H12-H4/5 | Y537S | 34 / 397 | 17 / 199 | 107.16 / 104.47 deg | -0.61 | 1.2e-12 | -0.58 | 9.5e-04 |  |
| cl7 | H12 centroid-H3 axis | Y537S | 34 / 397 | 17 / 199 | 10.71 / 11.40 A | +0.71 | 2.5e-09 | +0.91 | 1.7e-08 |  |
| cl7 | H11-H12 hinge | heavy atoms | 424 (rho) | 211 | - | -0.34 | 5.0e-09 | -0.47 | 1.0e-10 |  |
| cl7 | H11-H12 hinge | carboxylic acid | 409 / 15 | 203 / 8 | 116.11 / 113.53 deg | -0.60 | 8.7e-07 | -0.66 | 0.007 |  |
| cl7 | H11-H12 hinge | Y537S | 34 / 397 | 17 / 199 | 113.24 / 116.23 deg | +0.73 | 2.5e-36 | +0.78 | 2.1e-06 |  |
| cl7 | H11-H12 hinge | D538G | 427 / 4 | 214 / 2 | 115.90 / 111.52 deg | -0.91 | 1.9e-152 | -0.97 | 0.012 | yes |
| cl4 | angle H11-H3 | heavy atoms | 329 (rho) | 128 | - | +0.48 | 4.6e-14 | +0.55 | 1.3e-09 |  |
| cl4 | angle H11-H3 | peptide | 325 / 4 | 125 / 3 | 130.12 / 127.92 deg | -0.96 | 2.4e-80 | -0.99 | 0.019 | yes |
| cl4 | H11 kink | basic amine | 107 / 222 | 36 / 96 | 8.52 / 9.02 deg | +0.32 | 0.005 | +0.37 | 0.013 |  |
| cl4 | length H10/11 | L536S | 52 / 141 | 33 / 65 | 30.00 / 31.00 res | +0.41 | 2.9e-05 | +0.36 | 0.017 |  |
| cl4 | length H10/11 | other mutation | 26 / 167 | 16 / 82 | 29.00 / 31.00 res | +0.43 | 0.006 | +0.46 | 0.016 |  |
| cl4 | length H3 | basic amine | 64 / 118 | 33 / 75 | 24.00 / 22.00 res | -0.35 | 0.001 | -0.38 | 0.013 |  |
| cl4 | length H3 | L536S | 57 / 125 | 37 / 68 | 22.00 / 23.00 res | +0.49 | 1.2e-05 | +0.54 | 7.1e-05 |  |
| cl4 | length H3 | other mutation | 30 / 152 | 18 / 87 | 22.00 / 23.00 res | +0.46 | 8.6e-04 | +0.46 | 0.012 |  |
| cl4 | angle H12-H11 | basic amine | 95 / 212 | 34 / 93 | 125.14 / 127.16 deg | +0.39 | 4.7e-06 | +0.51 | 4.8e-04 |  |
| cl4 | angle H12-H3 | L536S | 77 / 230 | 41 / 82 | 42.63 / 47.29 deg | +0.60 | 2.9e-13 | +0.61 | 8.4e-07 |  |
| cl4 | angle H12-H3 | other mutation | 30 / 277 | 17 / 106 | 42.23 / 46.90 deg | +0.58 | 1.0e-06 | +0.55 | 0.003 |  |
| cl4 | angle H12-H4/5 | basic amine | 95 / 212 | 34 / 93 | 91.03 / 89.32 deg | -0.31 | 4.4e-04 | -0.40 | 0.008 |  |
| cl4 | angle H12-H4/5 | L536S | 77 / 230 | 41 / 82 | 92.29 / 89.22 deg | -0.63 | 3.9e-13 | -0.69 | 1.7e-08 |  |
| cl4 | angle H12-H4/5 | other mutation | 30 / 277 | 17 / 106 | 92.64 / 90.09 deg | -0.41 | 0.020 | -0.48 | 0.011 |  |
| cl4 | H12 centroid-H11 axis | heavy atoms | 307 (rho) | 123 | - | +0.38 | 3.7e-05 | +0.48 | 6.2e-07 |  |
| cl4 | H12 centroid-H11 axis | basic amine | 95 / 212 | 34 / 93 | 22.59 / 22.08 A | -0.30 | 0.006 | -0.31 | 0.045 |  |
| cl4 | length H12 | carboxylic acid | 116 / 13 | 74 / 10 | 8.00 / 11.00 res | +0.40 | 0.038 | +0.55 | 0.012 |  |
| cl4 | length H12 | L536S | 30 / 99 | 27 / 57 | 9.00 / 8.00 res | -0.45 | 7.8e-04 | -0.45 | 0.003 |  |

r = rank-biserial (positive = higher with the covariate), rho = Spearman for heavy atoms; small group = <= 5 entries in the smaller group (a flag for reading, not part of the rule).

## H12 window sensitivity (h12_window_sensitivity.csv, descriptive)

| descriptor | chains | Spearman rho 538-547 vs 538-544 | median abs diff | max abs diff |
|---|---|---|---|---|
| angle_H12_H3_deg | 748 | 0.981 | 0.35 | 7.80 |
| angle_H12_H11_deg | 748 | 0.969 | 2.38 | 7.39 |
| angle_H12_H45_deg | 748 | 0.966 | 3.04 | 7.73 |
| dist_H12c_H3axis_A | 748 | 0.455 | 0.41 | 1.51 |
| dist_H12c_H11axis_A | 748 | 0.904 | 0.37 | 1.82 |
| hinge_H11C_H12_deg | 748 | 0.924 | 2.36 | 7.45 |
| len_H12 | 761 | 1.000 | 0.00 | 0.00 |
| ulen_H12 | 510 | 1.000 | 0.00 | 0.00 |

Run ends do not depend on the window here: every H12 run overlapping 538-547 also overlaps 538-544, so lengths are identical; axes change because 545-547 add C-alpha.

## Figures

- figures/sv/SV03/SV03_h12_map.png/.svg: x = H12 centroid to H11 axis, y = angle H12-H3, colour = S08 cluster, marker = ligand_chemotype, one point per chain with an H12 axis.
- figures/sv/SV03/SV03_descriptors_by_cluster.png/.svg: the 16 tested descriptors by cluster group with the step-4a r and q per panel.

## Files

- results/sv/SV03/axis_method_check.csv (0.5 kB)
- results/sv/SV03/cl7_vs_cl4_wide.csv (2.8 kB)
- results/sv/SV03/env.txt (0.4 kB)
- results/sv/SV03/frame_check.json (1.0 kB)
- results/sv/SV03/geometry_per_chain.csv (896.7 kB)
- results/sv/SV03/h12_window_sensitivity.csv (0.7 kB)
- results/sv/SV03/helix_length_summary.csv (2.1 kB)
- results/sv/SV03/helix_per_chain.csv (2.13 MB)
- results/sv/SV03/numbers_geometry.json (3.5 kB)
- results/sv/SV03/numbers_tests.json (12.3 kB)
- results/sv/SV03/sensitivity_check.csv (67.3 kB)
- results/sv/SV03/tests.csv (192.7 kB)
- data/sv/SV03/core_v3_transforms.csv (257.6 kB)
- figures/sv/SV03/SV03_descriptors_by_cluster.png (741.5 kB)
- figures/sv/SV03/SV03_descriptors_by_cluster.svg (1.45 MB)
- figures/sv/SV03/SV03_h12_map.png (155.9 kB)
- figures/sv/SV03/SV03_h12_map.svg (151.6 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__angle_H11_H3_deg.png (69.2 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__angle_H12_H11_deg.png (56.3 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__angle_H12_H3_deg.png (59.5 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__angle_H12_H45_deg.png (66.1 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__dist_H12c_H11axis_A.png (63.5 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__dist_H12c_H3axis_A.png (61.5 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__hinge_H11C_H12_deg.png (60.1 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__kink_H11_deg.png (69.1 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H1.png (45.8 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H1011.png (54.8 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H12.png (46.8 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H3.png (50.8 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H45.png (44.2 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H7.png (46.1 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H8.png (40.2 kB)
- figures/panels/sv/SV03/SV03_descriptors_by_cluster__ulen_H9.png (48.8 kB)
- figures/panels/sv/SV03/SV03_h12_map__h12_map.png (154.9 kB)
- scripts/sv/SV03_geometry.py (20.3 kB)
- scripts/sv/SV03_tests.py (15.2 kB)
- scripts/sv/SV03_figures.py (4.7 kB)
- scripts/sv/SV03_summary.py (16.9 kB)
- scripts/sv/common/core_frame.py (2.8 kB)
- scripts/sv/common/helix_geometry.py (3.7 kB)
- scripts/sv/common/h12_map.py (5.7 kB)
- results/sv/SV03/SUMMARY.md (this file)
