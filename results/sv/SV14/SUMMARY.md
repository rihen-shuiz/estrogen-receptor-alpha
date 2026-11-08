Not done / substituted: (1) pocket 100 % = median pocket of wild-type EST chains (fallback, [DEFAULT] step 5): only 3 apo entries pass (2B23, 4Q13, 9VTU; 7Y8F, 7YMK fail because a PEG (7 heavy atoms) sits in the pocket, centroid 1.6-1.9 A from the E2 centroid), fewer than 5; the apo baseline (793 A^3) is reported but not used as 100 %. (2) No pocket program was used: the grid method is coded in numpy/scipy (scripts/sv/SV14_grid.py) as the step describes; receptor atoms read from the restored mmCIFs with esr_analysis.read_structure.read_cif (highest-occupancy altloc), not from data/clean/*.pdb, which keep every altloc atom (same atom selection otherwise: own polymer only, no ligand, water, ion, peptide or other chain). (3) Kept as computed, flagged: pocket_open is True for 788/788 chains (every component borders some non-buried cavity point), so the border count, not the flag, carries the information; the >= 10-ray pocket reaches the [DEFAULT] box edge in 203 chains (612 at >= 8, 4 at >= 12), so those volumes are bounded by the box. (4) Declared before the outcomes where the run file is open: ligand instances with different alternate-location ids -> the highest-occupancy alt id (31 chains); apo check over the asymmetric unit, any residue other than the chain's own polymer or water; E2 references = wild-type chains with EST alone, complete; rays marched on the 0.4 A grid (diagonals 14 steps = 9.70 A); ligand volume counted on the extended grid; several overlapping components -> largest overlap; class = SV12 action_class.csv row of (entry, analysed CCD), class analyses on primary_set rows; entry cluster group = >= 80 % rule; chain level = bootstrap over entries (10000, seed 42); BH over all 14 tested p-values (10 correlations + 4 filling tests). (5) Supporting files not named in the run file: filling_by_group.csv, filling_entry_level.csv, numbers.json, e2_reference_pocket_lining.csv, apo_blocking_residues.csv (results/sv/SV14/), data/sv/SV14/*; scripts/sv/common/class_style.py created (did not exist). Nothing installed.

# SV14 - Ligand size vs pocket size, and pocket filling

**Verdict**

- Ligand and pocket volume go together over all chains (Spearman rho 0.44 ligand level, 0.47 entry level, q 2.0e-14 / 1.4e-18) and within cl7 (0.41 / 0.46), not within cl4 (0.11 / 0.14, q >= 0.20) nor within agonist or SERM (q >= 0.15).
- Filling (median 33.8 %, IQR 24.8-39.2) does not differ cl7 vs cl4 (entry medians 33.0 vs 32.3 %, q 0.14; chain bootstrap q 0.19) or agonist vs SERM (q 0.94 / 0.60).
- References: ligand 100 % = E2 224 A^3; pocket 100 % = wild-type E2 pocket 956 A^3 (fallback: 3 apo entries pass, < 5); the E2 pocket is bimodal (505-579 vs 956-1027 A^3): in 9 of 17 chains it joins a buried region past E353/R394.
- cl4 ligands are larger (median 363 vs 250 A^3) and sit less completely inside the pocket (inside share 0.83 vs 0.98); filling stays similar because their pockets are larger too (929 vs 768 A^3).
- These are crystal-structure pocket geometries only; they say nothing about cellular activity or receptor degradation.

## Inputs (rule 1.10)

- SV00 card: results/sv/SV00/chain_table.csv (788 chains; cluster, construct_mutations = wild type 'none', ligand_chemotype), restored mmCIFs data/raw/, data/sv/lbd_extract/receptor_chains.csv and residue_mapping.csv.
- Q18 card: results/chain_ligand_of_interest.csv; SV11 card: results/sv/SV11/site_residues.csv (ligand pocket, in_site: 18 residues); SV03 card: data/sv/SV03/core_v3_transforms.csv (core_frame.py); SV12 card: results/sv/SV12/action_class.csv, class_coverage.csv, CCD records data/sv/rcsb_chemcomp/ (atom_count_heavy).

## Steps 1-4 - grid, pocket, ligand, filling

- Grid 0.4 A (voxel 0.064 A^3) over the box of the 18 ligand-pocket residues + 6 A in the core v3 frame; Bondi radii, no hydrogens; probe 1.4 A; buried = >= 10 of 14 rays hit a receptor vdW sphere within 10 A; pocket = 6-connected buried cavity component overlapping the ligand (apo: nearest to the E2 centroid).
- Atoms: receptor heavy atoms per chain median 1820; ligand instances 803 in 781 chains, 783 kept (16 chains with several Q18 instances).
- ligand_incomplete (fewer modelled heavy atoms than the CCD): 80 chains in 42 entries, out of the primary analysis; no ligand: 7 chains. Primary set: 701 chains, 331 entries, 287 ligands (cl7 393, cl4 294, minor 14).
- Primary set: pocket 855 A^3 (IQR 659-1030, range 464-1579); ligand 300 A^3; filling 33.8 % (IQR 24.8-39.2); inside share median 0.96, < 0.90 in 230 chains; ligand overlaps > 1 buried component in 269 chains; ligand vdW outside the box in 4 chains (counted, as outside the pocket).
- pocket_open: 788/788; border points median 1018 (IQR 359-1789); E2 reference chains 15-951.
- Box edge (>= 10 rays) in the primary set: 183 chains (cl7 81, cl4 92, minor 10).
- Buriedness sensitivity (primary set): Spearman rho with the >= 10-ray volume: >= 8 rays 0.71 chain / 0.75 entry; >= 12 rays 0.67 / 0.66.

## Step 5 - references (references.json)

- Ligand 100 % = 223.7 A^3 (IQR 218.9-224.8), median EST volume over 17 wild-type EST chains (7 entries).
- Apo entries: 5 entries (7 chains) with ligand_chemotype none; 3 pass (5 chains: 2B23, 4Q13, 9VTU); 7Y8F, 7YMK: PEG in the pocket (apo_blocking_residues.csv). Apo baseline = 793.2 A^3 (IQR of entry medians 724.4-805.6; entry medians 2B23 818, 4Q13 793, 9VTU 656).
- Pocket 100 % = wild-type EST chains (fallback: fewer than 5 apo entries pass): 956.2 A^3 (IQR 557.6-1011.2); E2-bound pocket (second reference) = the same 17 chains. The apo baseline is 83 % of it.
- The E2 reference pockets split: 8 chains 505-579 A^3 (<= 21 A^3 farther than 4 A from E2) and 9 chains 956-1027 A^3, where 405-495 A^3 lie farther than 4 A from E2, lined by residues 324, 326, 353, 357, 394 (among the 8 most frequent in every one of those chains; e2_reference_pocket_lining.csv). Same entry can fall in both groups (1ERE A, B, F vs C, D, E). E2 reference filling 21.0-43.0 %.
- E2 reference centroid: median of 17 EST centroids, all within 0.37 A of it.

## Step 6 - correlations (correlations.csv; Spearman, ligand = median per analysed CCD)

| subset | ligands | rho | q | entries | rho | q |
|---|---|---|---|---|---|---|
| all | 287 | 0.44 | 2.0e-14 | 331 | 0.47 | 1.4e-18 |
| cl7 | 183 | 0.41 | 4.1e-08 | 204 | 0.46 | 1.9e-11 |
| cl4 | 108 | 0.11 | 0.39 | 125 | 0.14 | 0.2 |
| class agonist | 4 | -0.40 | 0.7 | 27 | 0.09 | 0.72 |
| class SERM | 6 | 0.37 | 0.6 | 20 | 0.42 | 0.15 |

- Class gate (SV12 class_coverage.csv, >= 5 entries): agonist, SERM (entries in the SV14 primary set: agonist 27, SERM 20); SERD (4 entries), antagonist, partial agonist, other not testable.

## Step 6 - filling tests (tests.csv, BH with the correlations)

| comparison | level | medians | effect | p | q |
|---|---|---|---|---|---|
| cl7 vs cl4 | entry (202 vs 121) | 33.0 vs 32.3 % | rank-biserial 0.13 | 0.048 | 0.14 |
| cl7 vs cl4 | chain (393 vs 294) | 35.2 vs 33.0 % | 2.2 pp [-0.4, 5.7] | 0.097 | 0.19 |
| agonist vs SERM | entry (27 vs 20) | 31.7 vs 30.5 % | rank-biserial 0.01 | 0.94 | 0.94 |
| agonist vs SERM | chain (55 vs 53) | 38.4 vs 32.9 % | 5.5 pp [-12.4, 13.3] | 0.44 | 0.6 |

- By group (filling_by_group.csv, chain medians): cl7 pocket 768 A^3, ligand 250 A^3, inside share 0.98; cl4 929, 363, 0.83; agonist 560, 220, 0.98; SERM 878, 338, 0.82.

## Files

- results/sv/SV14/: pocket_per_chain.csv, references.json, correlations.csv, tests.csv, verdict.md, SUMMARY.md, env.txt, filling_by_group.csv, filling_entry_level.csv, numbers.json, e2_reference_pocket_lining.csv, apo_blocking_residues.csv
- figures/sv/SV14/: SV14_pocket_vs_ligand.png/.svg (main), SV14_filling.png/.svg; panels figures/panels/sv/SV14/ (3 PNG)
- data/sv/SV14/: atoms_core.npz, chains.csv, ligand_instances.csv, apo_check.csv, e2_centroid.json, grid_per_chain.csv, extract_stdout.txt, grid_stdout.txt
- scripts/sv/: SV14_extract.py, SV14_grid.py, SV14_analysis.py, SV14_e2_lining.py, SV14_apo_blockers.py, SV14_figures.py, SV14_summary.py; scripts/sv/common/class_style.py
