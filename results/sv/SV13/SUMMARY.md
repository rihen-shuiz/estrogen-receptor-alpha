Not done / substituted: (1) PyMOL not available (PRE addition; env.txt) -> Part C uses the step 12 fallback: C-alpha traces projected with matplotlib onto the saved view; the representative's cartoon is a thick C-alpha trace, its ligand of interest sticks with bonds from interatomic distances. (2) Part B: every comparison fails the coverage gate (0 of 22 rows testable), so no model or rank correlation was fitted and the forest plot shows the gate reason in place of estimates. (3) Representatives follow the [DEFAULT] as written: the SERM representative 2BJ4|A|1 has no modelled H12 (only 528 of 528-550 observed) and the SERD representative 7MSA|C|1 7 of 10 H12 residues, so their B z colouring is absent / partial (the existing RMSD matrix compares chains on their common positions). (4) Declared before any class comparison, where the run file is open: sensitivity set = all rows incl. connectivity_only (69 chains, all unlabelled, so counts are identical); C2b uses the ligands SV15 itself tests (relative same-assay dpX; E2 reference and absolute-only ligands excluded, SV15 steps 4 and 6) with a gate of >= 5 agonist entries with dpX; psi 537 outcome = wrapped difference from the circular mean of 435 cl7 chains (156.5 deg); representative ties by resolution, then chain_key; view = principal axes of the 5TM1|A|1 selection (PyMOL 'orient' convention), signs with its H12 centroid at +x/+y. (5) Shared code: optional edge_map parameter added to scripts/sv/common/h12_map.py (default output unchanged). (6) The model code was exercised once on seeded random pseudo-labels (data/sv/SV13/model_code_check.json: statuses only, no estimate used) so that the re-run with the teammate labels has a checked path. (7) Supporting files not named in the run file: chain_data.csv, gate.csv, numbers_tests.json, numbers.json (results/sv/SV13/), data/sv/SV13/*.

# SV13 - H12 by action class (provisional SV12 classes; H12 window 538-547)

**Verdict** (results/sv/SV13/verdict.md)

- C1 SERD vs SERM (H12 complete, fraction observed, mobility z; all chains and within cl4): **not testable** - SERD 4 entries (all in cl4) vs SERM 20 (19 in cl4); coverage gate needs >= 5 entries per side.
- C1 antagonist vs SERM: **not testable** (0 antagonist entries in SV12). C2a partial agonist vs agonist within cl7: **not testable** (0 partial-agonist entries; 25 agonist entries in cl7).
- C2b agonists vs dpX within cl7 (the VDR superagonist point): **not testable** (3 agonist entries with an SV15-tested dpX: DES 2, J3Z 1 entries; E2 is the reference and 4OH has an absolute value only).
- So no class-H12 link can be stated: 0 of 22 comparison x stratum x outcome rows pass the gate (sensitivity set with connectivity_only rows: same counts, 69 added chains carry no class).
- Drawn, not tested (panels, maps): chains with complete H12 538-547 agonist 53 of 55, SERM 19 of 53, SERD 5 of 9; classes follow the clusters (agonist 25 of 27 entries in cl7, SERM 19 of 20 and SERD 4 of 4 in cl4).
- Crystal structures show H12 order and position (H12 disorder, H12 mobility), not receptor degradation or cellular activity; the classes are provisional ChEMBL-mechanism labels (SV12).

## Inputs (rule 1.10)

- SV12 card: results/sv/SV12/chain_class.csv (class per chain, primary_set), class_coverage.csv (cross-check: 56 of 56 class x cluster cells equal).
- SV00 chain_table.csv; SV01a chain_outcomes.csv (complete, fraction observed); SV01b mobility_per_chain.csv, mobility_per_residue.csv, lattice_contacts_per_chain.csv, model_data.csv (rebuilt Q/X/lattice/mutation/peptide predictors equal model_data.csv: max |diff| cl7 0.0, cl4 0.0); SV02 phipsi_per_chain.csv; SV03 geometry_per_chain.csv, core_v3_transforms.csv; SV15 affinity_per_ligand.csv.
- Y03 card: results/ext/embedding_rmsd.csv (mds1/mds2); S04/S05 cards: results/rmsd_entry.npy + rmsd_entry_index.csv; S12: results/medoids.csv (cl7 medoid 5TM1|A|1); Q18: results/chain_ligand_of_interest.csv; restored mmCIFs data/raw/.

## Part A - classes on the existing maps (class_by_cluster.csv; primary set)

| class | chains / entries all | cl7 | cl4 | minor | complete H12 chains |
|---|---|---|---|---|---|
| agonist | 55 / 27 | 52 / 25 | 1 / 1 | 2 / 1 | 53 |
| partial_agonist | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 |
| SERM | 53 / 20 | 4 / 3 | 48 / 19 | 1 / 1 | 19 |
| SERD | 9 / 4 | 0 / 0 | 9 / 4 | 0 / 0 | 5 |
| antagonist | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 |
| other | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 |
| unlabelled | 602 / 271 | 350 / 174 | 242 / 94 | 10 / 4 | 427 |

- class_by_cluster.csv also lists every S08 cluster id (1-7) and the 'all' set; connectivity_only chains (69) are all unlabelled.
- Maps redrawn with plot_h12_map (SV03 function): fill = S08 cluster, marker shape and edge = SV12 class, unlabelled small grey-edged circles. MDS map: 788 chains; H12 map: 748 chains with an H12 axis.

## Part B - tests fixed in advance (tests.csv, gate.csv)

| comparison | stratum | side a (chains / entries) | side b (chains / entries) | gate | verdict |
|---|---|---|---|---|---|
| C1_SERD_vs_SERM | all | SERD 9 / 4 | SERM 53 / 20 | fail: SERD 4 entries < 5 | not testable |
| C1_SERD_vs_SERM | cl4 | SERD 9 / 4 | SERM 48 / 19 | fail: SERD 4 entries < 5 | not testable |
| C1_antagonist_vs_SERM | all | antagonist 0 / 0 | SERM 53 / 20 | fail: antagonist 0 entries < 5 | not testable |
| C1_antagonist_vs_SERM | cl4 | antagonist 0 / 0 | SERM 48 / 19 | fail: antagonist 0 entries < 5 | not testable |
| C2a_partial_agonist_vs_agonist | cl7 | partial_agonist 0 / 0 | agonist 52 / 25 | fail: partial_agonist 0 entries < 5 | not testable |
| C2b_agonist_rho_dpX | cl7 | agonist with dpX 6 / 3 | - | fail: agonist entries with a tested SV15 dpX: 3 < 5 (DES 2, J3Z 1; 25 agonist entries in cl7) | not testable |

- tests.csv: 88 rows (2 sets x 22 comparison x stratum x outcome rows x 2 levels), all 'not testable (coverage gate)'; models and C2b code are in scripts/sv/SV13_analysis.py (class + Q + X with lattice terms + mutations + peptide; chain level entry-clustered, entry level one chain per entry; C2b Spearman, chain level bootstrap over entries; BH within each comparison).
- Sensitivity set (connectivity_only included): identical gate counts; verdicts ['not testable'].
- Code-path check on seeded random pseudo-labels: 22 of 22 model / rho fits ran (statuses only, no estimate kept).

## Part C - pictures of H12 and the ligand (representatives.csv, view.txt)

| class | representative | ligand | cluster | median RMSD to class (A) | resolution (A) | H12 538-547 modelled | runner-up |
|---|---|---|---|---|---|---|---|
| agonist | 5WGD|B|1 | EST | cl7 | 0.382 | 1.80 | 10 of 10 | 5DXB|B|1 (0.426) |
| SERM | 2BJ4|A|1 | OHT | cl4 | 0.721 | 2.00 | 0 of 10 | 2JFA|A|1 (0.789) |
| SERD | 7MSA|C|1 | ZNM | cl4 | 0.560 | 2.24 | 7 of 10 | 6SBO|A|1 (0.705) |

- Ligand atoms (restored mmCIF, SV14 instance rule) vs data/sv/SV14/atoms_core.npz: 5WGD|B|1 20 atoms, max 1.9e-06 A; 2BJ4|A|1 29 atoms, max 1.9e-06 A; 7MSA|C|1 37 atoms, max 1.7e-06 A.
- H12 B z colour range [0.0, 3.0] (representatives' H12 B z 0.75 to 2.92; SV01b per-residue z vs the chain's core v3 C-alpha B).
- Panels: 6 class panels (partial agonist, antagonist and other say 0 chains) + 1 overlay of the 3 representatives; captions give chains, entries and chains with complete H12.
- view.txt:
  - frame: core v3 (data/sv/SV03/core_v3_transforms.csv; x_core = (x_deposited - m) @ R_fit + c, scripts/sv/common/core_frame.py)
  - selection: 5TM1|A|1 (cl7 medoid) C-alpha at the core v3 positions [[338, 362], [371, 396], [420, 439], [441, 456], [496, 531]] (122 atoms) + H12 538-547 (10 atoms); 132 atoms
  - method: principal axes of the selection (PyMOL 'orient' convention): screen x = PC1, screen y = PC2, screen z = x cross y (towards the viewer); signs: the selection's H12 centroid at +x and +y
  - projection: orthographic, [sx, sy, depth] = ROT @ (x_core - CENTRE); drawn with matplotlib (PyMOL not available: step 12 fallback)
  - CENTRE (A, core v3 frame): 8.9072 -6.2084 28.2263
  - ROT row 1 (screen x): 0.783534 -0.282425 0.553453
  - ROT row 2 (screen y): -0.345023 -0.938546 0.009519
  - ROT row 3 (screen z): 0.516753 -0.198412 -0.832826
  - limits (A, screen): x -32.16 47.55; y -22.82 32.59 (box of everything drawn in any panel + 2 A); equal aspect
  - H12 centroid of 5TM1|A|1 on screen (A): 13.98 6.24 2.97
  - used unchanged for every panel of figures/sv/SV13/SV13_h12_class_panels.png

## Figures

- figures/sv/SV13/SV13_h12_class_panels.png/.svg (main): per-class H12 panels + overlay, one saved view.
- figures/sv/SV13/SV13_forest.png/.svg (main): C1/C2a class effects and C2b rho; every row shows its gate reason.
- figures/sv/SV13/SV13_mds_map.png/.svg and SV13_h12_map.png/.svg (backup): the two redrawn maps.

## Files

- results/sv/SV13/class_by_cluster.csv (5.0 KB)
- results/sv/SV13/tests.csv (19.0 KB)
- results/sv/SV13/representatives.csv (1.0 KB)
- results/sv/SV13/view.txt (1.1 KB)
- results/sv/SV13/verdict.md (1.2 KB)
- results/sv/SV13/chain_data.csv (144.6 KB)
- results/sv/SV13/gate.csv (1.4 KB)
- results/sv/SV13/numbers_tests.json (6.2 KB)
- results/sv/SV13/numbers.json (2.6 KB)
- results/sv/SV13/env.txt (0.5 KB)
- data/sv/SV13/model_code_check.json (5.3 KB)
- data/sv/SV13/panel_ligand_atoms.csv (11.8 KB)
- data/sv/SV13/panel_ligand_bonds.csv (3.6 KB)
- data/sv/SV13/panel_meta.json (1.6 KB)
- data/sv/SV13/panel_traces.csv (518.6 KB)
- figures/sv/SV13/SV13_forest.png (435.2 KB)
- figures/sv/SV13/SV13_forest.svg (51.5 KB)
- figures/sv/SV13/SV13_h12_class_panels.png (801.5 KB)
- figures/sv/SV13/SV13_h12_class_panels.svg (215.1 KB)
- figures/sv/SV13/SV13_h12_map.png (139.5 KB)
- figures/sv/SV13/SV13_h12_map.svg (147.1 KB)
- figures/sv/SV13/SV13_mds_map.png (148.4 KB)
- figures/sv/SV13/SV13_mds_map.svg (154.9 KB)
- figures/panels/sv/SV13/SV13_forest__C1_logit.png (117.2 KB)
- figures/panels/sv/SV13/SV13_forest__C1_mobility_C2a.png (146.9 KB)
- figures/panels/sv/SV13/SV13_forest__C2b_rho.png (118.1 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__SERD.png (137.4 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__SERM.png (168.8 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__agonist.png (146.9 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__antagonist.png (18.5 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__other.png (18.1 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__overlay.png (175.3 KB)
- figures/panels/sv/SV13/SV13_h12_class_panels__partial_agonist.png (19.5 KB)
- figures/panels/sv/SV13/SV13_h12_map__h12_map.png (139.0 KB)
- figures/panels/sv/SV13/SV13_mds_map__mds_map.png (148.5 KB)
- scripts/sv/SV13_analysis.py (27.2 KB)
- scripts/sv/SV13_panels.py (14.6 KB)
- scripts/sv/SV13_figures.py (20.1 KB)
- scripts/sv/SV13_summary.py (15.4 KB)
- results/sv/SV13/SUMMARY.md (this file)
- scripts/sv/common/h12_map.py: optional edge_map parameter (SV13)
- docs/EXPLORATION.md: appended 'Suggested, not done (SV13, 2026-11-08)' (rule 1.1)
