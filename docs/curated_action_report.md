# Curated Action Class Analysis — ESR1 Ligand-Binding Domain

- **Base commit (`final`)**: `f66ae14b88f7026db4cfae72e0b6bd11148c7aeb`
- **Curated dataset**: `esr1_dataset_final.csv` (377 total rows, 368 real PDB x CCD ligand pairs)
- **Reproduction command**: `bash reproduce_curated.sh --dataset esr1_dataset_final.csv`
- **Outputs**:
  - Tables, audit logs, statistical summaries, and manifest: `results/12_curated_action/tables/`
  - 20 analytical figures (PNG + SVG + source CSV + captions): `results/12_curated_action/figures/`
  - 10 PyMOL 3D visualizations of H11, H12, and ligand + 2x5 contact sheet: `results/12_curated_action/pymol_ligands/`

## 1. Dataset Join and Coverage Audit

| Metric | Count |
|---|---|
| Real ligand pairs (`is_real_ligand=True`) | 368 (104 agonist, 79 partial_agonist, 119 antagonist, 64 unknown, 2 mixed) |
| Core 3 functional classes in CSV | 302 |
| Exact primary-CCD match to `chain_table.csv` (788 chains, 351 PDB) | 781 chains, 360 pairs |
| Strict `ligand_status == "single"` match | 765 chains, 355 pairs, 346 PDBs |
| Core 3 functional classes (strict single) | 631 chains, 291 pairs, 287 PDBs, 240 unique InChIKeys |
| `drug_classes` populated among real pairs | 54 / 368 (including 9 SERDs with `action_class=unknown`) |
| Structural vs CSV mutant flag mismatches | 47 pairs (e.g. `1A52` carries `Q502E` structurally while CSV marks `WT`) |
| Peptide flag mismatches (`peptide_bound` vs `has_peptide`) | 0 |
| InChIKey label conflicts | 1 (`EST` in `1QKT` labeled `partial_agonist` vs `agonist` elsewhere) |
| WT-assay on mutant-structure flagged rows | 19 pair x group rows |

## 2. Main Findings

1. **Conformational redundancy remains high across all functional classes.** Across all 351 PDB entries at the 1.5 A cut, effective PDB-weighted diversity is Hill `q1 = 2.11` (Shannon entropy 1.08 bits). Two complete-H12 representative chains cover 90% of all 788 chains within 1.0 A (48 representatives at 0.5 A). Re-clustering after strict quality filtering preserves the H12 partition (`ARI = 0.998` for known R-free <= 0.30; `ARI = 0.982` for all three QC criteria, Tier A).
2. **Manual curation breaks the near-1:1 class-to-conformation lock (revising historical claim L7).** On the 291 strict core-3 pairs, agonists sit almost exclusively in `cl7` (95/98), whereas `partial_agonist` splits across both conformations (65 in `cl7`, 13 in `cl4`), and 22 of 115 antagonists (19%) crystallize in the closed `cl7` state (predominantly in `Y537S` + coactivator peptide constructs).
3. **Within comparable strata, functional class does not shift H12 beyond ordinary crystal-copy variation (confirming L6 on a ~6x larger dataset).** In the only well-overlapped 3-class stratum (`cl7`, exact `Y537S` + peptide: 75 agonist, 63 partial_agonist, 20 antagonist), median `dist_H12c_H11axis_A` shifts are `+0.21 A` (antagonist - agonist) and `-0.17 A` (partial_agonist - agonist) versus the `+/-0.91 A` ordinary variation threshold. In `cl4` (`L536S` without peptide: 60 antagonist vs 9 partial_agonist), the shift is `-0.02 A` versus `+/-1.08 A` (`+0.43 A` in exact `L372S;L536S`, 18 vs 5). Results hold across all 11 sensitivity filters, leave-one-DOI-out checks, and direct within-DOI comparisons.
4. **Chemical scaffolds track function and H12 state, yet close analogues switch class or conformation with experimental context.** ECFP4 chemical clusters correlate strongly with functional class (`Cramer's V = 0.62` on pairs, `0.56` on unique ligands, `0.51` deduplicated by DOI) and H12 group (`V = 0.46`). However, close pairs (`Tc >= 0.55`) including stereoisomers (`AIH`/`AIU`, `86V`/`86Y`) and close analogues (`73I`/`OB9`/`OBT`, `Tc = 0.90-0.92`) switch functional class or H12 conformation across constructs.

## 3. Artifacts Index

- Full figure manifest: `results/12_curated_action/tables/figure_manifest.csv` (20 figures: `F01`–`F20`)
- Statistical effects and sensitivity tables: `results/12_curated_action/tables/effects_summary.csv`, `lodo_within_doi.csv`, `robustness_sensitivity.csv`
- 10 PyMOL H11/H12 ligand renders: `results/12_curated_action/pymol_ligands/ligands_10_selection.csv` and `00_contact_sheet_10_ligands.png`
