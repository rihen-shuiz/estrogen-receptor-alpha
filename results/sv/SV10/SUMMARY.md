Not done / substituted: nothing substituted (RDKit and statsmodels present; 306 of 306 SMILES parsed, 0 empty rows). (1) results/potency_by_cluster.csv is not an `out` of a card: it is evidence of the X16 DONE line in tracker/esr1.log and of tracker/prov/X15.json and X16.json (cards X15/X16 list figures/potency_by_cluster.png and results/potency_tests.csv). (2) Declared before running where the run file is open: PAINS = one count over PAINS_A + PAINS_B + PAINS_C (per-family counts in ligand_druglikeness.csv); QED ALERTS is treated as an alert count (Fisher); the size-check OLS ranks the count values themselves; entry group = the same >= 80 % rule on the entry's kept chains; entry ligand = the ligand_id_primary of most of its chains (9 ties resolved by the first chain); the md table shows rho at ligand level (entry rho in the csv); a no-variation rule was declared and not triggered. (3) 5 QED properties equal a Lipinski/Veber row for every ligand (qed_mw = mw, qed_alogp = clogp, qed_hbd = hbd, qed_psa = tpsa, qed_rotb = rotb); kept as written, so each BH family has 21 rows. (4) Supporting files not named in the run file: entry_druglikeness.csv, numbers.json.

# SV10 - Ligand drug-likeness, cl7 vs cl4 with a size check

**Verdict ([DEFAULT] size check): 10 of 21 parameters add beyond size (7 distinct quantities: hbd, hba, rotb, tpsa, qed_alerts, brenk_n, pains_n); 11 restate size.** These parameters describe oral drug-likeness of the crystallised ligands, not binding, conformation or activity.

## Inputs (rule 1.10)

- Q08 card out results/ligands.csv: canonical_isomeric_smiles for all 306 ligands (smiles fallback not needed).
- Y09 card out results/ext/ligand_descriptors.csv: ligand set, structural_majority (>= 80 % rule), heavy_atoms (RDKit heavy atoms of the SMILES equal it for 306 of 306). Groups recomputed from the SV00 chain table agree for 306 of 306 ligands.
- results/potency_by_cluster.csv (X15/X16 evidence): pchembl for 82 ligands, one value each.
- SV00 card out results/sv/SV00/chain_table.csv: entry_id, s08_cluster and ligand_id_primary (first row of Q18 results/chain_ligand_of_interest.csv per chain) for the entry level.

## Key numbers

- Ligands: 306 described; groups cl7 192, cl4 104, mixed cl7/cl4 6, none 2, minor 2. Compared: cl7 192, cl4 104; with pChEMBL cl7 27, cl4 52 (of 82).
- Entries: 351; groups cl7 219, cl4 124, minor 4, mixed cl7/cl4 4; 2 without a ligand of interest; compared cl7 217, cl4 124 (290 distinct ligands; 11 entries hold >1 ligand of interest; in 8 the entry group differs from its ligand's group).
- Size: heavy atoms median [IQR] cl7 24 [21-30] vs cl4 34 [32-37] (ligands, Mann-Whitney p 1.1e-21); entries 24 vs 34 (p 1.5e-25).
- Unadjusted: q < 0.05 for 16 of 21 parameters at ligand level, 17 at entry level. Size-adjusted q < 0.05: 12 ligand, 10 entry.
- Alerts (all 306 ligands): Brenk >= 1 in 156, PAINS >= 1 in 10 (indol_3yl_alk(461) x8, anil_di_alk_C(246) x1, anil_di_alk_D(198) x1).
- LE and LipE use pChEMBL, which mixes IC50, Ki and Kd from different assays.

| parameter | cl7 (ligands) | cl4 (ligands) | q ligand | q entry | rho HA | size-adj q lig / entry | direction at equal HA | label |
|---|---|---|---|---|---|---|---|---|
| mw | 344 [296-411] | 470 [436-512] | 1.5e-17 | 7.9e-22 | +0.97 | 0.027 / 0.11 |  | restates size |
| clogp | 4.58 [3.61-5.5] | 5.58 [4.87-5.95] | 4.4e-09 | 1.1e-14 | +0.59 | 0.69 / 0.14 |  | restates size |
| hbd | 2 [2-2] | 2 [1-2] | 3.4e-06 | 9e-10 | -0.14 | 0.00011 / 8.3e-08 | cl4 lower | adds beyond size |
| hba | 3 [2-5] | 4 [3-6] | 9.7e-05 | 5.5e-06 | +0.63 | 0.0025 / 0.0027 | cl4 lower | adds beyond size |
| lipinski_violations | 85/192 (44 %) | 80/104 (77 %) | 1.2e-07 | 1.4e-12 | +0.69 | 0.25 / 0.66 |  | restates size |
| rotb | 3 [2-5] | 6 [5-8] | 1.3e-21 | 5.7e-27 | +0.88 | 0.0013 / 1.1e-05 | cl4 higher | adds beyond size |
| tpsa | 58.2 [40.5-87] | 60.8 [53-71.4] | 0.27 | 0.41 | +0.47 | 2.1e-05 / 2.6e-07 | cl4 lower | adds beyond size |
| veber_violations | 7/192 (4 %) | 11/104 (11 %) | 0.029 | 0.03 | +0.38 | 0.14 / 0.11 |  | restates size |
| qed | 0.638 [0.489-0.757] | 0.473 [0.337-0.539] | 1.4e-09 | 1.2e-13 | -0.77 | 0.066 / 0.33 |  | restates size |
| qed_mw (= mw) | 344 [296-411] | 470 [436-512] | 1.5e-17 | 7.9e-22 | +0.97 | 0.027 / 0.11 |  | restates size |
| qed_alogp (= clogp) | 4.58 [3.61-5.5] | 5.58 [4.87-5.95] | 4.4e-09 | 1.1e-14 | +0.59 | 0.69 / 0.14 |  | restates size |
| qed_hba | 3 [2-5] | 4 [4-6] | 1.1e-07 | 2.8e-09 | +0.69 | 0.082 / 0.078 |  | restates size |
| qed_hbd (= hbd) | 2 [2-2] | 2 [1-2] | 3.4e-06 | 9e-10 | -0.14 | 0.00011 / 8.3e-08 | cl4 lower | adds beyond size |
| qed_psa (= tpsa) | 58.2 [40.5-87] | 60.8 [53-71.4] | 0.27 | 0.41 | +0.47 | 2.1e-05 / 2.6e-07 | cl4 lower | adds beyond size |
| qed_rotb (= rotb) | 3 [2-5] | 6 [5-8] | 1.3e-21 | 5.7e-27 | +0.88 | 0.0013 / 1.1e-05 | cl4 higher | adds beyond size |
| qed_arom | 3 [2-3] | 3 [3-3] | 8.1e-07 | 2.6e-11 | +0.54 | 0.86 / 0.4 |  | restates size |
| qed_alerts | 95/192 (49 %) | 40/104 (38 %) | 0.11 | 0.91 | +0.12 | 0.00028 / 0.0035 | cl4 lower | adds beyond size |
| brenk_n | 95/192 (49 %) | 53/104 (51 %) | 0.9 | 0.36 | +0.37 | 5e-06 / 1.7e-05 | cl4 lower | adds beyond size |
| pains_n | 0/192 (0 %) | 10/104 (10 %) | 3.5e-05 | 0.00012 | +0.14 | 0.00073 / 0.0029 | cl4 higher | adds beyond size |
| le | 0.441 [0.394-0.524] | 0.358 [0.329-0.39] | 0.0006 | 1.2e-10 | -0.79 | 0.75 / 0.78 |  | restates size |
| lipe | 3.04 [2.33-3.92] | 2.84 [2.19-3.68] | 0.27 | 4e-07 | -0.33 | 0.69 / 0.29 |  | restates size |

Full table with definitions, entry-level n and medians, p values and test statistics: druglikeness_by_group.csv / druglikeness_by_group.md.

## verdict.md

> # SV10 verdict ([DEFAULT] size check): 10 of 21 drug-likeness parameters add beyond size (size-adjusted q < 0.05 at ligand and entry level; 7 distinct quantities, since 3 are QED copies); 11 restate size.
> - Adds beyond size, direction at equal heavy atoms (ligand-level values; q ligand / entry): hbd = qed_hbd: cl4 lower (cl7 2 [2-2] vs cl4 2 [1-2]; size-adj q 0.00011 / 8.3e-08); hba: cl4 lower (cl7 3 [2-5] vs cl4 4 [3-6]; size-adj q 0.0025 / 0.0027); rotb = qed_rotb: cl4 higher (cl7 3 [2-5] vs cl4 6 [5-8]; size-adj q 0.0013 / 1.1e-05); tpsa = qed_psa: cl4 lower (cl7 58.2 [40.5-87] vs cl4 60.8 [53-71.4]; size-adj q 2.1e-05 / 2.6e-07); qed_alerts: cl4 lower (cl7 95/192 (49 %) vs cl4 40/104 (38 %); size-adj q 0.00028 / 0.0035); brenk_n: cl4 lower (cl7 95/192 (49 %) vs cl4 53/104 (51 %); size-adj q 5e-06 / 1.7e-05); pains_n: cl4 higher (cl7 0/192 (0 %) vs cl4 10/104 (10 %); size-adj q 0.00073 / 0.0029).
> - Restates size: mw, clogp, lipinski_violations, veber_violations, qed, qed_mw, qed_alogp, qed_hba, qed_arom, le, lipe.
> - Size context: heavy atoms median cl7 24 vs cl4 34 (ligands, Mann-Whitney p 1.1e-21).
> - These parameters describe oral drug-likeness of the crystallised ligands (hormone, fragments, tool compounds included), not binding, conformation or activity; LE/LipE use pChEMBL that mixes assay types.

## Files

- results/sv/SV10/ligand_druglikeness.csv (113.2 kB)
- results/sv/SV10/entry_druglikeness.csv (77.6 kB)
- results/sv/SV10/druglikeness_by_group.csv (13.4 kB)
- results/sv/SV10/druglikeness_by_group.md (5.9 kB)
- results/sv/SV10/verdict.md (1.3 kB)
- results/sv/SV10/numbers.json (2.7 kB)
- results/sv/SV10/env.txt (0.3 kB)
- scripts/sv/SV10_druglikeness.py (25.8 kB)
- scripts/sv/SV10_summary.py (7.6 kB)
- results/sv/SV10/SUMMARY.md (this file)
- No figure (run file: 'No figure').

Run order: .venv/bin/python -B -I scripts/sv/SV10_druglikeness.py; scripts/sv/SV10_summary.py
