Not done / substituted: (1) Correction after the first run, disclosed: the first run paired E2 per ligand x assay x standard type (my declaration), so assay CHEMBL1909145 (reports IC50 and Ki for both compounds; type-wise dpX differ by <= 0.01) counted twice for 7 ligands; step 3 asks for one dpX per binding assay, so the reported run uses one value per ligand x assay (median over its types). This moved the entry-level cl7 vs cl4 test from q 0.033 to 0.057 (ligand level q 0.15 -> 0.15; entry-level dpX vs ligand volume q 2.9e-04 -> 0.0021); first-run values kept in affinity_per_ligand.csv dpx_median_type_units and tests_first_run_type_units.csv, not used in the verdict. (2) Class comparison not tested: coverage gate [DEFAULT] (entries with relative dpX: agonist 3, SERM 18; agonist < 5). (3) Rule 1.4 chain level not run: the step fixes ligand (primary) and entry (secondary) levels and dpX is one value per ligand. (4) Declared before the outcomes where the run file is open: ligand = CCD (SV10/SV14 unit) carrying its ChEMBL parent's CHEMBL206 activities; a pair needs the same assay id and the same standard type; fallback = median over activities, Ki/Kd if any, else IC50; hollow plot position = fallback minus the E2 median of the same type family (drawing only); spread check on primary-set ligands with >= 2 paired assays, SD ddof 1, descriptive (no test); entry = SV10 entry ligand and entry group, class = SV12 row (entry, entry ligand); class test also needs >= 5 entries in the SV15 tested set; Spearman inputs = medians over SV14 primary-analysis chains whose ligand_ids_used equals the CCD; BH over the primary family, the sensitivity run (connectivity_only included, run file SV12 step 2) is its own BH family. (5) No network call (every field was in results/sv/SV12/chembl_activities.csv); nothing installed. (6) Supporting files not named in the run file: paired_assays.csv, entry_level.csv, tests_first_run_type_units.csv, numbers.json.

# SV15 - Affinities across assay systems, relative to E2 in the same assay

**Verdict**

- Same-assay E2 reference: 19 binding assays pair a ligand with E2 (31 ligand x assay pairs, 15 ChEMBL parents); 14 primary-set ligands (E2 excluded) get a relative dpX, 48 only an absolute value; median dpX -0.74 (IQR -1.90 to 0.00), 3 of 14 bind tighter than E2; 9 rest on one assay.
- Spread check (5 ligands with >= 2 paired assays): SD across assays median 0.86 (absolute pX) vs 0.56 (dpX), pooled 1.04 vs 0.67; smaller with the E2 reference in 4 of 5.
- cl7 vs cl4 (SV10 groups): median dpX -2.17 vs -0.23 (7 vs 6 ligands, q 0.15); entries -2.18 vs -0.23 (10 vs 18, q 0.057): not significant at q < 0.05 at either level. Classes: not tested (agonist 3 entries with relative dpX, < 5).
- dpX vs SV14 filling %: rho -0.03 ligands (q 0.92), -0.02 entries (q 0.92); vs ligand volume: rho 0.49 ligands (q 0.14, n 14), 0.60 entries (q 0.0021, n 31; entries repeat ligands: OHT 7, RAL 6, GEN 3, C3D 3).
- These are in-vitro binding affinities to the ERa target (CHEMBL206) relative to E2; they say nothing about cellular activity, receptor degradation or H12 order.

## Inputs (rule 1.10)

- SV12 card: results/sv/SV12/chembl_activities.csv (1170 CHEMBL206 activities, from the cached ChEMBL responses), ligand_ids.csv (CCD -> parent, match type; E2 = EST -> CHEMBL135), action_class.csv, class_coverage.csv (gate classes: agonist, SERM).
- SV10 card: results/sv/SV10/ligand_druglikeness.csv (ligand groups, compared = cl7/cl4) and entry_druglikeness.csv (entry group, entry ligand).
- SV14 card: results/sv/SV14/pocket_per_chain.csv (filling_pct, ligand_volume_A3, in_primary_analysis).
- results/potency_by_cluster.csv is not read or changed; SV15 supersedes it inside the SV section only.

## Step 1-2 - cleaning and assay systems

- Cleaning cascade (activities kept): raw 1170; target CHEMBL206 1170; relation '=' 968; pChEMBL present 551; no validity comment 551; not a potential duplicate 409 (removed 142).
- Systems: binding (B with Ki/Kd/IC50) 273 activities (IC50 223, Ki 45, Kd 5; 73 parents, 167 assays); functional (F) 49, listed only (20 ligands, column functional_types); dropped 87 (B with other types: EC50 72, AC50 9, ED50 1; A 5).
- E2 (CHEMBL135): 47 cleaned binding activities in 46 assays (IC50 34, Ki 9, Kd 4 assay values); E2 medians for the hollow plot position: Ki/Kd 8.70, IC50 8.71.

## Step 3-4 - relative affinity to E2 and fallback (affinity_per_ligand.csv, paired_assays.csv)

- Pairs: 38 same-assay same-type pairs (IC50 30, Ki 8) -> 31 ligand x assay dpX values in 19 assays; 0 values lost to a type mismatch inside an assay; 0 pairs in mutant-variant assays.
- Per ligand (306 CCDs): no_chembl_parent 146, no_binding_data 86, absolute 58, relative 15, reference (E2) 1. Primary set: no_chembl_parent 146, no_binding_data 66, absolute 48, relative 14, reference (E2) 1; connectivity_only: no_binding_data 20, absolute 10, relative 1 (relative: PTI).
- Absolute (hollow, never tested), primary set: 48 ligands (fallback from Ki/Kd 9, from IC50 39); groups cl4 35, cl7 11, mixed cl7/cl4 1, minor 1; classes unlabelled 42, SERD 4, SERM 1, agonist 1.

| CCD | SV10 group | class | dpX median | RBA % | paired assays | filling % (SV14) | entries |
|---|---|---|---|---|---|---|---|
| C3D | cl4 | SERM | +0.85 | 708 | 1 | 34.2 | 3 |
| DES | cl7 | agonist | +0.63 | 424 | 2 | 42.0 | 2 |
| RAL | cl4 | SERM | +0.18 | 151 | 5 | 35.8 | 6 |
| IOK | cl4 | unlabelled | +0.00 | 100 | 1 | 32.9 | 1 |
| ESL | cl7 | SERM | -0.02 | 95.5 | 1 | 22.3 | 1 |
| 53Q | cl4 | SERM | -0.45 | 35.1 | 1 | 40.4 | 1 |
| OHT | cl4 | SERM | -0.71 | 19.7 | 2 | 33.4 | 7 |
| J3Z | cl7 | agonist | -0.78 | 16.7 | 4 | 22.6 | 1 |
| IOG | cl4 | unlabelled | -0.90 | 12.6 | 1 | 34.5 | 1 |
| GW5 | mixed cl7/cl4 | unlabelled | -1.11 | 7.76 | 1 | 21.5 | 2 |
| B81 | cl7 | unlabelled | -2.17 | 0.676 | 1 | 40.9 | 1 |
| GEN | cl7 | unlabelled | -2.19 | 0.646 | 8 | 37.3 | 3 |
| 5VP | cl7 | unlabelled | -2.51 | 0.309 | 1 | 45.8 | 1 |
| 2OH | cl7 | unlabelled | -3.15 | 0.0708 | 1 | 26.2 | 1 |

## Step 5 - spread check (assay_spread.csv)

| CCD | paired assays | SD absolute pX | SD dpX | SD E2 pX | ratio |
|---|---|---|---|---|---|
| DES | 2 | 0.12 | 0.54 | 0.42 | 4.50 |
| GEN | 8 | 1.00 | 0.81 | 0.45 | 0.81 |
| J3Z | 4 | 0.74 | 0.35 | 0.56 | 0.47 |
| OHT | 2 | 2.38 | 0.84 | 1.53 | 0.35 |
| RAL | 5 | 0.86 | 0.56 | 0.57 | 0.65 |

- Median ratio SD(dpX)/SD(pX) 0.65; median reduction 0.30 log units. DES (2 assays): SD dpX 0.54 > SD pX 0.12, because its E2 values vary more (SD 0.42).

## Step 6 - comparisons (tests.csv; E2 excluded; BH per family)

| set | comparison | level | n | medians / rho | p | q |
|---|---|---|---|---|---|---|
| primary | dpX by SV10 ligand group | ligand | 7 vs 6 | -2.17 vs -0.23 (r_rb -0.57) | 0.1 | 0.15 |
| primary | dpX by SV10 ligand group | entry | 10 vs 18 | -2.18 vs -0.23 (r_rb -0.54) | 0.019 | 0.057 |
| primary | dpX by action class | ligand | - | not tested: coverage gate: fewer than 2 classes with >= 5 entries in the SV15 tested set (agonist: 3 entries, 2 ligands; SERM: 18 entries, 5 ligands) | - | - |
| primary | dpX by action class | entry | - | not tested: coverage gate: fewer than 2 classes with >= 5 entries in the SV15 tested set (agonist: 3 entries, 2 ligands; SERM: 18 entries, 5 ligands) | - | - |
| primary | dpX vs filling % (SV14) | ligand | 14 | rho -0.03 | 0.92 | 0.92 |
| primary | dpX vs filling % (SV14) | entry | 31 | rho -0.02 | 0.92 | 0.92 |
| primary | dpX vs ligand volume (SV14) | ligand | 14 | rho +0.49 | 0.072 | 0.14 |
| primary | dpX vs ligand volume (SV14) | entry | 31 | rho +0.60 | 3.5e-04 | 0.0021 |
| sensitivity | dpX by SV10 ligand group | ligand | 7 vs 7 | -2.17 vs +0.00 (r_rb -0.59) | 0.073 | 0.11 |
| sensitivity | dpX by SV10 ligand group | entry | 10 vs 19 | -2.18 vs +0.00 (r_rb -0.55) | 0.017 | 0.052 |
| sensitivity | dpX by action class | ligand | - | not tested: coverage gate: fewer than 2 classes with >= 5 entries in the SV15 tested set (agonist: 3 entries, 2 ligands; SERM: 18 entries, 5 ligands) | - | - |
| sensitivity | dpX by action class | entry | - | not tested: coverage gate: fewer than 2 classes with >= 5 entries in the SV15 tested set (agonist: 3 entries, 2 ligands; SERM: 18 entries, 5 ligands) | - | - |
| sensitivity | dpX vs filling % (SV14) | ligand | 15 | rho +0.04 | 0.89 | 1 |
| sensitivity | dpX vs filling % (SV14) | entry | 32 | rho -0.00 | 1 | 1 |
| sensitivity | dpX vs ligand volume (SV14) | ligand | 15 | rho +0.52 | 0.048 | 0.096 |
| sensitivity | dpX vs ligand volume (SV14) | entry | 32 | rho +0.60 | 2.6e-04 | 0.0016 |

- Class gate: SV12 gate classes agonist, SERM; entries with relative dpX in the SV15 tested set {'agonist': 3, 'SERM': 18}, ligands {'agonist': 2, 'SERM': 5}; passing ['SERM']. Classes are drawn, not tested.
- Sensitivity (connectivity_only included, +1 ligand): cl7 vs cl4 q 0.11 ligands / 0.052 entries; filling rho +0.04 / -0.00; ligand volume rho +0.52 (q 0.10) / +0.60 (q 0.0016); q < 0.05 calls equal to the primary set in all 6 tested rows.

## Figures

- figures/sv/SV15/SV15_dpx_by_group.png/.svg: dpX per ligand by SV10 group and by SV12 class (strip + box of relative ligands; hollow = absolute, plotted at fallback minus the E2 median of the same type).
- figures/sv/SV15/SV15_dpx_vs_filling.png/.svg: dpX vs SV14 filling % (relative filled and labelled, absolute hollow).
- Panels: figures/panels/sv/SV15/ (3 PNG).

## Files

- results/sv/SV15/affinity_per_ligand.csv (80.7 KB)
- results/sv/SV15/assay_spread.csv (1.2 KB)
- results/sv/SV15/tests.csv (3.6 KB)
- results/sv/SV15/verdict.md (1.0 KB)
- results/sv/SV15/SUMMARY.md (10.2 KB)
- results/sv/SV15/env.txt (0.4 KB)
- results/sv/SV15/paired_assays.csv (6.5 KB)
- results/sv/SV15/entry_level.csv (48.0 KB)
- results/sv/SV15/tests_first_run_type_units.csv (2.2 KB)
- results/sv/SV15/numbers.json (3.7 KB)
- figures/sv/SV15/SV15_dpx_by_group.png (178.8 KB)
- figures/sv/SV15/SV15_dpx_by_group.svg (100.1 KB)
- figures/sv/SV15/SV15_dpx_vs_filling.png (131.6 KB)
- figures/sv/SV15/SV15_dpx_vs_filling.svg (58.1 KB)
- figures/panels/sv/SV15/SV15_dpx_by_group__by_class.png (95.2 KB)
- figures/panels/sv/SV15/SV15_dpx_by_group__by_cluster.png (90.6 KB)
- figures/panels/sv/SV15/SV15_dpx_vs_filling__scatter.png (131.4 KB)
- scripts/sv/SV15_affinity.py (29.9 KB)
- scripts/sv/SV15_figures.py (8.9 KB)
- scripts/sv/SV15_summary.py (17.4 KB)
- docs/EXPLORATION.md: appended 'Suggested, not done (SV15, 2026-11-08)' (rule 1.1)
