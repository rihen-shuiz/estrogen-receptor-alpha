Not done / substituted: (1) step 8 (medoid H12 clash test) not run: the step runs it only if the verdict is 'mixed' and the re-issued verdict is 'artifact' (scripts/sv/SV01b_clash.py is written and gated, clash_test_not_run.md; it was run once outside results/ as a smoke test on the 235 partial/not-modelled chains and that output was discarded, not used); (2) choices the run file leaves open, fixed before fitting (docstrings of SV01b_models.py and common/lattice_contacts.py): residue B = C-alpha B (altloc: highest occupancy, then blank/A); mobility = mean z over every observed residue of 528-547; core SD with ddof 0 (existing regional-B code); core v3 used unchanged, so 528-531 are in both the core and the window; lattice terms = number of heavy atoms of 528-550 within 4.0 A of an ASU copy / a symmetry copy of another chain (polymer chains; waters and ligands removed as in SV00), with the chain's own bound peptide copy left out because it is 'peptide' in block B; (3) re-issued step verdict = 'artifact'/'biology' only if the [DEFAULT] condition holds for both outcomes (complete, mobility) at both levels, else 'mixed' (the run file names one verdict for two outcomes); (4) mobility models: chain level OLS with entry-clustered sandwich errors, entry level OLS on the first matrix chain of the entry among chains with the outcome, block tests cluster-robust Wald chi2 / nested F, Wald 95 % CI b +- 1.96 SE; (5) inherited unchanged from SV01a: median imputation + 'not reported' indicators, entry mean B as 'mean B', the drop rule extended to collinearity and joint separation, McFadden gains, one BH family.

# SV01b - H12 completeness and mobility with coordinates (lattice contacts)

**Re-issued verdict ([DEFAULT] rule, lattice terms in block X): artifact** (SV01a: mixed). Complete H12 (538-547): artifact (chain artifact, entry artifact); mobility (B z-score 528-547): artifact (chain artifact, entry artifact). Sensitivity: fraction observed artifact, mobility without TLS-refined entries artifact. Step 8 not run (verdict not mixed). The change from SV01a comes from the lattice terms, which for the completeness outcome are counted on modelled atoms (see Reading).

## Inputs

- SV00 (DONE): results/sv/SV00/chain_table.csv (788 chains, tls_refined, auth_asym_id), observed_mask.npy, data/sv/backbone/backbone.npz (C-alpha B), restored mmCIFs data/raw/, regenerated Q04 tables data/sv/lbd_extract/ (residue mapping, chain identities).
- SV01a (DONE): scripts/sv/SV01a_completeness.py imported (load, predictors, imputation, drop rule, fit); refitting the SV01a predictor set with the SV01b code reproduces results/sv/SV01a: 64/64 block rows, max |d gain| 9.0e-17, max |d p| 1.1e-16; 168 odds rows, status equal 168, max |d log OR| 4.4e-16 (reproduce_SV01a.json).
- Core v3: esr_analysis/resources/config/alignment_core_v3.json core_ranges (the S03 core) 338-362, 371-396, 420-439, 441-456, 496-531.
- C-alpha B in backbone.npz equals residue_mapping B_factor_A2 for 181829 observed C-alpha (max |d| 0.0).

## 6. Lattice contacts (scripts/sv/common/lattice_contacts.py, lattice_contacts_per_residue.csv / _per_chain.csv)

- 788 chains in 351 entries; crystal frame {'standard': 350, 'standardized': 1}; own bound peptide copy agrees with SV00 chain_peptide.csv for 788/788 chains ({'asu': 417, 'none': 369, 'symmetry': 2}).
- Residues 528-550 with atoms: 14340; with C-alpha observed (SV00 mask): 14336; 4 residues have atoms but no observed C-alpha (1UOM|A|1:547, 3ERD|A|1:550, 3ERD|B|1:550, 5UFX|A|1:530; lattice contact atoms on them: 0).
- Chains with any ASU-partner contact on 528-550: 359; with any symmetry-mate contact: 614.

By cluster and H12 status (lattice_contacts_by_group.csv; mean heavy atoms of 528-550 in contact, and where):

| group | H12 | chains | entries | any ASU | any sym | mean ASU atoms | mean sym atoms | sym 528-537 | sym 538-545 | sym 546-550 | modelled atoms 546-550 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cl4 | complete | 120 | 84 | 15 | 115 | 1.75 | 24.29 | 3.77 | 8.22 | 12.30 | 30.5 |
| cl4 | not_modelled | 14 | 10 | 0 | 3 | 0.00 | 0.36 | 0.36 | 0.00 | 0.00 | 0.0 |
| cl4 | partial | 195 | 84 | 47 | 156 | 1.80 | 9.11 | 2.36 | 6.66 | 0.10 | 2.9 |
| cl7 | complete | 431 | 216 | 293 | 327 | 1.15 | 5.86 | 3.46 | 2.35 | 0.06 | 24.2 |
| cl7 | not_modelled | 7 | 6 | 0 | 0 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.0 |
| cl7 | partial | 7 | 6 | 0 | 2 | 0.00 | 1.57 | 0.57 | 1.00 | 0.00 | 2.1 |
| cl4 | all | 329 | 128 | 62 | 274 | 1.71 | 14.28 | 2.79 | 6.95 | 4.54 | 12.8 |
| cl7 | all | 445 | 221 | 293 | 329 | 1.12 | 5.70 | 3.36 | 2.29 | 0.05 | 23.5 |

## 5. Local mobility (mobility_per_chain.csv, mobility_per_residue.csv)

- mobility_z defined for cl7 442/445 and cl4 325/329 chains (undefined = no residue of 528-547 observed: 5T92|A|1, 5T97|B|1, 5TM9|A|1, 6B0F|B|1, 7RS3|B|1, 7RS3|C|1, 7YMK|C|1); core with < 2 valid B: 0.
- TLS-refined (SV00 tls_refined): cl7 348, cl4 255 chains; the TLS sensitivity keeps 95 cl7 and 72 cl4 chains.

| group | H12 | chains with mobility | median mobility_z | median loop 528-537 z | median H12 538-547 z |
|---|---|---|---|---|---|
| cl4 | complete | 120 | 1.83 | 2.56 | 1.41 |
| cl4 | not_modelled | 10 | 3.82 | 3.82 | - |
| cl4 | partial | 195 | 1.84 | 2.16 | 1.49 |
| cl7 | complete | 431 | 1.86 | 2.65 | 1.14 |
| cl7 | not_modelled | 4 | 3.55 | 3.55 | - |
| cl7 | partial | 7 | 2.16 | 2.55 | 1.83 |

## 7. Refit with lattice terms (model_blocks.csv, odds_ratios.csv, verdict_table.csv)

Gain entered first / entered last (McFadden R2 for complete and fraction observed, OLS R2 for mobility); test and BH q of the block entered last:

| outcome | cluster | level | n | Q | X | B | Q+X | test Q+X (last) | q Q+X | test B (last) | q B |
|---|---|---|---|---|---|---|---|---|---|---|---|
| complete | cl7 | chain | 445 | 0.044 / 0.047 | 0.379 / 0.259 | 0.225 / 0.095 | 0.414 / 0.283 | cluster-robust Wald 44.1 (df 12), p 1.4e-05 | 0.00021 | cluster-robust Wald 11.3 (df 5), p 0.045 | 0.14 |
| fraction_observed | cl7 | chain | 445 | 0.033 / 0.045 | 0.532 / 0.327 | 0.352 / 0.149 | 0.558 / 0.356 | cluster-robust Wald 34.0 (df 12), p 0.00068 | 0.0076 | cluster-robust Wald 15.1 (df 5), p 0.01 | 0.045 |
| mobility | cl7 | chain | 442 | 0.010 / 0.025 | 0.161 / 0.150 | 0.076 / 0.054 | 0.184 / 0.162 | cluster-robust Wald 96.6 (df 14), p 2.2e-14 | 9.9e-13 | cluster-robust Wald 17.3 (df 6), p 0.0083 | 0.038 |
| mobility_no_tls | cl7 | chain | 95 | 0.050 / 0.028 | 0.408 / 0.333 | 0.144 / 0.109 | 0.419 / 0.384 | cluster-robust Wald 164.9 (df 13), p 2e-28 | 3.7e-26 | cluster-robust Wald 65.7 (df 4), p 1.8e-13 | 6.7e-12 |
| complete | cl7 | entry | 221 | 0.261 / 0.210 | 0.440 / 0.251 | 0.315 / 0.089 | 0.687 / 0.461 | LR 22.0 (df 12), p 0.038 | 0.12 | LR 4.2 (df 5), p 0.52 | 0.66 |
| fraction_observed | cl7 | entry | 221 | 0.259 / 0.105 | 0.671 / 0.454 | 0.156 / 0.100 | 0.746 / 0.690 | LR 28.2 (df 11), p 0.0031 | 0.021 | LR 4.1 (df 2), p 0.13 | 0.27 |
| mobility | cl7 | entry | 220 | 0.035 / 0.053 | 0.096 / 0.111 | 0.098 / 0.086 | 0.157 / 0.146 | F 2.7 (df 14), p 0.00098 | 0.0097 | F 3.8 (df 6), p 0.0014 | 0.012 |
| mobility_no_tls | cl7 | entry | 47 | 0.026 / 0.015 | 0.504 / 0.435 | 0.136 / 0.083 | 0.538 / 0.486 | F 2.9 (df 13), p 0.009 | 0.041 | F 1.6 (df 4), p 0.2 | 0.38 |
| complete | cl4 | chain | 329 | 0.077 / 0.025 | 0.506 / 0.399 | 0.079 / 0.060 | 0.523 / 0.504 | cluster-robust Wald 69.0 (df 14), p 2.9e-09 | 6e-08 | cluster-robust Wald 19.8 (df 6), p 0.003 | 0.021 |
| fraction_observed | cl4 | chain | 329 | 0.031 / 0.007 | 0.189 / 0.146 | 0.057 / 0.031 | 0.196 / 0.169 | cluster-robust Wald 105.9 (df 14), p 3.5e-16 | 2.6e-14 | cluster-robust Wald 44.0 (df 7), p 2.1e-07 | 3.8e-06 |
| mobility | cl4 | chain | 325 | 0.054 / 0.006 | 0.201 / 0.082 | 0.143 / 0.066 | 0.210 / 0.133 | cluster-robust Wald 49.8 (df 14), p 6.5e-06 | 9.9e-05 | cluster-robust Wald 49.7 (df 8), p 4.6e-08 | 9e-07 |
| mobility_no_tls | cl4 | chain | 72 | 0.300 / 0.053 | 0.468 / 0.206 | 0.266 / 0.100 | 0.565 / 0.398 | cluster-robust Wald 216.8 (df 14), p 2e-38 | 7.2e-36 | cluster-robust Wald 22.9 (df 6), p 0.00084 | 0.009 |
| complete | cl4 | entry | 128 | 0.130 / 0.040 | 0.722 / 0.557 | 0.117 / 0.076 | 0.748 / 0.707 | LR 123.6 (df 14), p 1.2e-19 | 1.5e-17 | LR 13.3 (df 6), p 0.039 | 0.12 |
| fraction_observed | cl4 | entry | 128 | 0.133 / 0.023 | 0.361 / 0.195 | 0.118 / 0.047 | 0.378 / 0.308 | LR 26.2 (df 14), p 0.025 | 0.088 | LR 4.0 (df 7), p 0.78 | 0.85 |
| mobility | cl4 | entry | 127 | 0.126 / 0.010 | 0.227 / 0.172 | 0.192 / 0.130 | 0.276 / 0.213 | F 2.7 (df 14), p 0.0023 | 0.018 | F 2.8 (df 8), p 0.0069 | 0.035 |
| mobility_no_tls | cl4 | entry | 34 | 0.251 / 0.040 | 0.715 / 0.481 | 0.259 / 0.112 | 0.777 / 0.630 | F 5.2 (df 14), p 0.0025 | 0.018 | F 2.2 (df 6), p 0.11 | 0.25 |

BH: one family over every p-value of SV01b (368 tests).

Q+X gain entered first, complete outcome, SV01a (no lattice terms) -> SV01b: cl7 chain 0.224 -> 0.414; cl4 chain 0.208 -> 0.523; cl7 entry 0.613 -> 0.687; cl4 entry 0.346 -> 0.748; fraction observed: cl7 chain 0.272 -> 0.558; cl4 chain 0.067 -> 0.196; cl7 entry 0.676 -> 0.746; cl4 entry 0.209 -> 0.378.

Lattice terms (per 1 SD; odds ratio for complete / fraction observed, coefficient in z units for mobility):

| outcome | cluster | level | term | estimate [95% CI] | q |
|---|---|---|---|---|---|
| complete | cl7 | chain | asu | dropped: joint separation |  |
| complete | cl7 | chain | sym | 12.6 [2, 79.2] | 0.035 |
| fraction_observed | cl7 | chain | asu | dropped: joint separation |  |
| fraction_observed | cl7 | chain | sym | 189 [1.23, 2.89e+04] | 0.13 |
| mobility | cl7 | chain | asu | 0.0165 [-0.0711, 0.104] | 0.8 |
| mobility | cl7 | chain | sym | -0.17 [-0.226, -0.114] | 5.4e-08 |
| mobility_no_tls | cl7 | chain | asu | -0.0532 [-0.123, 0.0167] | 0.28 |
| mobility_no_tls | cl7 | chain | sym | -0.232 [-0.348, -0.116] | 0.0012 |
| complete | cl7 | entry | asu | dropped: joint separation |  |
| complete | cl7 | entry | sym | 5.64 [0.0574, 554] | 0.61 |
| fraction_observed | cl7 | entry | asu | dropped: joint separation |  |
| fraction_observed | cl7 | entry | sym | 5.13e+05 [0.0354, 7.44e+12] | 0.26 |
| mobility | cl7 | entry | asu | 0.0688 [-0.0142, 0.152] | 0.24 |
| mobility | cl7 | entry | sym | -0.113 [-0.181, -0.0438] | 0.012 |
| mobility_no_tls | cl7 | entry | asu | 0.0119 [-0.156, 0.179] | 0.91 |
| mobility_no_tls | cl7 | entry | sym | -0.0766 [-0.193, 0.0392] | 0.38 |
| complete | cl4 | chain | asu | 3.89 [2.26, 6.7] | 1.8e-05 |
| complete | cl4 | chain | sym | 20.2 [8.05, 50.9] | 4.1e-09 |
| fraction_observed | cl4 | chain | asu | 2.25 [1.77, 2.86] | 1.1e-09 |
| fraction_observed | cl4 | chain | sym | 5.03 [3.28, 7.72] | 5.9e-12 |
| mobility | cl4 | chain | asu | -0.027 [-0.16, 0.106] | 0.79 |
| mobility | cl4 | chain | sym | -0.268 [-0.407, -0.13] | 0.002 |
| mobility_no_tls | cl4 | chain | asu | -0.208 [-0.341, -0.0749] | 0.017 |
| mobility_no_tls | cl4 | chain | sym | -0.272 [-0.527, -0.0164] | 0.12 |
| complete | cl4 | entry | asu | 3.92 [1.15, 13.4] | 0.1 |
| complete | cl4 | entry | sym | 476 [16.8, 1.35e+04] | 0.0036 |
| fraction_observed | cl4 | entry | asu | 2.34 [0.448, 12.2] | 0.49 |
| fraction_observed | cl4 | entry | sym | 13.6 [1.42, 130] | 0.085 |
| mobility | cl4 | entry | asu | 0.0345 [-0.141, 0.21] | 0.8 |
| mobility | cl4 | entry | sym | -0.224 [-0.388, -0.0594] | 0.036 |
| mobility_no_tls | cl4 | entry | asu | 0.513 [-0.00475, 1.03] | 0.15 |
| mobility_no_tls | cl4 | entry | sym | -0.176 [-0.445, 0.0942] | 0.38 |

All effects with BH q < 0.05 (full model): complete cl4 chain lattice contacts, symmetry mate (atoms of 528-550) (X) 20.2 [8.05, 50.9], q 4.1e-09; complete cl4 chain lattice contacts, ASU partner (atoms of 528-550) (X) 3.89 [2.26, 6.7], q 1.8e-05; complete cl4 chain space group P21 (X) 140 [7.44, 2.63e+03], q 0.0097; complete cl4 chain ESR1 chains in ASU (X) 0.275 [0.124, 0.608], q 0.012; complete cl4 chain basic amine (B) 5.14 [1.61, 16.3], q 0.031; complete cl4 entry lattice contacts, symmetry mate (atoms of 528-550) (X) 476 [16.8, 1.35e+04], q 0.0036; complete cl4 entry ESR1 chains in ASU (X) 0.0116 [0.000363, 0.371], q 0.049; complete cl7 chain space group P21 (X) 443 [7.35, 2.67e+04], q 0.023; complete cl7 chain lattice contacts, symmetry mate (atoms of 528-550) (X) 12.6 [2, 79.2], q 0.035; complete cl7 chain ligand heavy atoms (B) 0.28 [0.111, 0.707], q 0.035; complete cl7 chain growth pH (X) 0.259 [0.0918, 0.729], q 0.045; fraction_observed cl4 chain lattice contacts, symmetry mate (atoms of 528-550) (X) 5.03 [3.28, 7.72], q 5.9e-12; fraction_observed cl4 chain lattice contacts, ASU partner (atoms of 528-550) (X) 2.25 [1.77, 2.86], q 1.1e-09; fraction_observed cl4 chain L536S (B) 8.19 [2.64, 25.4], q 0.0034; fraction_observed cl4 chain Y537S (B) 7.62 [2.51, 23.1], q 0.0038; fraction_observed cl4 chain mean B, entry (A^2) (Q) 1.4 [1.09, 1.78], q 0.036; fraction_observed cl7 chain ligand heavy atoms (B) 0.12 [0.031, 0.464], q 0.017; fraction_observed cl7 chain space group P21 (X) 1.54e+03 [9.58, 2.46e+05], q 0.027; fraction_observed cl7 chain growth pH (X) 0.147 [0.0376, 0.577], q 0.032; fraction_observed cl7 chain dimer partner by symmetry (X) 2.57e+03 [7.05, 9.39e+05], q 0.041; mobility cl4 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.268 [-0.407, -0.13], q 0.002; mobility cl4 chain space group P21 (X) -0.953 [-1.53, -0.375], q 0.012; mobility cl4 chain peptide bound (B) 1.35 [0.487, 2.22], q 0.017; mobility cl4 chain L536S (B) -0.722 [-1.21, -0.238], q 0.023; mobility cl4 chain other mutations (B) 0.468 [0.146, 0.791], q 0.026; mobility cl4 chain data-collection temperature (K) (X) -0.897 [-1.54, -0.25], q 0.035; mobility cl4 entry lattice contacts, symmetry mate (atoms of 528-550) (X) -0.224 [-0.388, -0.0594], q 0.036; mobility cl4 entry growth T not reported (X) -0.729 [-1.27, -0.191], q 0.037; mobility cl7 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.17 [-0.226, -0.114], q 5.4e-08; mobility cl7 chain data-collection temperature (K) (X) 0.157 [0.0925, 0.222], q 3.3e-05; mobility cl7 chain ESR1 chains in ASU (X) -0.199 [-0.286, -0.113], q 9.7e-05; mobility cl7 chain space group P21 (X) -0.668 [-1.07, -0.264], q 0.011; mobility cl7 chain mean B, entry (A^2) (Q) -0.101 [-0.179, -0.0227], q 0.049; mobility cl7 entry mean B, entry (A^2) (Q) -0.156 [-0.239, -0.0725], q 0.0031; mobility cl7 entry lattice contacts, symmetry mate (atoms of 528-550) (X) -0.113 [-0.181, -0.0438], q 0.012; mobility cl7 entry space group P21 (X) -0.659 [-1.08, -0.242], q 0.017; mobility cl7 entry dimer partner by symmetry (X) -0.782 [-1.33, -0.234], q 0.029; mobility cl7 entry growth pH (X) 0.115 [0.0314, 0.198], q 0.035; mobility cl7 entry carboxylic acid (B) 0.523 [0.145, 0.902], q 0.035; mobility_no_tls cl4 chain ligand heavy atoms (B) -0.465 [-0.74, -0.19], q 0.0097; mobility_no_tls cl4 chain lattice contacts, ASU partner (atoms of 528-550) (X) -0.208 [-0.341, -0.0749], q 0.017; mobility_no_tls cl4 entry ESR1 chains in ASU (X) -1.1 [-1.8, -0.403], q 0.017; mobility_no_tls cl4 entry growth T not reported (X) -1.19 [-2, -0.369], q 0.026; mobility_no_tls cl4 entry ligand heavy atoms (B) -0.438 [-0.774, -0.103], q 0.045; mobility_no_tls cl7 chain ligand heavy atoms (B) 0.181 [0.131, 0.231], q 3.4e-11; mobility_no_tls cl7 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.232 [-0.348, -0.116], q 0.0012; mobility_no_tls cl7 chain data-collection temperature (K) (X) 0.164 [0.0384, 0.289], q 0.045; mobility_no_tls cl7 entry dimer partner by symmetry (X) -1.03 [-1.7, -0.363], q 0.018; mobility_no_tls cl7 entry data-collection temperature (K) (X) 0.157 [0.042, 0.272], q 0.036.

Dropped predictors (cluster/level): complete cl7 chain: lattice_asu [joint separation], ambient_temp_K_missing [separation], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation] | fraction_observed cl7 chain: lattice_asu [joint separation], ambient_temp_K_missing [separation], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation] | mobility cl7 chain: partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [collinear with earlier predictors] | mobility_no_tls cl7 chain: ambient_temp_K_missing [constant], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [constant], basic_amine [constant], carboxylic_acid [constant] | complete cl7 entry: lattice_asu [joint separation], ambient_temp_K_missing [separation], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation] | fraction_observed cl7 entry: dimer_symmetry [joint separation], lattice_asu [joint separation], ambient_temp_K_missing [separation], peptide [joint separation], partner_peptide [collinear with earlier predictors], mut_Y537S [joint separation], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation], carboxylic_acid [joint separation] | mobility cl7 entry: partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [collinear with earlier predictors] | mobility_no_tls cl7 entry: ambient_temp_K_missing [constant], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [constant], basic_amine [constant], carboxylic_acid [constant] | complete cl4 chain: peptide [separation], partner_peptide [separation], mut_D538G [separation] | fraction_observed cl4 chain: peptide [separation], partner_peptide [separation] | mobility cl4 chain: partner_peptide [collinear with earlier predictors] | mobility_no_tls cl4 chain: peptide [constant], partner_peptide [constant], mut_D538G [constant] | complete cl4 entry: peptide [separation], partner_peptide [separation], mut_D538G [separation] | fraction_observed cl4 entry: peptide [separation], partner_peptide [separation] | mobility cl4 entry: partner_peptide [collinear with earlier predictors] | mobility_no_tls cl4 entry: peptide [constant], partner_peptide [constant], mut_D538G [constant].

Rows / predictors / full R2: complete cl7 chain n 445 (14 incomplete), 17 predictors, R2 0.508; fraction_observed cl7 chain n 445 (14 incomplete), 17 predictors, R2 0.707; mobility cl7 chain n 442, 20 predictors, R2 0.238; mobility_no_tls cl7 chain n 95, 17 predictors, R2 0.528; complete cl7 entry n 221 (5 incomplete), 17 predictors, R2 0.776; fraction_observed cl7 entry n 221 (5 incomplete), 13 predictors, R2 0.846; mobility cl7 entry n 220, 20 predictors, R2 0.244; mobility_no_tls cl7 entry n 47, 17 predictors, R2 0.621; complete cl4 chain n 329 (209 incomplete), 20 predictors, R2 0.583; fraction_observed cl4 chain n 329 (209 incomplete), 21 predictors, R2 0.226; mobility cl4 chain n 325, 22 predictors, R2 0.276; mobility_no_tls cl4 chain n 72, 20 predictors, R2 0.664; complete cl4 entry n 128 (55 incomplete), 20 predictors, R2 0.824; fraction_observed cl4 entry n 128 (55 incomplete), 21 predictors, R2 0.426; mobility cl4 entry n 127, 22 predictors, R2 0.405; mobility_no_tls cl4 entry n 34, 20 predictors, R2 0.888.

## Verdict (verdict.md)

# SV01b verdict: **artifact** (re-issued; outcomes H12 complete 538-547 and local mobility 528-547, lattice terms in block X; [DEFAULT] rule, chain and entry level)
- Complete (artifact; SV01a mixed): chain (artifact) cl7 Q+X 0.414 vs B 0.225, cl4 Q+X 0.523 vs B 0.079; entry (artifact) cl7 Q+X 0.687 vs B 0.315, cl4 Q+X 0.748 vs B 0.117 (McFadden gain entered first).
- Mobility (artifact): chain (artifact) cl7 Q+X 0.184 vs B 0.076, cl4 Q+X 0.210 vs B 0.143; entry (artifact) cl7 Q+X 0.157 vs B 0.098, cl4 Q+X 0.276 vs B 0.192 (OLS R2 gain entered first).
- Sensitivity: fraction observed **artifact** (chain artifact, entry artifact); mobility without TLS-refined entries **artifact** (chain artifact, entry artifact).
- Entered last, Q+X vs B gives the same ordering as entered first (same rows in every model: both reduce to R2(Q+X) vs R2(B)); values for both orders are in model_blocks.csv and verdict_table.csv.
- Complete, Q+X gain entered first without -> with the lattice terms (SV01a -> SV01b): cl7 chain 0.224 -> 0.414; cl4 chain 0.208 -> 0.523; cl7 entry 0.613 -> 0.687; cl4 entry 0.346 -> 0.748.
- Lattice terms, complete outcome, OR per 1 SD (chain level): cl7 asu dropped, sym 12.6 [2, 79.2]; cl4 asu 3.89 [2.26, 6.7], sym 20.2 [8.05, 50.9].
- Lattice terms, mobility, coefficient per 1 SD (chain level): cl7 asu 0.0165 [-0.0711, 0.104], sym -0.17 [-0.226, -0.114]; cl4 asu -0.027 [-0.16, 0.106], sym -0.268 [-0.407, -0.13].
- Mobility defined for 767/774 cl7+cl4 chains; of them in entries without TLS refinement 167.
- Lattice contacts of 528-550 are counted on modelled atoms, so they are partly a consequence of the completeness outcome; step 8 tests the medoid H12 positions instead (not run: verdict not mixed).

## 8. Medoid H12 clash test

Step 8 not run: re-issued SV01b verdict is 'artifact', not 'mixed'.

## Reading

- Adding the lattice terms moves the completeness verdict from mixed (SV01a) to artifact. Apart from the lattice terms the fitted predictors are those of SV01a (except: fraction_observed cl7 entry loses dimer_symmetry, mut_Y537S, peptide to joint separation), so the rise of the Q+X gain is the lattice terms' share.
- For completeness this share is partly mechanical: contacts are counted on modelled atoms. In cl4, complete chains have 24.3 mean symmetry-contact atoms on 528-550, partial chains 9.1; 12.3 of the former sit on 546-550, where partial chains model 2.9 heavy atoms on average (complete 30.5). On 538-545, which both model, the means are 8.22 vs 6.66. The data cannot tell whether the contact holds the C-terminal end in place or exists only because the end is modelled; step 8 would separate the two but its condition (verdict mixed) is not met.
- Mobility is z-scored on the residues that are modelled, so it is not tied to missing atoms in the same way: more symmetry-mate contact goes with lower local B relative to the core in both clusters at chain level (lattice table above), the direction expected if contacts restrain the segment. The B effects with q < 0.05 in cl4 mobility (chain level) are peptide bound (4 chains) +1.35, L536S (235 chains) -0.72, other mutations (290 chains) +0.47.
- Small samples: cl7 has 14 incomplete chains (5 entries) for 17 fitted predictors; the TLS-free mobility models at entry level fit 17 predictors on 47 cl7 entries (R2 0.62) and 20 predictors on 34 cl4 entries (R2 0.89), so their gains are optimistic.

## Files

- results/sv/SV01b/clash_test_not_run.md (0.1 kB)
- results/sv/SV01b/env.txt (0.3 kB)
- results/sv/SV01b/lattice_checks.json (0.6 kB)
- results/sv/SV01b/lattice_contacts_by_group.csv (2.5 kB)
- results/sv/SV01b/lattice_contacts_per_chain.csv (56.3 kB)
- results/sv/SV01b/lattice_contacts_per_residue.csv (714.0 kB)
- results/sv/SV01b/mobility_per_chain.csv (114.5 kB)
- results/sv/SV01b/mobility_per_residue.csv (513.7 kB)
- results/sv/SV01b/model_blocks.csv (48.6 kB)
- results/sv/SV01b/model_data.csv (136.6 kB)
- results/sv/SV01b/numbers.json (1.1 kB)
- results/sv/SV01b/odds_ratios.csv (102.9 kB)
- results/sv/SV01b/predictor_scaling.csv (1.5 kB)
- results/sv/SV01b/reproduce_SV01a.json (0.3 kB)
- results/sv/SV01b/verdict.md (1.7 kB)
- results/sv/SV01b/verdict_table.csv (3.4 kB)
- figures/sv/SV01b/SV01b_forest_completeness.png (350.6 kB)
- figures/sv/SV01b/SV01b_forest_completeness.svg (164.6 kB)
- figures/sv/SV01b/SV01b_forest_mobility.png (333.2 kB)
- figures/sv/SV01b/SV01b_forest_mobility.svg (131.8 kB)
- figures/panels/sv/SV01b/SV01b_forest_completeness__forest_cl4.png (157.6 kB)
- figures/panels/sv/SV01b/SV01b_forest_completeness__forest_cl7.png (161.5 kB)
- figures/panels/sv/SV01b/SV01b_forest_completeness__gains_cl4.png (53.4 kB)
- figures/panels/sv/SV01b/SV01b_forest_completeness__gains_cl7.png (53.0 kB)
- figures/panels/sv/SV01b/SV01b_forest_mobility__forest_cl4.png (144.9 kB)
- figures/panels/sv/SV01b/SV01b_forest_mobility__forest_cl7.png (147.1 kB)
- figures/panels/sv/SV01b/SV01b_forest_mobility__gains_cl4.png (53.0 kB)
- figures/panels/sv/SV01b/SV01b_forest_mobility__gains_cl7.png (52.9 kB)
- scripts/sv/SV01b_lattice.py (4.1 kB)
- scripts/sv/SV01b_models.py (30.2 kB)
- scripts/sv/SV01b_clash.py (10.9 kB)
- scripts/sv/SV01b_figures.py (6.4 kB)
- scripts/sv/SV01b_summary.py (16.9 kB)
- scripts/sv/common/lattice_contacts.py (12.2 kB)
- docs/EXPLORATION.md: appended 'Suggested, not done (SV01b, 2026-11-08)'
- results/sv/SV01b/SUMMARY.md (this file)
