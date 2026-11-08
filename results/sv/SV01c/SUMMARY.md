Window propagation, fixed before fitting (docstring of scripts/sv/SV01c_window.py): H12 complete = all 7 C-alpha of 538-544 observed (fraction observed over 538-544; status complete / partial / not modelled re-derived on 538-544); mobility = mean C-alpha B z-score over the observed residues of 528-544 (SV01b: 528-547; loop 528-537, core v3 reference, B and altloc rules unchanged); lattice terms keep 528-550, the SV01b counts reused unchanged (step 6 defines them on 528-550, not on the H12 window). Not done / substituted: no figure (step 9 asks for the verdict side by side, given as tables); step 8 run because the SV01c verdict is mixed, with scripts/sv/SV01b_clash.py unchanged except that the distance it reports is taken to the nearest lattice image (the shared lattice_contacts.lattice_hits reports the distance to the neighbour's stored image, which is wrong for 196/1058 hits found across a lattice translation, up to 146.9 A; hit sets identical, so clash counts and 'blocked' are those of the shared code; the module is not edited); everything else inherited unchanged from SV01b (blocks Q/X/B with lattice terms in X, median imputation + 'not reported' indicators, drop rule, McFadden / OLS R2 gains, block tests, one BH family over the 338 SV01c p-values, [DEFAULT] verdict rule with the SV01b step rule).

# SV01c - H12 completeness and mobility on the data-derived H12 window 538-544

**SV01c verdict (538-544): mixed; primary SV01b verdict (538-547): artifact.** Complete H12: artifact -> mixed (chain artifact -> mixed, entry artifact -> artifact); mobility: artifact -> artifact (chain artifact -> artifact, entry artifact -> artifact). Sensitivity: fraction observed artifact -> mixed, mobility without TLS-refined entries artifact -> artifact. Step 8 run (verdict mixed).

## Inputs

- SV02 (DONE): data-derived H12 window 538-544 = results/sv/SV02/helix_boundaries.csv row H12 (global consensus; proposal, not adopted; results/sv/SV02/h12_window.md).
- SV01b (DONE): scripts/sv/SV01b_models.py imported with SV01a_completeness.H12 and SV01b_models.MOB set to the window (no file edited); lattice terms from results/sv/SV01b/lattice_contacts_per_chain.csv (scripts/sv/common/lattice_contacts.py output, not recomputed).
- Check (reproduce_SV01b.json): the same code on 538-547 reproduces results/sv/SV01b: 128/128 block rows, max |d gain| 9.7e-17, max |d q| 9.9e-17; 368/368 effect rows, status equal 368, max |d coef| 1.8e-15; BH tests 368 = 368; mobility z for 788 chains, max |d| 4.4e-16. Only the window differs.
- C-alpha B cross-check (backbone.npz vs residue mapping): 181829 C-alpha, max |d| 0.0.

## Outcomes on the two windows (chain_outcomes.csv, status_transitions.csv)

| group | level | 538-547 status | 538-544 status | n |
|---|---|---|---|---|
| cl4 | chain | complete | complete | 120 |
| cl4 | chain | not_modelled | not_modelled | 14 |
| cl4 | chain | partial | complete | 190 |
| cl4 | chain | partial | partial | 5 |
| cl4 | entry | complete | complete | 73 |
| cl4 | entry | not_modelled | not_modelled | 6 |
| cl4 | entry | partial | complete | 49 |
| cl7 | chain | complete | complete | 431 |
| cl7 | chain | not_modelled | not_modelled | 7 |
| cl7 | chain | partial | complete | 4 |
| cl7 | chain | partial | partial | 3 |
| cl7 | entry | complete | complete | 216 |
| cl7 | entry | not_modelled | not_modelled | 4 |
| cl7 | entry | partial | complete | 1 |
| minor | chain | complete | complete | 2 |
| minor | chain | not_modelled | not_modelled | 2 |
| minor | chain | partial | complete | 10 |
| minor | entry | complete | complete | 1 |
| minor | entry | not_modelled | not_modelled | 1 |
| minor | entry | partial | complete | 5 |
| all | chain | complete | complete | 553 |
| all | chain | not_modelled | not_modelled | 23 |
| all | chain | partial | complete | 204 |
| all | chain | partial | partial | 8 |
| all | entry | complete | complete | 289 |
| all | entry | not_modelled | not_modelled | 8 |
| all | entry | partial | complete | 54 |

Incomplete H12 (model events) 538-547 -> 538-544: cl7 14 -> 10 chains (5 -> 4 entries at entry level); cl4 209 -> 19 chains (55 -> 6 entries). Mobility defined for cl7 442, cl4 325 chains (no TLS: 95, 72); undefined: 5T92|A|1, 5T97|B|1, 5TM9|A|1, 6B0F|B|1, 7RS3|B|1, 7RS3|C|1, 7YMK|C|1.

## Verdict side by side (side_by_side.csv; gain = McFadden R2 for complete / fraction observed, OLS R2 for mobility)

| outcome | level | cluster | order | Q+X 538-547 | B 538-547 | Q+X 538-544 | B 538-544 | level verdict 538-547 -> 538-544 |
|---|---|---|---|---|---|---|---|---|
| complete | chain | cl7 | first | 0.414 | 0.225 | 0.130 | 0.386 | artifact -> mixed |
| complete | chain | cl7 | last | 0.283 | 0.095 | 0.072 | 0.327 | artifact -> mixed |
| complete | chain | cl4 | first | 0.523 | 0.079 | 0.506 | 0.216 | artifact -> mixed |
| complete | chain | cl4 | last | 0.504 | 0.060 | 0.427 | 0.137 | artifact -> mixed |
| complete | entry | cl7 | first | 0.687 | 0.315 | 0.549 | 0.167 | artifact -> artifact |
| complete | entry | cl7 | last | 0.461 | 0.089 | 0.598 | 0.216 | artifact -> artifact |
| complete | entry | cl4 | first | 0.748 | 0.117 | 0.416 | 0.214 | artifact -> artifact |
| complete | entry | cl4 | last | 0.707 | 0.076 | 0.308 | 0.106 | artifact -> artifact |
| fraction_observed | chain | cl7 | first | 0.558 | 0.352 | 0.272 | 0.443 | artifact -> mixed |
| fraction_observed | chain | cl7 | last | 0.356 | 0.149 | 0.219 | 0.391 | artifact -> mixed |
| fraction_observed | chain | cl4 | first | 0.196 | 0.057 | 0.646 | 0.354 | artifact -> mixed |
| fraction_observed | chain | cl4 | last | 0.169 | 0.031 | 0.483 | 0.191 | artifact -> mixed |
| fraction_observed | entry | cl7 | first | 0.746 | 0.156 | 0.549 | 0.167 | artifact -> artifact |
| fraction_observed | entry | cl7 | last | 0.690 | 0.100 | 0.598 | 0.216 | artifact -> artifact |
| fraction_observed | entry | cl4 | first | 0.378 | 0.118 | 0.416 | 0.214 | artifact -> artifact |
| fraction_observed | entry | cl4 | last | 0.308 | 0.047 | 0.308 | 0.106 | artifact -> artifact |
| mobility | chain | cl7 | first | 0.184 | 0.076 | 0.215 | 0.081 | artifact -> artifact |
| mobility | chain | cl7 | last | 0.162 | 0.054 | 0.192 | 0.058 | artifact -> artifact |
| mobility | chain | cl4 | first | 0.210 | 0.143 | 0.200 | 0.156 | artifact -> artifact |
| mobility | chain | cl4 | last | 0.133 | 0.066 | 0.117 | 0.073 | artifact -> artifact |
| mobility | entry | cl7 | first | 0.157 | 0.098 | 0.180 | 0.124 | artifact -> artifact |
| mobility | entry | cl7 | last | 0.146 | 0.086 | 0.164 | 0.108 | artifact -> artifact |
| mobility | entry | cl4 | first | 0.276 | 0.192 | 0.249 | 0.205 | artifact -> artifact |
| mobility | entry | cl4 | last | 0.213 | 0.130 | 0.180 | 0.136 | artifact -> artifact |
| mobility_no_tls | chain | cl7 | first | 0.419 | 0.144 | 0.435 | 0.130 | artifact -> artifact |
| mobility_no_tls | chain | cl7 | last | 0.384 | 0.109 | 0.389 | 0.084 | artifact -> artifact |
| mobility_no_tls | chain | cl4 | first | 0.565 | 0.266 | 0.551 | 0.243 | artifact -> artifact |
| mobility_no_tls | chain | cl4 | last | 0.398 | 0.100 | 0.402 | 0.094 | artifact -> artifact |
| mobility_no_tls | entry | cl7 | first | 0.538 | 0.136 | 0.526 | 0.121 | artifact -> artifact |
| mobility_no_tls | entry | cl7 | last | 0.486 | 0.083 | 0.462 | 0.057 | artifact -> artifact |
| mobility_no_tls | entry | cl4 | first | 0.777 | 0.259 | 0.729 | 0.230 | artifact -> artifact |
| mobility_no_tls | entry | cl4 | last | 0.630 | 0.112 | 0.618 | 0.118 | artifact -> artifact |

## Models on 538-544 (model_blocks.csv, odds_ratios.csv)

| outcome | cluster | level | n | events | fitted predictors | full R2 | test Q+X (last) | q Q+X | test B (last) | q B |
|---|---|---|---|---|---|---|---|---|---|---|
| complete | cl7 | chain | 445 | 10 | 15 | 0.458 | Wald 15.1 (df 10), p 0.127 | 0.306 | Wald 18.4 (df 5), p 0.00252 | 0.0258 |
| fraction_observed | cl7 | chain | 445 | 10 | 15 | 0.663 | Wald 15.5 (df 9), p 0.0786 | 0.221 | Wald 10.6 (df 6), p 0.102 | 0.266 |
| mobility | cl7 | chain | 442 | - | 20 | 0.273 | Wald 97.3 (df 14), p 1.54e-14 | 1.58e-12 | Wald 21.9 (df 6), p 0.00125 | 0.0162 |
| mobility_no_tls | cl7 | chain | 95 | - | 17 | 0.519 | Wald 127.0 (df 13), p 8.23e-21 | 1.39e-18 | Wald 51.4 (df 4), p 1.88e-10 | 9.1e-09 |
| complete | cl7 | entry | 221 | 4 | 13 | 0.765 | LR 23.9 (df 10), p 0.00782 | 0.0562 | LR 8.6 (df 3), p 0.0346 | 0.138 |
| fraction_observed | cl7 | entry | 221 | 4 | 13 | 0.765 | LR 23.9 (df 10), p 0.00782 | 0.0562 | LR 8.6 (df 3), p 0.0346 | 0.138 |
| mobility | cl7 | entry | 220 | - | 20 | 0.288 | F 3.3 (df 14), p 0.0001 | 0.00188 | F 5.0 (df 6), p 8.31e-05 | 0.00165 |
| mobility_no_tls | cl7 | entry | 47 | - | 17 | 0.583 | F 2.5 (df 13), p 0.0209 | 0.0979 | F 1.0 (df 4), p 0.431 | 0.617 |
| complete | cl4 | chain | 329 | 19 | 18 | 0.643 | Wald 11.1 (df 12), p 0.521 | 0.669 | Wald 25.2 (df 6), p 0.000313 | 0.00473 |
| fraction_observed | cl4 | chain | 329 | 19 | 18 | 0.837 | Wald 22.8 (df 12), p 0.0291 | 0.12 | Wald 17.3 (df 6), p 0.00824 | 0.058 |
| mobility | cl4 | chain | 325 | - | 22 | 0.273 | Wald 44.8 (df 14), p 4.34e-05 | 0.00105 | Wald 54.9 (df 8), p 4.69e-09 | 1.98e-07 |
| mobility_no_tls | cl4 | chain | 72 | - | 20 | 0.645 | Wald 190.4 (df 14), p 4.98e-33 | 1.68e-30 | Wald 19.1 (df 6), p 0.00407 | 0.0344 |
| complete | cl4 | entry | 128 | 6 | 12 | 0.522 | LR 14.9 (df 8), p 0.0609 | 0.186 | LR 5.1 (df 4), p 0.275 | 0.474 |
| fraction_observed | cl4 | entry | 128 | 6 | 12 | 0.522 | LR 14.9 (df 8), p 0.0609 | 0.186 | LR 5.1 (df 4), p 0.275 | 0.474 |
| mobility | cl4 | entry | 127 | - | 22 | 0.385 | F 2.2 (df 14), p 0.0134 | 0.0807 | F 2.9 (df 8), p 0.00608 | 0.0467 |
| mobility_no_tls | cl4 | entry | 34 | - | 20 | 0.848 | F 3.8 (df 14), p 0.011 | 0.07 | F 1.7 (df 6), p 0.202 | 0.409 |

Dropped predictors, complete outcome (cluster / level: predictor [reason]): cl7 chain: p21 [joint separation], lattice_asu [joint separation], lattice_sym [joint separation], ambient_temp_K_missing [separation], partner_peptide [collinear with earlier predictors], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation] | cl7 entry: dimer_symmetry [joint separation], lattice_asu [joint separation], lattice_sym [joint separation], ambient_temp_K_missing [separation], partner_peptide [collinear with earlier predictors], mut_Y537S [joint separation], mut_L536S [constant], mut_D538G [separation], basic_amine [joint separation], carboxylic_acid [joint separation] | cl4 chain: dimer_symmetry [separation], lattice_asu [joint separation], peptide [separation], partner_peptide [separation], mut_D538G [separation] | cl4 entry: esr1_chains_asu [joint separation], dimer_symmetry [separation], lattice_asu [joint separation], lattice_sym [joint separation], ph_missing [separation], ambient_temp_K_missing [separation], peptide [separation], partner_peptide [separation], mut_Y537S [separation], mut_L536S [separation], mut_D538G [separation].

Effects with BH q < 0.05 (full model; odds ratio per 1 SD or present vs absent for complete / fraction observed, coefficient in z units for mobility): complete cl4 chain Y537S (B) 283 [14.8, 5.39e+03], q 0.0031; complete cl4 chain ligand heavy atoms (B) 0.145 [0.0508, 0.414], q 0.00473; complete cl7 chain ligand heavy atoms (B) 0.135 [0.0354, 0.512], q 0.0297; fraction_observed cl4 chain ligand heavy atoms (B) 0.0285 [0.00368, 0.221], q 0.00929; fraction_observed cl7 chain basic amine (B) 0.00148 [2.1e-05, 0.104], q 0.0265; mobility cl4 chain space group P21 (X) -1.07 [-1.71, -0.435], q 0.0134; mobility cl4 chain peptide bound (B) 1.42 [0.54, 2.29], q 0.0186; mobility cl4 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.219 [-0.361, -0.0771], q 0.0258; mobility cl4 chain data-collection temperature (K) (X) -1.01 [-1.7, -0.324], q 0.0342; mobility cl4 chain L536S (B) -0.737 [-1.24, -0.23], q 0.0358; mobility cl4 chain other mutations (B) 0.524 [0.158, 0.89], q 0.0404; mobility cl7 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.232 [-0.291, -0.172], q 1.58e-12; mobility cl7 chain data-collection temperature (K) (X) 0.154 [0.0886, 0.219], q 0.000101; mobility cl7 chain ESR1 chains in ASU (X) -0.183 [-0.27, -0.0966], q 0.000869; mobility cl7 chain space group P21 (X) -0.618 [-1.03, -0.205], q 0.0301; mobility cl7 entry lattice contacts, symmetry mate (atoms of 528-550) (X) -0.173 [-0.242, -0.104], q 2.58e-05; mobility cl7 entry mean B, entry (A^2) (Q) -0.136 [-0.22, -0.0525], q 0.0177; mobility cl7 entry carboxylic acid (B) 0.601 [0.222, 0.981], q 0.0215; mobility cl7 entry space group P21 (X) -0.645 [-1.06, -0.227], q 0.0258; mobility cl7 entry growth pH (X) 0.127 [0.0438, 0.211], q 0.0273; mobility_no_tls cl4 chain lattice contacts, ASU partner (atoms of 528-550) (X) -0.242 [-0.373, -0.11], q 0.00473; mobility_no_tls cl4 chain ligand heavy atoms (B) -0.484 [-0.787, -0.182], q 0.02; mobility_no_tls cl4 entry ESR1 chains in ASU (X) -1.27 [-2.11, -0.428], q 0.0295; mobility_no_tls cl7 chain ligand heavy atoms (B) 0.176 [0.115, 0.237], q 5.59e-07; mobility_no_tls cl7 chain lattice contacts, symmetry mate (atoms of 528-550) (X) -0.282 [-0.418, -0.146], q 0.00109.

Lattice contacts of 528-550 by H12 status on the new window (lattice_contacts_by_group.csv; mean heavy atoms in contact):

| group | H12 (538-544) | chains | entries | any ASU | any sym | mean ASU atoms | mean sym atoms | modelled heavy atoms 528-550 |
|---|---|---|---|---|---|---|---|---|
| cl4 | complete | 310 | 123 | 62 | 267 | 1.81 | 15.09 | 126.2 |
| cl4 | not_modelled | 14 | 10 | 0 | 3 | 0.00 | 0.36 | 6.8 |
| cl4 | partial | 5 | 4 | 0 | 4 | 0.00 | 3.00 | 65.6 |
| cl7 | complete | 435 | 217 | 293 | 329 | 1.14 | 5.83 | 154.9 |
| cl7 | not_modelled | 7 | 6 | 0 | 0 | 0.00 | 0.00 | 6.3 |
| cl7 | partial | 3 | 3 | 0 | 0 | 0.00 | 0.00 | 29.7 |
| minor | complete | 12 | 6 | 4 | 10 | 8.17 | 12.92 | 136.2 |
| minor | not_modelled | 2 | 2 | 0 | 1 | 0.00 | 3.00 | 15.5 |
| cl4 | all | 329 | 128 | 62 | 274 | 1.71 | 14.28 | 120.2 |
| cl7 | all | 445 | 221 | 293 | 329 | 1.12 | 5.70 | 151.7 |
| minor | all | 14 | 7 | 4 | 11 | 7.00 | 11.50 | 119.0 |

## 8. Medoid H12 clash test on 538-544 (clash_test.csv, clash_summary.csv, clash_checks.json)

- Targets: 31 chains partial or not modelled on 538-544 (all clusters: cl4 not_modelled 14, cl4 partial 5, cl7 not_modelled 7, cl7 partial 3, minor not_modelled 2); moved atoms = medoid residues 538-544 (cl7 57, cl4 54 heavy atoms); clash = 3.0 A; blocked = >= 1 backbone atom in clash.
- Core fit: >= 112 core v3 C-alpha per fit, RMSD median 0.48 A, max 1.56 A; direct fit vs SV00 gate transforms on the medoid H12 C-alpha: median 0.08 A, max 0.39 A (n 62).

| group | H12 | chains | cl7 position blocked | cl4 position blocked | only cl7 | only cl4 | both | neither |
|---|---|---|---|---|---|---|---|---|
| cl4 | not_modelled | 14 | 2 | 2 | 0 | 0 | 2 | 12 |
| cl4 | partial | 5 | 0 | 0 | 0 | 0 | 0 | 5 |
| cl7 | not_modelled | 7 | 0 | 0 | 0 | 0 | 0 | 7 |
| cl7 | partial | 3 | 0 | 0 | 0 | 0 | 0 | 3 |
| minor | not_modelled | 2 | 1 | 1 | 0 | 0 | 1 | 1 |
| all | all | 31 | 3 | 3 | 0 | 0 | 3 | 28 |

Blocked chain x medoid rows: 2BJ4|A|1 (cl4, not_modelled) cl7 medoid: 1 backbone / 5 heavy atoms in clash, min 0.78 A, partner B@1 (receptor, symmetry); 2BJ4|A|1 (cl4, not_modelled) cl4 medoid: 16 backbone / 28 heavy atoms in clash, min 0.47 A, partner B@1 (receptor, symmetry); 2BJ4|B|1 (cl4, not_modelled) cl7 medoid: 19 backbone / 37 heavy atoms in clash, min 0.22 A, partner B@2 (receptor, symmetry); 2BJ4|B|1 (cl4, not_modelled) cl4 medoid: 11 backbone / 22 heavy atoms in clash, min 0.28 A, partner B@2 (receptor, symmetry); 6PFM|B|1 (minor, not_modelled) cl7 medoid: 20 backbone / 38 heavy atoms in clash, min 0.30 A, partner D@16 (receptor, symmetry); 6PFM|B|1 (minor, not_modelled) cl4 medoid: 28 backbone / 54 heavy atoms in clash, min 0.19 A, partner D@16 (receptor, symmetry).

## Reading

- On 538-544, 204 chains that are partial on 538-547 become complete (cl4 190, cl7 4, minor 10); what stays incomplete is mostly H12 not modelled at all (cl7 7 of 10, cl4 14 of 19 incomplete chains). The completeness outcome on the new window is therefore a different, much rarer event than on the primary window.
- The completeness verdict changes from artifact to mixed at chain level only: in cl7 the B gain exceeds the Q+X gain in both orders, in cl4 the reverse. In cl7 (chain level) both lattice terms and space group P21 are dropped by joint separation (not-modelled cl7 chains have no lattice contact on 528-550: 0 of 7 with a symmetry contact), so block X loses the terms that carried its SV01b gain; the cl7 B block is then led by ligand heavy atoms, OR 0.135 per 1 SD [0.0354, 0.512], q 0.0297.
- The events are few: cl7 10 incomplete chains (4 entries at entry level), cl4 19 chains (6 entries), against 12-18 fitted predictors (table above), so the gains on the new window are unstable and the entry-level 'artifact' rests on 4 (cl7) and 6 (cl4) incomplete entries. At entry level every incomplete row is H12 not modelled, so fraction observed gives the same fit as the binary outcome there.
- The mobility verdict does not depend on the window: it stays artifact at both levels and without TLS-refined entries; dropping 545-547 from the mobility mean changes the Q+X and B gains by at most 0.047 (side-by-side table).
- Step 8: of 31 partial / not-modelled chains on 538-544, 3 have both medoid H12 positions blocked (by symmetry copies of receptor chains), 0 only the cl7 position, 0 only the cl4 position, and 28 neither. The crystal lattice leaves room for either H12 position in most chains whose H12 is missing on the new window, and where it blocks, it blocks both, so it does not favour one position over the other.
- Wording: crystal structures show H12 order and position (modelled or not, local B), not receptor activity.

## Verdict (verdict.md)

# SV01c verdict (H12 window 538-544, SV02 proposal): **mixed**; primary SV01b (538-547): **artifact** ([DEFAULT] rule, lattice terms in X, chain and entry level)
- Complete H12: artifact -> **mixed** (chain artifact -> mixed, entry artifact -> artifact).
- Mobility (528-547 -> 528-544): artifact -> **artifact** (chain artifact -> artifact, entry artifact -> artifact).
- Sensitivity: fraction observed artifact -> mixed; mobility without TLS artifact -> artifact.
- Complete, gain entered first 538-547 -> 538-544, chain: cl7 Q+X 0.414 vs B 0.225 -> 0.130 vs 0.386; cl4 Q+X 0.523 vs B 0.079 -> 0.506 vs 0.216.
- Complete, entry: cl7 Q+X 0.687 vs B 0.315 -> 0.549 vs 0.167; cl4 Q+X 0.748 vs B 0.117 -> 0.416 vs 0.214.
- Mobility, chain: cl7 Q+X 0.184 vs B 0.076 -> 0.215 vs 0.081; cl4 Q+X 0.210 vs B 0.143 -> 0.200 vs 0.156; entry: cl7 Q+X 0.157 vs B 0.098 -> 0.180 vs 0.124; cl4 Q+X 0.276 vs B 0.192 -> 0.249 vs 0.205.
- Entered last gives the same Q+X vs B ordering as entered first (same rows in every model); both orders in side_by_side.csv.
- Lattice terms keep 528-550 (SV01b counts reused); only the outcomes moved. Step 8 (verdict mixed): of 31 chains partial / not modelled on 538-544, medoid H12 blocked at both positions 3, only cl7 0, only cl4 0, neither 28 (clash_summary.csv).

## Files

- results/sv/SV01c/verdict.md (1.3 kB)
- results/sv/SV01c/side_by_side.csv (5.1 kB)
- results/sv/SV01c/verdict_table.csv (3.3 kB)
- results/sv/SV01c/model_blocks.csv (46.9 kB)
- results/sv/SV01c/odds_ratios.csv (100.1 kB)
- results/sv/SV01c/chain_outcomes.csv (76.8 kB)
- results/sv/SV01c/status_transitions.csv (1.0 kB)
- results/sv/SV01c/mobility_per_chain.csv (114.2 kB)
- results/sv/SV01c/lattice_contacts_by_group.csv (2.3 kB)
- results/sv/SV01c/predictor_scaling.csv (1.5 kB)
- results/sv/SV01c/model_data.csv (136.7 kB)
- results/sv/SV01c/reproduce_SV01b.json (0.6 kB)
- results/sv/SV01c/numbers.json (1.9 kB)
- results/sv/SV01c/env.txt (0.4 kB)
- results/sv/SV01c/clash_test.csv (8.1 kB)
- results/sv/SV01c/clash_summary.csv (2.2 kB)
- results/sv/SV01c/clash_checks.json (1.2 kB)
- scripts/sv/SV01c_window.py (15.6 kB)
- scripts/sv/SV01c_clash.py (6.9 kB)
- scripts/sv/SV01c_summary.py (17.4 kB)
- docs/EXPLORATION.md: appended 'Suggested, not done (SV01c, 2026-11-08)'
- results/sv/SV01c/SUMMARY.md (this file)
