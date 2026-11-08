Not done / substituted: (1) docs/sv/briefs/SV01a.md carried rules 1-3 but not the SV01 text, so the step text was read from docs/sv/RUN_SV1028.md section 4 (SV01); (2) the run file gives no rule for missing predictor values: pre-declared median imputation within the cluster plus a 'not reported' indicator in block X (pH missing for 44 cl7 and 234 cl4 chains), so every model of a cluster uses the same rows; (3) 'mean B' = the entry mean B of Q05/K05 (the QC-tier and Y14 value), not the chain C-alpha mean B, which averages over the outcome residues; (4) entry level = one chain per entry (the entry's first matrix chain in the cluster, Y14 convention); (5) the step's drop rule was extended, before fitting, to exact collinearity and to joint separation in the full model (drop the predictor with the largest |log OR| while the fit separates); (6) the [DEFAULT] verdict rule names no level, so it was applied at chain and at entry level and the step verdict is 'artifact'/'biology' only when both levels agree; (7) chain-level errors: logit MLE with entry-clustered sandwich errors (= GEE, independence working correlation); the fraction outcome uses the fractional-logit QMLE (statsmodels Logit). Lattice contacts and the mobility outcome belong to SV01b and were not attempted here.

# SV01a - H12 completeness, first pass

**Verdict ([DEFAULT] rule, H12 complete on 538-547): mixed** - chain level mixed, entry level artifact. Sensitivity (fraction of 538-547 observed): mixed (chain mixed, entry artifact).

## Inputs

- results/sv/SV00/chain_table.csv and observed_mask.npy (tracker card SV00 out; SV00 DONE in tracker/esr1.log); 788 chains, 351 entries, remediated chains excluded: 0.
- Outcome from the mask (all 10 C-alpha of 538-547 observed) agrees with the S16 h12_label for 788/788 chains.

## 1. Distribution

H12 status (h12_status_counts.csv; entry level = first matrix chain of each entry in the group):

| level | group | n | complete | partial | not modelled | fraction complete |
|---|---|---|---|---|---|---|
| chain | cl7 | 445 | 431 | 7 | 7 | 0.969 |
| chain | cl4 | 329 | 120 | 195 | 14 | 0.365 |
| chain | minor | 14 | 2 | 10 | 2 | 0.143 |
| chain | Cys->Ser/Ala only | 30 | 17 | 4 | 9 | 0.567 |
| chain | D538G | 12 | 4 | 8 | 0 | 0.333 |
| chain | L536S | 243 | 81 | 159 | 3 | 0.333 |
| chain | Y537S | 443 | 407 | 32 | 4 | 0.919 |
| chain | other | 2 | 0 | 2 | 0 | 0.000 |
| chain | wild type | 58 | 44 | 7 | 7 | 0.759 |
| chain | all | 788 | 553 | 212 | 23 | 0.702 |
| entry | cl7 | 221 | 216 | 1 | 4 | 0.977 |
| entry | cl4 | 128 | 73 | 49 | 6 | 0.570 |
| entry | minor | 7 | 1 | 5 | 1 | 0.143 |
| entry | Cys->Ser/Ala only | 16 | 10 | 2 | 4 | 0.625 |
| entry | D538G | 3 | 2 | 1 | 0 | 0.667 |
| entry | L536S | 84 | 44 | 40 | 0 | 0.524 |
| entry | Y537S | 216 | 208 | 7 | 1 | 0.963 |
| entry | other | 1 | 0 | 1 | 0 | 0.000 |
| entry | wild type | 31 | 25 | 3 | 3 | 0.806 |
| entry | all | 351 | 289 | 54 | 8 | 0.823 |

Within cluster by construct (chain level, fraction complete): cl7|Cys->Ser/Ala only 5/6; cl7|D538G 4/4; cl7|Y537S 394/405; cl7|wild type 28/30; cl4|Cys->Ser/Ala only 12/24; cl4|D538G 0/8; cl4|L536S 80/237; cl4|Y537S 13/37; cl4|wild type 15/23.

Partial H12 (212 chains, partial_ends.csv): C-side 207, both ends 5; by cluster cl4: C-side 193, both ends 2; cl7: C-side 4, both ends 3; minor: C-side 10. No partial chain lacks only N-side residues. Missing at the C end: median 1, max 5 residues; internal gaps in 1 chain(s). Last observed residue of partial chains: 542: 1, 543: 3, 544: 24, 545: 65, 546: 119.

Per-residue observation (completeness_profile.csv, chain fraction): cl7: 538 0.98, 545 0.98, 546 0.98, 547 0.97, 550 0.09; cl4: 538 0.95, 545 0.89, 546 0.71, 547 0.36, 550 0.17. H11-H12 loop minimum (525-537): cl4 0.35 at 530; cl7 0.82 at 533.

## 2. Within-entry concordance (concordance.csv, concordance_summary.csv)

| measure | n | k | fraction |
|---|---|---|---|
| multi-chain entries | 317 | 317 | 1.000 |
| H12 label differs (complete/partial/not modelled) | 317 | 59 | 0.186 |
| binary outcome differs (complete vs not) | 317 | 56 | 0.177 |
| S08 cluster differs | 317 | 5 | 0.016 |
| label and cluster both differ | 317 | 1 | 0.003 |
| label differs between chains of the same cluster (entries with >=2 chains in one cluster) | 313 | 59 | 0.188 |
| cl7: binary outcome differs between the entry's cl7 chains (entries with >=2 cl7 chains) | 212 | 7 | 0.033 |
| cl4: binary outcome differs between the entry's cl4 chains (entries with >=2 cl4 chains) | 97 | 48 | 0.495 |

## 3. Block models (model_blocks.csv, odds_ratios.csv)

McFadden R2 gain, outcome H12 complete (Q quality, X experiment, B biology; Q+X combined):

| cluster | level | n | incomplete | block | first | last | test (entered last) | q (BH) |
|---|---|---|---|---|---|---|---|---|
| cl7 | chain | 445 | 14 | Q | 0.044 | 0.056 | cluster-robust Wald 7.9 (df 3), p 0.048 | 0.16 |
| cl7 | chain | 445 | 14 | X | 0.175 | 0.129 | cluster-robust Wald 13.3 (df 8), p 0.1 | 0.25 |
| cl7 | chain | 445 | 14 | B | 0.225 | 0.155 | cluster-robust Wald 13.7 (df 5), p 0.018 | 0.088 |
| cl7 | chain | 445 | 14 | Q+X | 0.224 | 0.153 | cluster-robust Wald 14.8 (df 11), p 0.19 | 0.4 |
| cl7 | entry | 221 | 5 | Q | 0.261 | 0.276 | LR 13.2 (df 3), p 0.0042 | 0.026 |
| cl7 | entry | 221 | 5 | X | 0.247 | 0.230 | LR 11.0 (df 8), p 0.2 | 0.42 |
| cl7 | entry | 221 | 5 | B | 0.315 | 0.142 | LR 6.8 (df 5), p 0.24 | 0.46 |
| cl7 | entry | 221 | 5 | Q+X | 0.613 | 0.440 | LR 21.0 (df 11), p 0.033 | 0.13 |
| cl4 | chain | 329 | 209 | Q | 0.077 | 0.080 | cluster-robust Wald 17.3 (df 3), p 0.00061 | 0.0074 |
| cl4 | chain | 329 | 209 | X | 0.122 | 0.096 | cluster-robust Wald 67.6 (df 9), p 4.5e-11 | 3.8e-09 |
| cl4 | chain | 329 | 209 | B | 0.079 | 0.072 | cluster-robust Wald 20.7 (df 6), p 0.0021 | 0.016 |
| cl4 | chain | 329 | 209 | Q+X | 0.208 | 0.200 | cluster-robust Wald 80.3 (df 12), p 3.5e-12 | 6e-10 |
| cl4 | entry | 128 | 55 | Q | 0.130 | 0.139 | LR 24.4 (df 3), p 2.1e-05 | 0.00051 |
| cl4 | entry | 128 | 55 | X | 0.207 | 0.157 | LR 27.4 (df 9), p 0.0012 | 0.011 |
| cl4 | entry | 128 | 55 | B | 0.117 | 0.078 | LR 13.7 (df 6), p 0.033 | 0.13 |
| cl4 | entry | 128 | 55 | Q+X | 0.346 | 0.307 | LR 53.7 (df 12), p 3.1e-07 | 1.3e-05 |

BH: one family over every p-value of the step (169 tests).

Sensitivity outcome (fraction observed), gain entered first, Q+X vs B: cl7 chain 0.272 vs 0.352; cl4 chain 0.067 vs 0.057; cl7 entry 0.676 vs 0.405; cl4 entry 0.209 vs 0.118.

Odds ratios of the full model with BH q < 0.05 (complete outcome): cl4 chain ESR1 chains in ASU (X, per 1 SD) OR 0.333 [0.209, 0.53], q 0.00012; cl4 chain space group P21 (X, present vs absent) OR 78.9 [8.97, 695], q 0.0014; cl4 chain other mutations (B, present vs absent) OR 0.0827 [0.0187, 0.366], q 0.011; cl4 chain R-free (Q, per 1 SD) OR 2.2 [1.35, 3.57], q 0.014; cl4 chain resolution (A) (Q, per 1 SD) OR 0.348 [0.177, 0.684], q 0.016; cl4 chain ligand heavy atoms (B, per 1 SD) OR 0.514 [0.329, 0.802], q 0.023; cl4 entry ESR1 chains in ASU (X, per 1 SD) OR 0.104 [0.0343, 0.318], q 0.0014; cl4 entry resolution (A) (Q, per 1 SD) OR 0.139 [0.0427, 0.449], q 0.011; cl7 chain ligand heavy atoms (B, per 1 SD) OR 0.288 [0.131, 0.636], q 0.016; cl7 chain space group P21 (X, present vs absent) OR 280 [6.36, 1.23e+04], q 0.023.

Dropped predictors (complete outcome), cluster/level: cl7 chain: ambient_temp_K_missing [X] separation (with = 1: 4 chains, all complete), partner_peptide [B] collinear with earlier predictors, mut_L536S [B] constant, mut_D538G [B] separation (with = 1: 4 chains, all complete), basic_amine [B] joint separation (no convergence;|log OR|>15; |log OR| 82.8); cl7 entry: ambient_temp_K_missing [X] separation (with = 1: 2 chains, all complete), partner_peptide [B] collinear with earlier predictors, mut_L536S [B] constant, mut_D538G [B] separation (with = 1: 2 chains, all complete), basic_amine [B] joint separation (no convergence;|log OR|>15; |log OR| 161.8); cl4 chain: peptide [B] separation (with = 1: 4 chains, none complete), partner_peptide [B] separation (with = 1: 4 chains, none complete), mut_D538G [B] separation (with = 1: 8 chains, none complete); cl4 entry: peptide [B] separation (with = 1: 3 chains, none complete), partner_peptide [B] separation (with = 1: 3 chains, none complete), mut_D538G [B] separation (with = 1: 1 chains, none complete).

## Verdict (verdict.md)

# SV01a verdict: **mixed** (H12 complete, window 538-547; [DEFAULT] rule, chain and entry level)
- Chain level (mixed): McFadden gain entered first, cl7 Q+X 0.224 vs B 0.225; cl4 Q+X 0.208 vs B 0.079.
- Entry level (artifact): cl7 Q+X 0.613 vs B 0.315; cl4 Q+X 0.346 vs B 0.117.
- With the same rows in every model, the 'entered last' comparison of Q+X vs B equals the 'entered first' one (both reduce to R2(Q+X) vs R2(B)); values for both orders are in model_blocks.csv.
- Sensitivity, fraction of 538-547 observed: **mixed** (chain mixed, entry artifact).
- Events vs fitted predictors (complete outcome): cl7 chain 14 incomplete / 16 predictors (full R2 0.378); cl7 entry 5 incomplete / 16 predictors (full R2 0.755); cl4 chain 209 incomplete / 18 predictors (full R2 0.280); cl4 entry 55 incomplete / 18 predictors (full R2 0.424).
- Dropped at chain level (complete outcome; full list per level in odds_ratios.csv): cl7 ambient_temp_K_missing [X, separation], partner_peptide [B, collinear with earlier predictors], mut_L536S [B, constant], mut_D538G [B, separation], basic_amine [B, joint separation]; cl4 peptide [B, separation], partner_peptide [B, separation], mut_D538G [B, separation].
- First pass without lattice contacts and mobility; SV01b refits with them and re-issues the verdict.

Reading: the chain-level cl7 comparison is a near tie (B exceeds Q+X by 0.0017), cl7 has few incomplete chains for the number of predictors, and the drop rule removed several B predictors that separate the outcome (their information is not credited to B). SV01b refits with lattice contacts.

## Files

- results/sv/SV01a/chain_outcomes.csv (42.0 kB)
- results/sv/SV01a/completeness_profile.csv (230.0 kB)
- results/sv/SV01a/concordance.csv (18.7 kB)
- results/sv/SV01a/concordance_summary.csv (0.6 kB)
- results/sv/SV01a/env.txt (0.3 kB)
- results/sv/SV01a/h12_status_counts.csv (2.1 kB)
- results/sv/SV01a/model_blocks.csv (19.8 kB)
- results/sv/SV01a/model_data.csv (98.4 kB)
- results/sv/SV01a/numbers.json (0.3 kB)
- results/sv/SV01a/odds_ratios.csv (36.2 kB)
- results/sv/SV01a/partial_ends.csv (11.2 kB)
- results/sv/SV01a/predictor_scaling.csv (1.2 kB)
- results/sv/SV01a/verdict.md (1.3 kB)
- results/sv/SV01a/verdict_table.csv (1.6 kB)
- figures/sv/SV01a/SV01a_completeness_profile.png (121.7 kB)
- figures/sv/SV01a/SV01a_completeness_profile.svg (17.4 kB)
- figures/sv/SV01a/SV01a_forest_completeness.png (328.8 kB)
- figures/sv/SV01a/SV01a_forest_completeness.svg (152.5 kB)
- figures/panels/sv/SV01a/SV01a_completeness_profile__profile.png (121.3 kB)
- figures/panels/sv/SV01a/SV01a_forest_completeness__forest_cl4.png (143.8 kB)
- figures/panels/sv/SV01a/SV01a_forest_completeness__forest_cl7.png (146.7 kB)
- figures/panels/sv/SV01a/SV01a_forest_completeness__gains_cl4.png (53.7 kB)
- figures/panels/sv/SV01a/SV01a_forest_completeness__gains_cl7.png (51.7 kB)
- scripts/sv/SV01a_completeness.py (31.3 kB)
- scripts/sv/SV01a_figures.py (9.8 kB)
- scripts/sv/SV01a_summary.py (9.3 kB)
- scripts/sv/common/sv_fig.py (2.1 kB)
- results/sv/SV01a/SUMMARY.md (this file)
