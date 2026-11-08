# SV01a verdict: **mixed** (H12 complete, window 538-547; [DEFAULT] rule, chain and entry level)
- Chain level (mixed): McFadden gain entered first, cl7 Q+X 0.224 vs B 0.225; cl4 Q+X 0.208 vs B 0.079.
- Entry level (artifact): cl7 Q+X 0.613 vs B 0.315; cl4 Q+X 0.346 vs B 0.117.
- With the same rows in every model, the 'entered last' comparison of Q+X vs B equals the 'entered first' one (both reduce to R2(Q+X) vs R2(B)); values for both orders are in model_blocks.csv.
- Sensitivity, fraction of 538-547 observed: **mixed** (chain mixed, entry artifact).
- Events vs fitted predictors (complete outcome): cl7 chain 14 incomplete / 16 predictors (full R2 0.378); cl7 entry 5 incomplete / 16 predictors (full R2 0.755); cl4 chain 209 incomplete / 18 predictors (full R2 0.280); cl4 entry 55 incomplete / 18 predictors (full R2 0.424).
- Dropped at chain level (complete outcome; full list per level in odds_ratios.csv): cl7 ambient_temp_K_missing [X, separation], partner_peptide [B, collinear with earlier predictors], mut_L536S [B, constant], mut_D538G [B, separation], basic_amine [B, joint separation]; cl4 peptide [B, separation], partner_peptide [B, separation], mut_D538G [B, separation].
- First pass without lattice contacts and mobility; SV01b refits with them and re-issues the verdict.
