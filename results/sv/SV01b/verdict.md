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
