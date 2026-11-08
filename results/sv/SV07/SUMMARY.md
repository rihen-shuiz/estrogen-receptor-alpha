Not done / substituted: (1) PARTLY BLOCKED - the 'sparse-sample effect (gone when restricted)' rule is degenerate for cluster-mean spikes: Y10 already averages, at each position, only the chains that observe it (nan-mean), so the cluster means recomputed on the chains observing the position in both clusters equal Y10 at all 246 positions (max |diff| 2e-14 A) and no cluster-mean spike can be 'gone when restricted'; the rule could still act through the medoid measure (complete-H12 medoids re-picked among the chains observing the position; 6 spikes exceed on the medoid). Tried: the other reading of 'in both clusters' (entries holding chains of both clusters) - only 2 of 347 cl7/cl4 entries do, too few to recompute a profile; not used. (2) The run file gives no operational rule for 'driven by <= 3 entries', 'numbering error', 'persists across many entries' or the covariate tests; declared in scripts/sv/SV07_spikes.py before any diagnostic was computed: greedy removal of up to 3 entries (each time the entry whose removal lowers d_mean most); numbering flag = residue differs from UniProt with no deposited explanation (residue_status unexplained_difference) or UniProt-author offset not constant over the LBD / not from SIFTS (SV00); covariate tests = OLS of the signed C-alpha coordinate along the cl7->cl4 displacement on cluster + covariate, entry-clustered Wald chi2 (chain) and nested F (entry x cluster units), one BH family over 84 tests, 'associated' = q < 0.05 at both levels; primary verdict at chain level (the Y10 unit), entry level reported. (3) SV01b saved B z only for 528-547; the SV01b formula was applied to all positions and checked against that file (13575 residues, max |diff| 9e-16). (4) SSE context is SV02's pydssp fallback ('H' = any helix). (5) Supporting files not named in the run file: top_entries.csv, restricted_profile.csv, covariate_tests.csv, lattice_contacts_spikes_per_residue.csv, lattice_contacts_spikes_per_chain.csv, lattice_checks.json, numbers.json.

# SV07 - Spike audit of the Y10 displacement profile (slide 5)

**Verdict:** 19 positions outside 528-550 meet the [DEFAULT] spike rule (Y10 cluster-mean or medoid-pair C-alpha distance > 1.0 A). **14 are genuine loop differences** (331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 464, 466, 467, 468): the cl7-cl4 distance stays > 1.0 A after removing the 3 most influential entries, at chain and entry level. **5 are driven by <= 3 entries**: 330, 397, 420, 459 exceed only on the medoid pair (cluster means 0.18-0.95 A), and 462 drops to 0.96 A without one entry (4TUZ). No spike is a numbering error (0 unexplained residue differences at the spike positions; removing the 11 chains with a non-constant or non-SIFTS offset changes d_mean by <= 0.08 A) and none is a sparse-sample effect under the rule as written (see line 1). Chain- and entry-level verdicts agree for 19 of 19 spikes. Among the genuine spikes, ligand chemotype is associated (q < 0.05 at chain and entry level) at 331, 332, 334, 335, 336, 337, 338, 339, dimer-partner source at 334, 335, 336, 337, 338, 339, 340, 467, 468 and peptide_bound at 331, 466, 467, 468; 333, 464 show no association. At 335 the cl4-cl7 coefficient falls from 2.74 to 0.51 A with chemotype in the model (mean signed coordinate cl4/acid: 1.12; cl4/amine: 3.75; cl4/both: 4.97; cl4/other: 0.44; cl7/acid: 0.36; cl7/amine: 1.43; cl7/none: -0.95; cl7/other: -0.00). At 467 it falls from 2.19 to 0.59 A with peptide_bound (mean signed coordinate cl4/False: 2.19; cl4/True: 2.06; cl7/False: 1.86; cl7/True: -0.31; chains cl7: False 22, True 132; cl4: False 240, True 3). These covariates are unevenly spread over the clusters, so the adjusted coefficients rest on the few chains whose value is rare in their cluster.

## Inputs (rule 1.10)

- Y10 card out results/ext/cluster_displacement_profile.csv (spike rule applied to it).
- Core v3 C-alpha: SV00 gate fit data/sv/gate/aligned_coordinates.npz = S03 results/alignment/aligned_coordinates.npz (max |diff| 4e-14 A, same gaps True); Y10 reproduced with SV05's y10_profile: max |diff| 2e-14, same NaN True, medoids 5TM1|A|1 / 7UJ7|C|1.
- SV00 observed_mask.npy, chain_table.csv (cluster, peptide_bound, ligand_chemotype, dimer_partner_source, offsets), chain_offsets.csv; residue identity data/sv/lbd_extract/residue_mapping.csv (181829 of 181829 observed C-alpha mapped).
- SV02 sse_per_position.csv, sse_codes_per_chain.csv, helix_boundaries.csv; SV01b B z formula; SV01b module lattice_contacts.py (10423 spike residues found = 10423 observed in the SV00 mask, 0 mismatches).

## Spikes (spikes.csv, verdict.csv)

d = cl7-cl4 C-alpha distance (A). Entry = entry-level means. -3 = after greedy removal of 3 entries. Medoid restr. = complete-H12 medoids re-picked among chains observing the position. SSE = SV02 consensus helix cl7 / cl4 and H share among observing chains. B z = median per cluster. Lattice = share of observing chains with a symmetry-copy contact; partner = ASU copy of another chain.

| pos | exceeds | obs cl7 / cl4 | d mean Y10 | d entry | d -3 (chain / entry) | d medoid Y10 / restr. | SSE cl7 / cl4 | B z cl7 / cl4 | lattice cl7 / cl4 | partner cl7 / cl4 | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 330 | medoid | 98% (436) / 95% (313) | 0.95 | 0.96 | 0.92 / 0.94 | 1.69 / 1.69 | loop 0.00 / loop 0.00 | 1.7 / 1.4 | 0.12 / 0.50 | 0.00 / 0.03 | driven by <= 3 entries |
| 331 | mean | 89% (395) / 77% (253) | 2.15 | 2.21 | 2.12 / 2.18 | - / 3.61 | loop 0.00 / loop 0.00 | 2.1 / 1.9 | 0.65 / 0.61 | 0.00 / 0.03 | genuine loop difference |
| 332 | mean | 78% (347) / 63% (206) | 2.53 | 2.68 | 2.45 / 2.60 | - / 2.56 | loop 0.00 / loop 0.00 | 2.9 / 2.5 | 0.44 / 0.37 | 0.00 / 0.00 | genuine loop difference |
| 333 | mean | 73% (324) / 57% (187) | 3.58 | 3.60 | 3.46 / 3.48 | - / 4.81 | loop 0.00 / loop 0.00 | 3.3 / 2.4 | 0.82 / 0.66 | 0.01 / 0.05 | genuine loop difference |
| 334 | mean | 72% (320) / 54% (178) | 3.36 | 3.42 | 3.19 / 3.27 | - / 6.11 | loop 0.00 / loop 0.00 | 3.6 / 2.5 | 0.71 / 0.82 | 0.01 / 0.05 | genuine loop difference |
| 335 | mean | 74% (330) / 55% (180) | 2.74 | 2.85 | 2.53 / 2.73 | - / 6.13 | loop 0.00 / loop 0.00 | 3.4 / 2.5 | 0.69 / 0.84 | 0.00 / 0.06 | genuine loop difference |
| 336 | mean | 81% (362) / 57% (188) | 2.78 | 2.80 | 2.58 / 2.70 | - / 4.98 | loop 0.00 / loop 0.00 | 3.2 / 2.2 | 0.80 / 0.80 | 0.01 / 0.04 | genuine loop difference |
| 337 | mean | 90% (400) / 57% (189) | 2.63 | 2.80 | 2.47 / 2.72 | - / 5.29 | loop 0.00 / loop 0.07 | 2.7 / 2.5 | 0.48 / 0.79 | 0.00 / 0.03 | genuine loop difference |
| 338 | mean | 95% (424) / 54% (177) | 2.28 | 2.59 | 2.14 / 2.45 | - / 4.33 | loop 0.00 / loop 0.08 | 2.5 / 2.2 | 0.54 / 0.59 | 0.00 / 0.01 | genuine loop difference |
| 339 | mean + medoid | 95% (424) / 54% (177) | 2.70 | 3.34 | 2.55 / 3.17 | 4.40 / 4.40 | H3 0.94 / loop 0.68 | 2.0 / 2.1 | 0.03 / 0.10 | 0.00 / 0.00 | genuine loop difference |
| 340 | mean + medoid | 98% (435) / 60% (199) | 1.91 | 2.37 | 1.75 / 2.25 | 3.31 / 3.31 | H3 0.94 / loop 0.59 | 1.7 / 1.9 | 0.16 / 0.06 | 0.00 / 0.00 | genuine loop difference |
| 397 | medoid | 100% (445) / 100% (329) | 0.18 | 0.23 | 0.18 / 0.22 | 1.03 / 1.03 | loop 0.05 / loop 0.06 | 1.3 / 1.3 | 0.06 / 0.16 | 0.00 / 0.00 | driven by <= 3 entries |
| 420 | medoid | 92% (411) / 84% (277) | 0.22 | 0.21 | 0.13 / 0.17 | 1.05 / 1.05 | loop 0.00 / loop 0.03 | 2.1 / 2.5 | 0.00 / 0.00 | 0.00 / 0.00 | driven by <= 3 entries |
| 459 | medoid | 97% (431) / 92% (302) | 0.41 | 0.40 | 0.39 / 0.38 | 1.04 / 1.04 | loop 0.23 / loop 0.45 | 2.3 / 1.5 | 0.02 / 0.11 | 0.97 / 0.87 | driven by <= 3 entries |
| 462 | mean | 12% (52) / 30% (100) | 1.04 | 1.14 | 0.83 / 0.86 | - / 1.49 | loop 0.00 / loop 0.00 | 4.7 / 3.1 | 0.15 / 0.16 | 0.19 / 0.18 | driven by <= 3 entries |
| 464 | mean | 15% (68) / 45% (148) | 1.52 | 1.70 | 1.16 / 1.29 | - / 4.26 | loop 0.00 / loop 0.00 | 4.6 / 2.9 | 0.63 / 0.62 | 0.01 / 0.16 | genuine loop difference |
| 466 | mean | 28% (126) / 67% (221) | 1.86 | 2.19 | 1.57 / 1.84 | - / 4.87 | loop 0.52 / H9 0.95 | 4.7 / 1.9 | 0.41 / 0.61 | 0.05 / 0.11 | genuine loop difference |
| 467 | mean | 35% (154) / 74% (243) | 2.19 | 2.50 | 2.02 / 2.29 | - / 6.09 | loop 0.53 / H9 0.89 | 4.6 / 1.8 | 0.25 / 0.49 | 0.06 / 0.09 | genuine loop difference |
| 468 | mean | 39% (175) / 80% (263) | 1.14 | 1.42 | 1.07 / 1.33 | - / 1.26 | loop 0.58 / H9 0.90 | 4.0 / 1.6 | 0.01 / 0.02 | 0.00 / 0.00 | genuine loop difference |

- Restricted recomputation: cluster means equal Y10 (max |diff| 2e-14 A). Re-picked medoids differ from Y10 at 1 position where Y10 has a value and give a value at 38 positions where the Y10 medoids do not observe the residue; there they exceed 1.0 A at 331, 332, 333, 334, 335, 336, 337, 338, 462, 463, 464, 465, 466, 467, 468, 470, 471 (not spikes: the step-1 rule uses the Y10 values).
- Entry-level cluster means exceed 1.0 A outside 528-550 at 331, 332, 333, 334, 335, 336, 337, 338, 339, 340, 462, 464, 465, 466, 467, 468 (max 3.60 A); 465 exceeds at entry level only (1.17 A vs Y10 0.95 A).
- Largest single-entry influence per spike (top_entries.csv; 5 per spike with offsets, residue identity, SSE code, B z, peptide, ligand): 330 7UJ7 +0.01 A; 331 8DUS +0.01 A; 332 8DUS +0.04 A; 333 8DUS +0.06 A; 334 8DUS +0.08 A; 335 7UJO +0.07 A; 336 7UJO +0.07 A; 337 7N9O +0.05 A; 338 6VPF +0.05 A; 339 7RS4 +0.06 A; 340 2I0J +0.06 A; 397 6ZOR +0.00 A; 420 7UJO +0.04 A; 459 8DUK +0.01 A; 462 4TUZ +0.08 A; 464 9EBG +0.19 A; 466 3Q95 +0.10 A; 467 2B1Z +0.06 A; 468 5DP0 +0.03 A.
- Numbering of the top entries (top_entries.csv): residue status at the spike position ['match'] for every listed chain; chains with a non-zero, non-constant or non-SIFTS offset among them: 2IOK|B|1 at 340 (offset -1000, constant True, sifts_pdbe_api); 9VT5|A|1 at 420 (offset 304, constant False, sifts_pdbe_api); 9VT5|B|1 at 420 (offset 0, constant False, sifts_pdbe_api). A constant non-zero offset (e.g. -1000 for author numbering 1305-1550) is a SIFTS-mapped constant, not a flag. Flagged chains overall (offset not constant or not SIFTS): 6IAR|A|1, 9VT5|A|1, 9VT5|B|1, 9VTU|A|1, 9VTU|B|1, 9W11|A|1, 9W11|B|1, 9W12|A|1, 9W12|B|1, 9W13|A|1, 9W13|B|1.

## Genuine spikes against peptide, ligand chemotype and dimer source (covariate_tests.csv)

q = BH over all covariate tests (chain: OLS, entry-clustered Wald chi2; entry: entry x cluster units, nested F). Coef = cl4-minus-cl7 difference of the signed coordinate along the displacement (A) without -> with the covariate.

| pos | peptide_bound q chain / entry | coef | ligand_chemotype q chain / entry | coef | dimer source q chain / entry | coef | associated (q < 0.05 both) |
|---|---|---|---|---|---|---|---|
| 331 | 0.028 / 0.00014 | 2.15 -> 1.24 | 0.015 / 0.0032 | 2.15 -> 2.14 | 0.24 / 0.18 | 2.15 -> 2.18 | ligand_chemotype, peptide_bound |
| 332 | 0.84 / 0.98 | 2.53 -> 2.44 | 0.00014 / 4.3e-09 | 2.53 -> 1.80 | 0.06 / 0.081 | 2.53 -> 2.48 | ligand_chemotype |
| 333 | 0.47 / 0.19 | 3.58 -> 3.25 | 0.01 / 0.23 | 3.58 -> 3.26 | 0.84 / 0.98 | 3.58 -> 3.59 | none |
| 334 | 0.7 / 0.31 | 3.36 -> 3.58 | 1.4e-14 / 2.7e-12 | 3.36 -> 2.01 | 1.2e-05 / 8.6e-05 | 3.36 -> 3.17 | dimer_partner_source, ligand_chemotype |
| 335 | 0.1 / 0.17 | 2.74 -> 3.43 | 6.3e-08 / 3.7e-19 | 2.74 -> 0.51 | 2.3e-07 / 4.3e-09 | 2.74 -> 2.38 | dimer_partner_source, ligand_chemotype |
| 336 | 0.096 / 0.38 | 2.78 -> 3.24 | 1.3e-07 / 8.1e-22 | 2.78 -> 0.66 | 1.8e-08 / 3e-11 | 2.78 -> 2.47 | dimer_partner_source, ligand_chemotype |
| 337 | 0.089 / 0.27 | 2.63 -> 2.99 | 2.5e-06 / 1.5e-19 | 2.63 -> 0.96 | 7.5e-08 / 3.6e-08 | 2.63 -> 2.39 | dimer_partner_source, ligand_chemotype |
| 338 | 0.11 / 0.11 | 2.28 -> 2.78 | 0.00036 / 8.8e-13 | 2.28 -> 1.03 | 0.00036 / 7.4e-07 | 2.28 -> 2.11 | dimer_partner_source, ligand_chemotype |
| 339 | 0.034 / 0.19 | 2.70 -> 3.36 | 0.0024 / 9.9e-09 | 2.70 -> 1.33 | 2.6e-06 / 1.2e-21 | 2.70 -> 2.34 | dimer_partner_source, ligand_chemotype |
| 340 | 0.87 / 0.84 | 1.91 -> 1.82 | 0.081 / 7e-07 | 1.91 -> 1.06 | 8.8e-06 / 4.8e-19 | 1.91 -> 1.67 | dimer_partner_source |
| 464 | 0.86 / 0.98 | 1.52 -> 1.64 | 0.84 / 0.84 | 1.52 -> 1.42 | 0.12 / 0.5 | 1.52 -> 1.50 | none |
| 466 | 0.0027 / 0.034 | 1.86 -> 0.90 | 0.66 / 0.98 | 1.86 -> 1.93 | 0.02 / 0.27 | 1.86 -> 1.82 | peptide_bound |
| 467 | 3.4e-09 / 0.00032 | 2.19 -> 0.59 | 0.86 / 0.99 | 2.19 -> 2.26 | 8.9e-06 / 0.038 | 2.19 -> 2.13 | dimer_partner_source, peptide_bound |
| 468 | 1.1e-08 / 0.00036 | 1.14 -> -0.04 | 0.7 / 0.98 | 1.14 -> 1.20 | 5.2e-05 / 0.041 | 1.14 -> 1.10 | dimer_partner_source, peptide_bound |


## Figure

- figures/sv/SV07/SV07_spike_profile.png/.svg: (a) Y10 profile, all chains and chains observing the position in both clusters (means coincide; medoids re-picked), 1.0 A threshold, spikes labelled with their verdict; (b) share of cl7 / cl4 chains observing each position. Panels in figures/panels/sv/SV07/.

## Files

- results/sv/SV07/spikes.csv (18.6 kB)
- results/sv/SV07/verdict.csv (6.7 kB)
- results/sv/SV07/top_entries.csv (29.5 kB)
- results/sv/SV07/restricted_profile.csv (38.2 kB)
- results/sv/SV07/covariate_tests.csv (26.6 kB)
- results/sv/SV07/lattice_contacts_spikes_per_residue.csv (528.4 kB)
- results/sv/SV07/lattice_contacts_spikes_per_chain.csv (32.6 kB)
- results/sv/SV07/lattice_checks.json (0.4 kB)
- results/sv/SV07/numbers.json (3.8 kB)
- results/sv/SV07/env.txt (0.5 kB)
- figures/sv/SV07/SV07_spike_profile.png (513.4 kB)
- figures/sv/SV07/SV07_spike_profile.svg (76.3 kB)
- figures/panels/sv/SV07/SV07_spike_profile__profile.png (378.1 kB)
- figures/panels/sv/SV07/SV07_spike_profile__observation.png (102.6 kB)
- scripts/sv/SV07_spikes.py (32.6 kB)
- scripts/sv/SV07_lattice.py (3.4 kB)
- scripts/sv/SV07_figures.py (8.3 kB)
- scripts/sv/SV07_summary.py (16.4 kB)
