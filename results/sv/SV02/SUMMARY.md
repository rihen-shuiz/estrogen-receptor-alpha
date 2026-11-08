Not done / substituted: (1) mkdssp >= 4 is not installed (no conda/Homebrew) - the step's fallback pydssp 0.9.1 (numpy backend) ran instead; pydssp gives H / E / '-' only and its 'H' is the union of 3-, 4- and 5-turn helices, so the 'alpha (H)' fractions are 'any helix' and 3_10 (G) and pi (I) fractions are 'not available'; pydssp also places the amide H on the C(i-1)-N / CA-N bisector, not along the previous C=O as DSSP does. (2) pydssp has no chain-break handling: each chain ran as 2380 segments in total (consecutive label_seq_id, all of N/CA/C/O, C-N <= 2.5 A = mkdssp's peptide-bond limit), prolines as non-donors, 24 segments shorter than 6 residues padded with far-away dummies (pydssp raises below 6); E bridges between segments are therefore not detected (H and 'other' unaffected). (3) Step 7 mkdssp phi/psi comparison on 5 chains: not possible (mkdssp did not run); substitute sanity check only: gemmi.calculate_phi_psi on the restored mmCIFs, 5 chains, 1156 residues, max |dphi| 2.8e-14 and |dpsi| 2.8e-14 deg where residues i-1..i+1 have no altloc; 55 residues with an altloc differ by up to 8.2 / 11.2 deg (gemmi reads the first altloc, SV00 kept the highest-occupancy one). (4) Helix names: the helix table (K02, results/helix_table.csv) names only H12, so names come from the project's mapping of its rows, scripts/ext_common.HELIX_NAMES (H1, H3, H4/5, H6, H7, H8, H9, H10/11, H12) - H4/5 and H10/11 stay combined as there; no H2 and no H6 consensus helix (414-415 = 2 positions). (5) Implementation choices fixed before the output: low-support positions are trimmed from run ends (never a boundary, allowed inside; trimming occurred 0 times); per-chain censoring tests the two flanking residues - a one-residue test is False by construction under DSSP (the residue next to an H run closes its H-bond, so it is always observed), as a first run showed (0 of 6277) - this equals SV03's 'real end' rule; rama overlap at entry level uses entry x cluster units and entry-level cluster means. (6) Supporting files not named in the run file: sse_codes_per_chain.csv, consensus_runs_all.csv, dssp_run.json, numbers_consensus.json, numbers_rama.json, phipsi_check_gemmi.csv, data/sv/SV02/. The 537 psi check passed (no new values).

# SV02 - Consensus secondary structure and per-residue Ramachandran plots

**Verdict:** 8 global consensus helices ([DEFAULT] >= 4 positions with pydssp H in >= 80 % of observing chains). **Data-derived H12 window = 538-544** (proposal only; primary window stays 538-547); it equals none of 538-547, 537-546, 536-545 (Jaccard 0.70 / 0.70 / 0.70). Consensus H12 is 538-545 in cl7 and 537-544 in cl4. Ramachandran at 537: psi circular mean cl7 +156.6, cl4 -46.8 deg (EXPLORATION.md +157 / -47: within 10 deg, check passed); cl7-cl4 separation 157.2 deg, overlap 0 at chain and entry level.

## Steps 1-2 - DSSP (pydssp fallback) and numbering

- 788 kept chains, 182293 residues with all four backbone atoms passed (full chain, incl. residues outside 305-550; 9 residues lack an atom), 2380 segments (median 3 per chain).
- UniProt numbering: SV00 esr1_position = auth_seq_id + SV00 uniprot_author_offset for 728475 / 728475 mapped atoms.
- Positions 305-550 with a pydssp state: 181824 (SV00 C-alpha mask 181829; 5 C-alpha-observed residues lack another backbone atom). States: H 124973, - 50879, E 5972.
- Observing a position = the residue has N, CA, C and O (a pydssp state).

## Steps 3-4 - consensus helices (helix_boundaries.csv)

| helix | helix table | global | cl7 | cl4 | global - table (start, end) |
|---|---|---|---|---|---|
| H1 | AB4 311-322 | 312-321 | 312-321 | 312-321 | +1, -1 |
| H3 | AB5 338-362 | 342-361 | 339-361 | 342-361 | +4, -1 |
| H4/5 | AB7 371-396 | 372-394 | 372-394 | 372-394 | +1, -2 |
| H7 | AB9 420-439 | 422-437 | 422-437 | 422-437 | +2, -2 |
| H8 | AC1 441-456 | 442-455 | 442-455 | 442-455 | +1, -1 |
| H9 | AC2 465-493 | 473-492 | 473-492 | 466-492 | +8, -1 |
| H10/11 | AC3 496-531 | 497-526 | 497-530 | 497-525 | +1, -5 |
| H12 | AC4 537-546 | 538-544 | 538-545 | 537-544 | +1, -2 |

- Runs that failed the 4-position rule: global 339-340; global 414-415; global 471-471; cl7 414-415; cl4 414-415.
- Low-support positions (< 50 % of the group's chains observe them): global [462, 463, 464, 465, 466, 549, 550]; cl7 [462, 463, 464, 465, 466, 467, 468, 469, 549, 550]; cl4 [305, 462, 463, 464, 529, 530, 531, 532, 547, 548, 549, 550].
- Global helices: min H fraction 0.847 (H10/11), low-support positions inside helices 0.

## Step 5 - per-chain table (sse_per_chain.csv, chain x global consensus helix)

| helix | chains with an H run | not observed | observed, not helical | start censored | end censored |
|---|---|---|---|---|---|
| H1 312-321 | 788 | 0 | 0 | 4 | 0 |
| H3 342-361 | 788 | 0 | 0 | 189 | 0 |
| H4/5 372-394 | 788 | 0 | 0 | 2 | 1 |
| H7 422-437 | 788 | 0 | 0 | 116 | 0 |
| H8 442-455 | 788 | 0 | 0 | 0 | 2 |
| H9 473-492 | 788 | 0 | 0 | 434 | 0 |
| H10/11 497-526 | 788 | 0 | 0 | 0 | 172 |
| H12 538-544 | 761 | 23 | 4 | 55 | 223 |

Start/end = first/last position of the chain's H run(s) overlapping the consensus range; censored = residue start-1 or start-2 (end+1 or end+2) not observed or across a chain break.

## Step 6 - data-derived H12 window (h12_window.md)

- Proposal: 538-544 (7 positions). Not adopted: 538-547 stays primary; a human signs off. It differs from 538-547, so SV01c runs (SV01c step 9).

| window | shared | only in window | only in proposal | Jaccard | start, end shift |
|---|---|---|---|---|---|
| 538-547 (settings h12_range, primary) | 7 | [545, 546, 547] | - | 0.70 | +0, -3 |
| 537-546 (helix table AC4, 7NFB struct_conf) | 7 | [537, 545, 546] | - | 0.70 | +1, -2 |
| 536-545 | 7 | [536, 537, 545] | - | 0.70 | +2, -1 |

## Steps 7-11 - phi/psi and Ramachandran plots (rama_summary.csv)

- phipsi_per_chain.csv: 181827 chain x position rows (788 chains; residues with N, CA, C), phi in 179576, psi in 179510, both in 177265; bonded = C-N <= 2.0 A (2 of 179629 observed consecutive links are longer).
- Positions ([DEFAULT], fixed in SV02_phipsi.py before plotting): [537, 533, 534, 535, 536, 538, 539, 540, 329, 359, 457, 333, 467] (Y12 top_positions [537, 329, 359, 457]).

| position | context (helix table) | chains cl7 / cl4 / minor | phi, psi cl7 | phi, psi cl4 | separation | overlap chain | overlap entry |
|---|---|---|---|---|---|---|---|
| 537 | H12 (helix table AC4 537-546) | 429 / 304 / 11 | -68.7, 156.6 | -54.5, -46.8 | 157.2 | 0.00 | 0.00 |
| 533 | loop H10/11-H12 (helix table) | 351 / 144 / 10 | -81.3, -41.5 | -91.6, -23.0 | 21.2 | 0.30 | 0.26 |
| 534 | loop H10/11-H12 (helix table) | 360 / 211 / 11 | -121.7, 134.4 | -84.9, 132.1 | 36.9 | 0.29 | 0.23 |
| 535 | loop H10/11-H12 (helix table) | 382 / 254 / 11 | -70.7, 130.6 | -64.1, 129.6 | 6.7 | 0.34 | 0.31 |
| 536 | loop H10/11-H12 (helix table) | 408 / 284 / 11 | -99.0, 139.5 | -72.6, 168.0 | 38.8 | 0.09 | 0.05 |
| 538 | H12 (helix table AC4 537-546) | 435 / 311 / 12 | -60.5, -39.3 | -64.7, -40.0 | 4.3 | 0.32 | 0.32 |
| 539 | H12 (helix table AC4 537-546) | 435 / 313 / 12 | -65.1, -46.1 | -62.7, -47.0 | 2.5 | 0.43 | 0.44 |
| 540 | H12 (helix table AC4 537-546) | 435 / 313 / 12 | -63.4, -38.2 | -62.6, -41.5 | 3.4 | 0.37 | 0.32 |
| 329 | loop H1-H3 (helix table) | 436 / 313 / 14 | -88.4, 167.3 | -73.0, 151.9 | 21.8 | 0.17 | 0.12 |
| 359 | H3 (helix table AB5 338-362) | 444 / 329 / 14 | -72.1, -31.0 | -65.0, -38.1 | 10.0 | 0.26 | 0.21 |
| 457 | loop H8-H9 (helix table) | 438 / 313 / 14 | -112.3, 9.8 | -106.8, 3.5 | 8.3 | 0.35 | 0.37 |
| 333 | loop H1-H3 (helix table) | 306 / 178 / 9 | -70.1, -2.3 | -72.0, -6.9 | 5.0 | 0.47 | 0.49 |
| 467 | H9 (helix table AC2 465-493) | 126 / 219 / 12 | -68.7, -28.1 | -60.8, -43.5 | 17.3 | 0.21 | 0.26 |

Circular means in degrees; circular SDs, entry counts and per-cluster overlaps are in rama_summary.csv. Overlap = share of cl7 + cl4 units closer to the other cluster's mean (wrapped phi/psi distance).
- Step 11 check (537, circular-mean psi): cl7 +156.63 vs +157 (|d| 0.37), cl4 -46.82 vs -47 (|d| 0.18): pass.
- Step 12 (resume) not applicable: SV02 ran in full.

## Files

- results/sv/SV02/consensus_cl4.csv (0.7 kB)
- results/sv/SV02/consensus_cl7.csv (0.7 kB)
- results/sv/SV02/consensus_global.csv (0.7 kB)
- results/sv/SV02/consensus_runs_all.csv (1.3 kB)
- results/sv/SV02/helix_boundaries.csv (0.6 kB)
- results/sv/SV02/phipsi_check_gemmi.csv (111.1 kB)
- results/sv/SV02/phipsi_per_chain.csv (8.65 MB)
- results/sv/SV02/rama_summary.csv (3.9 kB)
- results/sv/SV02/sse_codes_per_chain.csv (383.8 kB)
- results/sv/SV02/sse_per_chain.csv (430.5 kB)
- results/sv/SV02/sse_per_position.csv (70.2 kB)
- results/sv/SV02/h12_window.md (3.0 kB)
- results/sv/SV02/dssp_run.json (0.7 kB)
- results/sv/SV02/numbers_consensus.json (2.9 kB)
- results/sv/SV02/numbers_rama.json (1.2 kB)
- results/sv/SV02/env.txt (0.5 kB)
- data/sv/SV02/dssp_residues.csv.gz (703.5 kB)
- data/sv/SV02/dssp_segments.csv (93.6 kB)
- figures/sv/SV02/SV02_helix_track.png (162.9 kB)
- figures/sv/SV02/SV02_helix_track.svg (61.2 kB)
- figures/sv/SV02/SV02_sse_heatmap.png (100.6 kB)
- figures/sv/SV02/SV02_sse_heatmap.svg (78.1 kB)
- figures/sv/SV02/rama/SV02_rama_329.png (66.2 kB)
- figures/sv/SV02/rama/SV02_rama_329.svg (96.5 kB)
- figures/sv/SV02/rama/SV02_rama_333.png (68.7 kB)
- figures/sv/SV02/rama/SV02_rama_333.svg (67.7 kB)
- figures/sv/SV02/rama/SV02_rama_359.png (52.5 kB)
- figures/sv/SV02/rama/SV02_rama_359.svg (99.8 kB)
- figures/sv/SV02/rama/SV02_rama_457.png (62.8 kB)
- figures/sv/SV02/rama/SV02_rama_457.svg (97.1 kB)
- figures/sv/SV02/rama/SV02_rama_467.png (55.8 kB)
- figures/sv/SV02/rama/SV02_rama_467.svg (53.0 kB)
- figures/sv/SV02/rama/SV02_rama_533.png (72.1 kB)
- figures/sv/SV02/rama/SV02_rama_533.svg (69.0 kB)
- figures/sv/SV02/rama/SV02_rama_534.png (72.8 kB)
- figures/sv/SV02/rama/SV02_rama_534.svg (76.6 kB)
- figures/sv/SV02/rama/SV02_rama_535.png (65.5 kB)
- figures/sv/SV02/rama/SV02_rama_535.svg (83.9 kB)
- figures/sv/SV02/rama/SV02_rama_536.png (69.4 kB)
- figures/sv/SV02/rama/SV02_rama_536.svg (89.9 kB)
- figures/sv/SV02/rama/SV02_rama_537.png (57.9 kB)
- figures/sv/SV02/rama/SV02_rama_537.svg (94.6 kB)
- figures/sv/SV02/rama/SV02_rama_538.png (55.7 kB)
- figures/sv/SV02/rama/SV02_rama_538.svg (96.6 kB)
- figures/sv/SV02/rama/SV02_rama_539.png (51.2 kB)
- figures/sv/SV02/rama/SV02_rama_539.svg (96.8 kB)
- figures/sv/SV02/rama/SV02_rama_540.png (50.5 kB)
- figures/sv/SV02/rama/SV02_rama_540.svg (96.8 kB)
- figures/panels/sv/SV02/SV02_helix_track__helix_fraction.png (111.1 kB)
- figures/panels/sv/SV02/SV02_helix_track__observing.png (49.9 kB)
- figures/panels/sv/SV02/SV02_rama_329__rama.png (65.9 kB)
- figures/panels/sv/SV02/SV02_rama_333__rama.png (68.4 kB)
- figures/panels/sv/SV02/SV02_rama_359__rama.png (52.3 kB)
- figures/panels/sv/SV02/SV02_rama_457__rama.png (62.5 kB)
- figures/panels/sv/SV02/SV02_rama_467__rama.png (55.6 kB)
- figures/panels/sv/SV02/SV02_rama_533__rama.png (71.8 kB)
- figures/panels/sv/SV02/SV02_rama_534__rama.png (72.5 kB)
- figures/panels/sv/SV02/SV02_rama_535__rama.png (65.2 kB)
- figures/panels/sv/SV02/SV02_rama_536__rama.png (69.1 kB)
- figures/panels/sv/SV02/SV02_rama_537__rama.png (57.7 kB)
- figures/panels/sv/SV02/SV02_rama_538__rama.png (55.4 kB)
- figures/panels/sv/SV02/SV02_rama_539__rama.png (51.0 kB)
- figures/panels/sv/SV02/SV02_rama_540__rama.png (50.3 kB)
- figures/panels/sv/SV02/SV02_sse_heatmap__heatmap.png (99.9 kB)
- scripts/sv/SV02_consensus.py (17.0 kB)
- scripts/sv/SV02_dssp.py (8.2 kB)
- scripts/sv/SV02_figures.py (10.6 kB)
- scripts/sv/SV02_phipsi.py (13.7 kB)
- scripts/sv/SV02_summary.py (12.1 kB)
- results/sv/SV02/SUMMARY.md (this file)
