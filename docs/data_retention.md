# Data retention: what was removed before, what is kept now

This note answers the professor's notes of 2026-10-18. Every number comes from a results file, which is named in brackets.

## In short

Nothing is removed for quality any more. Every ESR1 chain keeps its row and carries labels: the QC tier, a warning, an H12 B-factor group and its cluster. The main analysis now uses all 788 eligible chains from 351 entries, instead of 351 chains (one per entry). Only 113 of the 901 ESR1 chains are outside the RMSD matrix, and each has a mathematical reason: it has no LBD residues, or too few, to measure (results/chain_status.csv).

## Old strict set vs now

| | Old strict run | Now (keep-all) |
|---|---|---|
| Chains in the RMSD matrix | 351 (one per entry) | 788 (every eligible chain) |
| Entries | 351 | 351 |
| Alignment core | 84 residues (340-360, 400-420, 450-470, 500-520) | 123 residues, core v3 (see note 1) |
| Clusters at 1.0 / 1.5 / 2.0 / 2.5 Å | 18 / 9 / 5 / 4 | 16 / 7 / 4 / 3 |

Sources: results/archive_v2_strict/cluster_stats.csv, results/cluster_stats.csv.

## What each former filter removed (K01, results/attrition_funnel.csv)

Necessary steps (still applied, because the chain is not an ESR1 LBD or has nothing to measure):

- 149 chain-models from the 39 ESR2 validation entries. These are analysed separately as ESR2.
- 462 peptide chains, 146 partner chains and 8 DNA chains in ESR1 entries. These are not the receptor.
- 109 ESR1 chains with no LBD overlap (constructs without the LBD).

Project choices that removed chains before and do not now:

- 4 chains below the LBD/core C-alpha screen (at least 180 LBD C-alpha and 60 core C-alpha). The screen is unchanged, so these chains stay outside the matrix. Their rows are kept, with the reason.
- 437 chains dropped by "one chain per entry" (788 to 351). All 437 are now in the matrix.
- 8 entries dropped for R-free above 0.30, in the old sensitivity comparison only. These are now labelled, not dropped.
- Distance pruning removed 0 chains in both runs. The smallest pairwise overlap is 177 shared C-alpha, above the minimum of 150.

## QC labels instead of removal (K05, K06)

Of the 788 chains in the matrix, 544 are tier A, 240 are tier B and 4 are tier C. The 244 tier B or C chains carry a written warning, such as "resolution 2.70 A (tier B)". The other 544 have warning "none". Over all 901 ESR1 chains, 8 are labelled "unknown" because the entry has no R-free (results/chain_status.csv). The tier limits are in esr_analysis/resources/config/qc_tiers.json.

## Effect of the strict subsets (S13, R06)

The table compares the clustering when only some chains are reclustered. The subsets are sliced from the same saved matrix (results/qc_clustering_delta.csv).

| Subset | Chains | Clusters at 1.0 / 1.5 / 2.0 / 2.5 Å | ARI vs all chains |
|---|---|---|---|
| All chains | 788 | 16 / 7 / 4 / 3 | 1.0 |
| Tier A + B | 784 | 16 / 7 / 4 / 3 | 1.0 at every cut |
| Tier A only | 544 | 14 / 4 / 3 / 2 | 0.91 / 0.98 / 0.99 / 0.99 |

Removing tier C changes nothing. Keeping only tier A merges some small clusters at 1.5 Å, and 1 of the 7 all-chain clusters has no tier-A member. The large clusters stay the same. Entropy changes by at most 0.07 bits under either strict subset (results/entropy_qc_compare.csv).

## Entry-level vs chain-level view (K07, results/entry_vs_chain_clusters.csv)

Chains of one entry come from the same crystal, so they are not independent. The entry-level view takes one chain per entry: the chain with the most mapped LBD C-alpha, then the lowest mean B. These are 351 chains, sliced from the same matrix. This view gives 18 / 6 / 3 / 3 clusters at 1.0 / 1.5 / 2.0 / 2.5 Å, against 16 / 7 / 4 / 3 for all chains. The ARI between the two views is 0.84 / 0.97 / 0.98 / 0.99. Both views are kept.

## Note 1: which helices are most stable by B-factor group, and the core used

Each chain's C-alpha B-factors were divided by the chain's own median B. The chains were then split into low, mid and high tertiles of mean B (docs/helix_stability.md, results/helix_stability_by_bgroup.csv). In all three groups, the five most stable helices are AC1 441-456 (rank 1 everywhere), AB7 371-396 (rank 2 everywhere), AC3 496-531, AB5 338-362 and AB9 420-439. The least stable helices are AB8 412-416 and AB3 306-310. These five helices form the alignment core now used: core v3, 123 residues, without H12 (esr_analysis/resources/config/alignment_core_v3.json). The median H12/core B ratio is 1.255 over 765 chains, with tertile medians of 1.30, 1.22 and 1.25 (results/h12_core_ratio_by_bgroup.csv). H12 is therefore about 1.25 times the core in mean B, not 12 times.

Changing the core barely moves the primary result. Against the old 84-residue core, the pairwise RMSDs correlate at Spearman 0.988, and the 1.5 Å clustering agrees at ARI 0.96, with 7 clusters in both. At the finer 1.0 Å cut the core matters more: ARI 0.77, with 20 clusters on the old core against 16 (results/core_sensitivity.csv).

## Note 2: label a structure you do not like, keep as much data as possible

This is done. No chain is removed for resolution, R-free, B-factor, H12 disorder, a large core RMSD, or because it is an apo structure. Each of these is a label column in results/chain_status.csv:

- qc_label and warning;
- h12_label: 553 complete, 212 partial, 23 not modelled (results/h12_metric.csv);
- b_group and h12_b_group;
- core_rmsd_label.

H12 stays in the RMSD. The strict 351-chain outputs are archived in results/archive_v2_strict/ and are not deleted. The only chains outside the matrix are the 113 with no LBD residues, or too few, to measure. Each is listed in results/chain_status.csv with in_matrix false and its reason.
