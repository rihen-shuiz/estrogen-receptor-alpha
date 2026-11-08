# Limitations added by RUN_SV1028 (supervisor requests)

Added by step SV09 of `docs/sv/RUN_SV1028.md`. The numbers come from `results/sv/<ID>/` and are checked in
`results/sv/SV09/number_check.csv`; the results themselves are summarised at the end of `docs/RESULTS_SUMMARY.md`.

## Secondary structure from pydssp, not DSSP

mkdssp could not be installed, so SV02 and everything that uses its helix table (helix ends and lengths in SV03, helix
completeness in SV04, candidate positions in SV05, secondary-structure context in SV07) rests on pydssp. Its 'H' covers
every helix type, it places the amide hydrogen differently from DSSP, and it has no chain-break handling: chains were
run in segments, so beta bridges across breaks are missed. Helix ends next to unobserved residues are censored; the
H12 end is censored in ⟦223‖S02_nc‖per_chain_status/2/end_cens⟧ of the ⟦761‖S02_nc‖per_chain_status/2/n_helix⟧ chains with an H12 run.

## The H12 window is not signed off

Primary analyses use H12 ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧; the data-derived window of SV02 is ⟦538‖S01c_n‖window_proposed/0⟧-⟦544‖S01c_n‖window_proposed/1⟧. The completeness verdict depends on
the choice: artifact on ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ (SV01b), mixed on ⟦538‖S01c_n‖window_proposed/0⟧-⟦544‖S01c_n‖window_proposed/1⟧ (SV01c), where an incomplete H12 becomes rare (cl4 ⟦209‖S01c_n‖n_incomplete_chains/538_547/cl4⟧ -> ⟦19‖S01c_n‖n_incomplete_chains/538_544/cl4⟧
chains) and the entry-level verdict rests on ⟦4‖S01c_n‖n_incomplete_entries/538_544/cl7⟧ (cl7) and ⟦6‖S01c_n‖n_incomplete_entries/538_544/cl4⟧ (cl4) incomplete entries. Until the team signs off a
window, statements about H12 completeness hold for ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧ only.

## Lattice contacts and completeness are partly circular

SV01b counts lattice contacts on the modelled atoms of residues ⟦528‖S01c_n‖lattice_range/0⟧-⟦550‖S01c_n‖lattice_range/1⟧, so a fully modelled H12 has more atoms
that can touch a neighbour, and the completeness verdict 'artifact' partly restates that modelled residues make
contacts. The medoid clash test that could separate cause from consequence ran only on the ⟦31‖S01c_cl‖n_chains‖cluster_group=all&measure=both positions blocked⟧ chains that are
incomplete on ⟦538‖S01c_n‖window_proposed/0⟧-⟦544‖S01c_n‖window_proposed/1⟧ (SV01c), not as a test on the ⟦235‖S01b_sum‖on the (\d+) partial/not-modelled chains⟧ partial or not-modelled chains of ⟦538‖S01c_n‖window_primary/0⟧-⟦547‖S01c_n‖window_primary/1⟧. The H12 mobility
outcome is z-scored on modelled residues and is not tied to missing atoms in this way.

## Few events in cl7

cl7 has ⟦14‖S01c_n‖n_incomplete_chains/538_547/cl7⟧ chains without a complete H12, in ⟦5‖S01c_n‖n_incomplete_entries/538_547/cl7⟧ entries, against ⟦17‖S01b_mb‖n_predictors_full‖outcome=complete&cluster=7&level=chain&block=Q&order=first⟧ fitted predictors in SV01b. Predictors
that separate the outcome are dropped, and their information is credited to no block. The cl7 block gains, and those
of the entry-level models without TLS-refined entries, are optimistic.

## Construct mutations and crystal recipes go with the clusters

⟦443‖S01a_h‖n‖level=chain&group=Y537S⟧ chains carry Y537S and ⟦243‖S01a_h‖n‖level=chain&group=L536S⟧ carry L536S; L536S does not occur in cl7, and the opened literature describes both as
mutations introduced to drive the conformation or to aid crystallisation (SV08). Effects of mutations, bound peptide
and ligand chemistry found within a cluster (SV01, SV03, SV07) are therefore associations inside particular crystal
systems, not causal effects. No apo wild-type ERalpha structure exists (SV08), so H12 disorder is never compared with an
unliganded baseline.

## Midpoint clustering is not independent evidence

The agreement of the SV04 midpoint clustering with S08 (ARI ⟦0.954‖S04_ari‖ari‖distance=midpoint_rms&level=chain&K=2⟧) rests on the same observed H12 C-alpha that S08 used,
and the H12 truncation curve rests on the ⟦75‖S04_ve‖n_complete_chains‖helix=H12&end=N&k_removed=0⟧ chains whose whole H12 is helical.

## One spike-audit rule could not be applied

The 'sparse-sample effect' rule of SV07 is degenerate for cluster means (BLOCKED part), so it cannot exclude a spike
caused by which chains observe a loop. The covariate-adjusted spike coefficients rest on the few chains whose covariate
value is rare in their cluster.

## Action classes are provisional and coincide with the clusters

The SV12 classes are ChEMBL mechanism labels for ⟦51‖S12_sum‖covers (\d+) of 362 rows⟧ of ⟦362‖S12_n‖n_rows⟧ entry x ligand rows (⟦14‖S12_sum‖\((\d+) of 306 ligands⟧ ligands); the DOI-based route failed its
validation gate and the teammate's labels were not available. Agonist entries lie almost all in cl7 and SERM and SERD
entries in cl4, so class and cluster cannot be separated, and every SV13 class comparison failed the coverage gate.
The SV13 SERM representative has no modelled H12, so its panel shows the ligand but not H12. Crystal structures show
H12 order and position; none of the SV steps measures receptor degradation or cellular activity, and the opened
literature reports the link from H12 destabilisation to degradation as asserted, not shown (SV08b).

## Pocket volumes rest on an in-house grid and an E2 reference

No pocket program was used. The SV14 grid (probe ⟦1.4‖S14_ref‖grid/probe_A⟧ A, buriedness by ray casting) has no published ERalpha
counterpart with a stated method (SV08b). The ⟦100‖S14_sum‖pocket (100) % = wild-type⟧ % pocket is the wild-type E2 pocket because only ⟦3‖S14_ref‖apo_baseline/entries_passing⟧ apo entries pass,
the E2 pocket itself is bimodal, the pocket reaches the grid box in ⟦203‖S14_n‖touches_box_b10_chains⟧ of ⟦788‖S14_n‖chains⟧ chains (those volumes are bounded by
the box), and the 'pocket open' flag is true for every chain, so it carries no information. Filling percentages hold
only for this method and reference.

## Affinities cover few ligands

Same-assay affinities relative to E2 (dpX) exist for ⟦14‖S15_n‖tested_ligands_primary⟧ ligands; ⟦48‖S15_n‖basis_counts_primary_set/absolute⟧ more have only absolute values from different
assays, which are drawn but not compared. pChEMBL values mix IC50, Ki and Kd (also in the LE and LipE of SV10). All are
in-vitro binding to the ERalpha target, not cellular activity.

## Pocket detection saw no partner chains

P2Rank (SV11) ran on files holding only the receptor chain and its ligand, so pockets formed together with a bound
peptide or a dimer partner are not seen as such, and the AF-2 groove and dimer-interface pockets it reports are
low-probability pockets.

## Literature evidence is uneven

SV08 and SV08b record claims only from opened text, read through a reader that caps verbatim quotes at ⟦125‖S08_sum‖quotes at (\d+) characters⟧ characters.
Some claims come from abstracts or preprints, the primary sources of the published ERalpha volumes could not be opened,
and no opened source gives the residue-level effect of L536S.
