# Results summary: ESR1 structural diversity and redundancy

What the project did and found. Every number below is rendered from result files by
`scripts/collect_numbers_ext.py`; each finding names its figure (`figures/...`) and its result file (`results/...`).
All figures, with captions and sources: `docs/gallery.html` and `docs/figures/INDEX.md`. Methods: `docs/METHODOLOGY.md`.

## The data

- {{R.n_chains_all_esr1}} ESR1 chains from the PDB; {{R.n_chains}} chains from {{R.n_entries}} entries pass the
  completeness screen and are all used (keep-all; {{R.tier_A_chains}} tier A, {{R.tier_B_chains}} tier B,
  {{R.tier_C_chains}} tier C). All {{Y02.n_xray}} entries are X-ray structures; {{R.n_peg_entries}} were grown with PEG;
  {{R.n_peptide_entries}} carry a co-crystallised peptide. Figure `figures/ext/Y02_metadata_overview.png`; file
  `results/ext/entry_metadata.csv`.
- Construct mutations are common: {{Y02.chains_mutation_status.Y537S}} chains carry Y537S and
  {{Y02.chains_mutation_status.L536S}} carry L536S; only {{Y02.chains_mutation_status.wild type}} are wild type in the
  LBD (`results/ext/chain_metadata.csv`).

## 1. There are two conformations, and many methods agree

- Average linkage at {{R.set_primary_structural_threshold_A}} A gives {{R.nclusters_cut1p5}} clusters: cl7
  ({{R.cl7_size}} chains), cl4 ({{R.cl4_size}}) and five minor clusters ({{R.minor_size}} chains). The median RMSD
  between cl7 and cl4 is {{R.median_cl7_cl4_A}} A (`results/entry_clusters.csv`; `figures/ext/Y07_rmsd_heatmap_annotated.png`).
- All five embeddings separate the two (`figures/ext/Y03_embeddings.png`).
- Clustering zoo: {{R.zoo_n_runs}} runs in {{R.zoo_n_families}} families. The declared rule selects {{Y04.best_run}}
  (silhouette {{Y04.best_sil|.2f}}, bootstrap stability {{Y04.best_stab|.2f}}, ARI with S08 {{Y04.best_ari_s08|.2f}}).
  The best runs of the top-7 families agree with S08 at ARI {{Y06.ari_top7_min|.2f}}-{{Y06.ari_top7_max|.2f}};
  density-based, medoid-based, single-linkage and affinity-propagation methods do not
  (`figures/ext/Y04_zoo_metrics.png`, `Y04_zoo_grid.png`, `Y06_agreement.png`, `Y06_alluvial.png`; `results/ext/zoo_runs.csv`).
- Compactness: cl7 mean pairwise RMSD {{Y05.cl7.mean_intra_A|.2f}} A (diameter {{Y05.cl7.diameter_A|.1f}} A), cl4
  {{Y05.cl4.mean_intra_A|.2f}} A (diameter {{Y05.cl4.diameter_A|.1f}} A); cl4 is looser (`figures/ext/Y05_compactness.png`).
- Robustness: on core v3 and on the old core, every method in every subset (all, tier A, R-free <= 0.25, one chain
  per entry, complete H12) recovers the split with agreement {{Y28.real_min|.2f}}-{{Y28.real_max|.2f}}
  (`figures/ext/Y28_robustness_grid.png`).
- Caveat found on the way (`figures/ext/Y42_metric_check.png`): the RMSD over pairwise-common positions is not a metric
  ({{Y42.full.triangle_violation_frac|.1%}} of triples break the triangle inequality); chains without H12 look close
  to everyone. This is why PAM picked two H12-less chains as medoids.

## 2. Why two clusters: two positions of helix 12

- On H12 alone the median cl7-cl4 distance is {{Y42.h12only.median_cl7_cl4_A|.1f}} A; without residues 528-550 it is
  {{Y42.noh12.median_cl7_cl4_A|.2f}} A and the clusters merge (`figures/ext/Y42_region_rmsd.png`). Positions 528-550 carry
  {{Y10.share_528_550|.0%}} of the squared difference between the cluster means (`figures/ext/Y10_displacement_profile.png`).
- Geometry (`figures/ext/Y11_features_by_cluster.png`): H12 sits {{R.h12_shift_cl4}} A from its cl7 position in cl4
  (cl7 {{R.h12_shift_cl7}} A); the H12-H3 axis angle is {{R.h12_h3_angle_cl7|.0f}} vs {{R.h12_h3_angle_cl4|.0f}} deg; the
  K362-E542 C-alpha distance is {{R.k362_e542_cl7}} vs {{R.k362_e542_cl4}} A. The two distributions do not overlap.
  This is the agonist (H12 capping the pocket, coactivator groove free) vs antagonist (H12 in the groove) geometry.
- Without any superposition, backbone dihedrals separate the clusters (silhouette {{Y12.silhouette_2d|.2f}}); dihedral
  PC1 is essentially one residue, 537, the H11-H12 hinge (`figures/ext/Y12_phipsi_pca.png`).
- What decides the state (`figures/ext/Y13_feature_importance.png`, `Y14_split_hypotheses.png`):
  {{R.cl7_with_peptide}} of {{R.cl7_size}} cl7 chains have a co-crystallised peptide vs {{R.cl4_with_peptide}} of
  {{R.cl4_size}} cl4 chains (entry-level log odds ratio {{R.or_peptide_entry}}); Y537S (log OR {{R.or_Y537S_entry}}) goes
  with cl7, L536S (log OR {{R.or_L536S_entry}}; {{R.cl7_L536S}} cl7 chains) with cl4; ligands with a basic amine go with
  cl4 (ligand-level log OR {{R.or_amine_ligand}}). Quality adds little once peptide and mutation are known (McFadden R2
  gain {{R.quality_gain_r2}}), and cl4 structures are, if anything, better (median resolution {{R.res_cl4_entry|.2f}} vs
  {{R.res_cl7_entry|.2f}} A).
- Inside cl4, the 1.0 A sub-partition separates L536S constructs from the others; they differ at the H11-H12 loop
  ({{Y45.top_positions}}) rather than H12 (`figures/ext/Y45_cl4_substructure.png`).

## 3. Outliers: quality or a different conformation?

- Five detectors flag {{Y15.n_any}} chains at least once; {{Y15.n_consensus}} chains ({{Y15.n_consensus_entries}} entries) are
  flagged by >= 3 (`figures/ext/Y15_outliers.png`).
- Consensus outliers are not worse in resolution, R-free, B or missing residues; {{Y16.n_quality}} of
  {{Y16.n_consensus}} have an individual quality reason, {{Y16.n_conformation}} look like real conformations. They differ
  in experimental conditions (fewer peptides, Y537S and P21 crystals) (`figures/ext/Y16_outlier_quality.png`;
  `results/ext/outlier_verdicts.csv`).
- Experimental metadata (`figures/ext/Y17_metadata_effects.png`): crystal form, construct and peptide track the
  cluster strongly; precipitant class, PEG, glycerol, DMSO and bound additives do not - the solvent does not explain
  the clusters. Crystallisation pH differs (median {{R.ph_cl7_entry}} vs {{R.ph_cl4_entry}}), as part of the different
  crystal recipes.
- Crystal contacts at H12 (`figures/ext/Y44_h12_contacts.png`): median {{Y44.median_cl7|.0f}} lattice atoms touch H12
  in cl7 (the peptide sits there) vs {{Y44.median_cl4|.0f}} in cl4; in cl4 more contacts go with a fully ordered H12
  ({{Y44.median_cl4_complete}} vs {{Y44.median_cl4_partial}}) and a larger H12 shift (Spearman {{Y44.rho_cl4_shift|.2f}}).
- Colour-by maps (`figures/ext/Y18_colour_by_quality.png`, `Y18_colour_by_experiment.png`, `Y18_colour_by_biology.png`,
  interactive `docs/interactive_embedding.html`) show the same thing visually.

## 4. B-factors and the alignment core

- The helix ranking by normalised B is the same in every B-group, cluster and tier; core v3 = the five most stable
  helices ({{R.core_v3_ranges}}). Per chain, normalised B on core v3 is {{Y21.core_v3_median|.2f}} vs
  {{Y21.old_core_median|.2f}} on the old core (`figures/ext/Y21_helix_b.png`, `Y21_core_selection.png`).
- H12 is relatively more mobile in cl4 (H12 / core-v3 B ratio {{Y21.h12_core_v3_ratio_cl4|.2f}} vs
  {{Y21.h12_core_v3_ratio_cl7|.2f}}); PyMOL putty `figures/ext/Y30_pymol_bputty.png`.

## 5. Ligand chemistry and conformation

- Chemistry zoo: {{R.chem_n_runs}} runs; auto-selected {{Y08.balanced_run}} (K = {{Y08.balanced_k}}, silhouette
  {{Y08.balanced_sil|.2f}}, stability {{Y08.balanced_stab|.2f}}) after the declared rule's winner proved degenerate
  (`figures/ext/Y08_chem_zoo_metrics.png`, `Y08_chem_map.png`; C08 needs a human sign-off).
- Ligands of cl4 chains are larger ({{R.heavy_cl4|.0f}} vs {{R.heavy_cl7|.0f}} heavy atoms, MW {{R.mw_cl4}} vs
  {{R.mw_cl7}}) and carry a basic amine far more often ({{R.amine_cl4}}% vs {{R.amine_cl7}}%)
  (`figures/ext/Y09_physchem_by_cluster.png`, `Y09_scaffolds.png`).
- Chemistry and conformation are linked at the ligand level: permutation p {{Y25.p_perm|.1g}}, Cramer's V
  {{Y25.cramers_v|.2f}}; {{Y25.nn_agree|.0%}} of ligands share the H12 state of their nearest chemical neighbour (null
  {{Y25.nn_null|.0%}}) (`figures/ext/Y25_chem_struct_link.png`, `Y25_nn_agreement.png`).
- cl4 ligands are more potent at ESR1 (+{{Y25.potency_raw|.2f}} pChEMBL; +{{Y25.potency_adj|.2f}} within chemistry
  clusters, p {{Y25.potency_adj_p|.2f}}).

## 6. Redundancy: how many structures capture the diversity?

- Effective number of conformations (Hill q1): {{Y26.hill_s08_15.q1|.2f}} at 1.5 A (q0 = {{Y26.hill_s08_15.q0}} clusters);
  q1 stays about 2 for every method, core and subset (`figures/ext/Y26_hill_numbers.png`).
- Two complete-H12 representatives place 90% of all {{R.n_chains}} chains within 1.0 A of one of them;
  at 0.5 A it takes {{Y26.n90_r05|.0f}} (`figures/ext/Y26_set_cover.png`). The earlier cluster-count definition needs
  {{R.n90_cut1p5}} random chains to see 90% of the 1.5 A clusters, because of rare singletons
  (`figures/ext/Y26_rarefaction.png`).
- Copies in one crystal differ by a median {{Y43.med_same_entry|.2f}} A, the same ligand in another entry by
  {{Y43.med_same_ligand|.2f}} A (`figures/ext/Y43_ncs_copies.png`).
- Time: since 2010 only {{Y27.frac_novel_since2010|.1%}} of new chains are > 1.5 A from every earlier one (median
  novelty {{Y27.median_novelty_since2010|.2f}} A); both states have existed since {{Y27.first_year.7}}
  (`figures/ext/Y27_timeline.png`, `Y27_novelty.png`).

## 7. ESR2

- ESR2 ({{Y22.n_chains}} chains, {{Y22.n_entries}} entries) shows the same split: main basin
  {{Y22.peptide_frac_main|.0%}} peptide-bound, second basin {{Y22.peptide_frac_second|.0%}} (`figures/ext/Y22_esr2_parallel.png`).
- In a joint frame, ESR2 agonist-state chains sit {{Y23.block_medians.ESR1cl7-ESR2main|.2f}} A from ESR1 cl7, and the ESR2
  second basin {{Y23.block_medians.ESR1cl4-ESR2second|.2f}} A from ESR1 cl4, while the two ESR1 states are
  {{Y23.block_medians.ESR1cl7-ESR1cl4|.2f}} A apart: the H12 state matters more than the receptor
  (`figures/ext/Y23_joint_embedding.png`). The draft pocket selectivity residues (384, 421) shift < 0.25 A in backbone
  (`figures/ext/Y23_receptor_residue_diff.png`).
- ChEMBL: {{Y24.n_both}} compounds measured at both receptors, Spearman {{Y24.spearman|.2f}}, {{Y24.frac_within_10fold|.0%}}
  within 10-fold (`figures/ext/Y24_potency_esr1_vs_esr2.png`).
- The ESR2 reference could not meet core fraction 1.0 because ESR1 position 466 has no ESR2 counterpart; changing the
  reference does not change the ESR2 clusters (`docs/esr2_core_fraction.md`).

## Storyline: ten figures for the professor

1. **`figures/ext/Y07_rmsd_heatmap_annotated.png`** - The {{R.n_chains}} x {{R.n_chains}} RMSD matrix has two dark blocks.
   Everything that follows explains these two blocks.
2. **`figures/ext/Y04_zoo_grid.png`** - We tried {{R.zoo_n_families}} clustering families; the good ones all find the
   same two clusters. Density and medoid methods fail, and that failure taught us something about the distance itself.
3. **`figures/ext/Y05_compactness.png`** - Both clusters are tight (mean pairwise RMSD < 1 A) and far apart
   ({{Y05.mean_cl7_cl4_A|.1f}} A). cl7 is the more compact one.
4. **`figures/ext/Y10_displacement_profile.png`** - The difference sits in helix 12 and the loop before it
   ({{Y10.share_528_550|.0%}} of it). Elsewhere only two poorly ordered loops move (up to {{Y10.max_d_outside_528|.1f}} A).
5. **`figures/ext/Y30_pymol_medoids.png`** - In 3-D: H12 caps the ligand pocket in cl7 and lies in the coactivator groove
   in cl4. These are the classic agonist and antagonist conformations.
6. **`figures/ext/Y12_phipsi_pca.png`** - Even without any superposition, one backbone angle (residue 537) separates the
   clusters. The switch is a hinge, and residue 537 is where the Y537S construct mutation sits.
7. **`figures/ext/Y14_split_hypotheses.png`** - Coactivator peptide, Y537S vs L536S construct and a basic-amine ligand
   predict the cluster; quality does not. The experiment design decides which state a crystal shows.
8. **`figures/ext/Y16_outlier_quality.png`** - The odd structures are not the bad ones: only {{Y16.n_quality}} of
   {{Y16.n_consensus}} consensus outliers have a quality problem. Solvent and additives explain nothing
   (`figures/ext/Y17_metadata_effects.png`).
9. **`figures/ext/Y25_nn_agreement.png`** - Similar ligands give the same H12 state ({{Y25.nn_agree|.0%}} vs
   {{Y25.nn_null|.0%}} by chance). Chemistry and conformation are linked.
10. **`figures/ext/Y26_set_cover.png`** - Two structures cover 90% of all chains at 1 A; the PDB holds about two
    effective conformations. Most of the {{R.n_chains}} chains are repetition, and ESR2 repeats the same two states
    (`figures/ext/Y23_joint_embedding.png`).

## What is still open

Human steps (ligand labels Q10/Q11 by a teammate, C08 sign-off, E01 citations, W-tasks) are listed in
`docs/HUMAN_STEPS.md`. Tracker state when this file was generated: {{R.tasks_DONE}} DONE, {{R.tasks_MANUAL}} MANUAL,
{{R.tasks_TODO}} TODO of {{R.tasks_total}} tasks.
