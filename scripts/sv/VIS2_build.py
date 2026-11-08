"""VIS2 - collect the valuable visuals into visuals/ (RUN_SV1028).

Pool (step 1): every row of figures/panels/manifest.csv (VIS1): the panels (png_gutter, sv_script, incl. the two
shared legend files), the `single` images and the `not_cut` originals.

[DEFAULT] valuable (step 2), declared here before copying: a graph is kept when it supports a number or claim stated
in docs/RESULTS_SUMMARY.md (incl. the section "Supervisor requests 2026-10-28") or in an SV verdict file; one graph
per claim (the clearest one); keep-all (788-chain) graphs over September (351-entry) ones; the declared 1.5 A
threshold over alternatives; the 10 main-figure files that SV09 lists (incl. the SV13 per-class H12 panels, the SV13
forest plot and the SV14 main figure of the 2026-11-04 addition) are always kept. Excluded: QC/diagnostic graphs tied
to no claim, duplicates, superseded versions and legend-only files; a needed legend goes next to its graph with the
same number and the suffix b. Every pool item gets a decision and a reason in results/sv/VIS2/selection.csv.

Copies (never moves) into visuals/ as NN_<slug>.png, numbered in the story order of RESULTS_SUMMARY.md (The data,
sections 1-7, Storyline, then the SV paragraphs in their order), and writes visuals/INDEX.md.

Claim cells hold one or more items  "quote" (document, location) ; the quote is verbatim text of the document
(line breaks of the document become single spaces; a literal | is written \\| inside the table).

Run: .venv/bin/python -B scripts/sv/VIS2_build.py   (then scripts/sv/VIS2_check.py)
"""
from __future__ import annotations

import csv
import fnmatch
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "figures/panels/manifest.csv"
VIS = ROOT / "visuals"
OUT = ROOT / "results/sv/VIS2"
RS = "docs/RESULTS_SUMMARY.md"
P = "figures/panels/"          # VIS1 PNG-gutter panels
PS = "figures/panels/sv/"      # SV-script panels

# ----------------------------------------------------------------------------------------------------------------
# Selection in story order: (slug, pool key, what it shows, [(quote, document, location)]).
# A slug starting with "+" is a legend placed next to the previous graph (same number, suffix b).
# ----------------------------------------------------------------------------------------------------------------
SEL = [
    # ---- The data
    ("data_precipitant_peg", P + "figures_ext__Y02_metadata_overview__02_main_precipitant_class_n_351.png",
     "Bar chart of the main crystallisation precipitant class of the 351 entries: PEG 323, other classes 2-17 each.",
     [("323 were grown with PEG", RS, "The data")]),
    ("data_construct_mutations", P + "figures_ext__Y02_metadata_overview__05_lbd_mutation_status_per_chain.png",
     "Bar chart of LBD mutation status of the 788 chains: Y537S 443, L536S 243, wild type 58, Cys->Ser/Ala only 30, "
     "D538G 12, other 2.",
     [("443 chains carry Y537S and 243 carry L536S; only 58 are wild type in the LBD", RS, "The data")]),
    # ---- 1. two conformations
    ("dendrogram_1p5_cut", "figures/dendrogram.png",
     "Average-linkage dendrogram of the 788 chains (351 entries) with the 1.5 A cut, which gives 7 clusters.",
     [("Average linkage at 1.5 A gives 7 clusters", RS, "§1")]),
    ("rmsd_heatmap_two_blocks", "figures/ext/Y07_rmsd_heatmap_annotated.png",
     "788 x 788 C-alpha RMSD matrix in leaf order with cluster, QC tier, H12, peptide and construct bars: two blocks, "
     "cl7 445 and cl4 329 chains, median 3.5 A between them.",
     [("The median RMSD between cl7 and cl4 is 3.49 A", RS, "§1"),
      ("The 788 x 788 RMSD matrix has two dark blocks.", RS, "Storyline")]),
    ("embedding_pca", P + "figures_ext__Y03_embeddings__01_pca_trust_0_783_rho_0_87.png",
     "PCA map of the RMSD matrix (trust 0.783, rho 0.87) coloured by S08 cluster, cl7 and cl4 apart; the one of the "
     "five embeddings that carries the cluster key.",
     [("All five embeddings separate the two", RS, "§1")]),
    ("zoo_best_per_family", P + "figures_ext__Y04_zoo_metrics__02_best_setting_of_each_family.png",
     "Selection score of the best setting of each clustering family; spectral sigma=2.0,k=2 ranks first "
     "(K=2, silhouette 0.74, stability 1.00).",
     [("The declared rule selects spectral|sigma=2.0,k=2 (silhouette 0.74, bootstrap stability 1.00", RS, "§1")]),
    ("zoo_pairwise_ari", P + "figures_ext__Y06_agreement__01_pairwise_ari_1_identical_partitions.png",
     "Pairwise ARI between S08 cuts and the best run of each family: 0.89-0.98 for the top families, near 0 for "
     "dbscan, hdbscan, k-medoids, single linkage, affinity propagation and optics.",
     [("The best runs of the top-7 families agree with S08 at ARI 0.89-0.98; density-based, medoid-based, "
       "single-linkage and affinity-propagation methods do not", RS, "§1")]),
    ("cluster_compactness", P + "figures_ext__Y05_compactness__01_s08_1_5_a_clusters_cl7_is_the_tightest.png",
     "Mean pairwise RMSD, mean RMSD to medoid and diameter per S08 cluster at 1.5 A: cl7 0.70 A, cl4 0.91 A mean "
     "pairwise.",
     [("cl7 mean pairwise RMSD 0.70 A (diameter 3.1 A), cl4 0.91 A (diameter 2.4 A); cl4 is looser", RS, "§1")]),
    ("robustness_split_agreement", P + "figures_ext__Y28_robustness_grid__01_cl7_cl4_split_agreement.png",
     "cl7/cl4 split agreement per subset x method on core v3 and the old core (0.98-1.00) and on the fixed-position "
     "and no-H12 controls (0.53-0.97).",
     [("recovers the split with agreement 0.98-1.00", RS, "§1")]),
    ("metric_triangle_inequality", P + "figures_ext__Y42_metric_check__01_s04_breaks_the_triangle_inequality_in_5.png",
     "Histogram of the triangle-inequality excess of the violating triples of the pairwise-common-residue RMSD "
     "(5.4% of triples).",
     [("5.4% of triples break the triangle inequality", RS, "§1")]),
    ("metric_h12less_hubs", P + "figures_ext__Y42_metric_check__02_chains_without_h12_sit_close_to.png",
     "Mean RMSD of each chain to the 787 others vs its observed LBD C-alpha count, by H12 label: chains without "
     "modelled H12 sit lowest.",
     [("chains without H12 look close to everyone", RS, "§1")]),
    # ---- 2. two positions of helix 12
    ("region_rmsd_h12_only", P + "figures_ext__Y42_region_rmsd__04_h12only_10_positions.png",
     "MDS of the RMSD over the 10 H12 positions only: median cl7-cl4 13.97 A, clusters far apart.",
     [("On H12 alone the median cl7-cl4 distance is 14.0 A", RS, "§2")]),
    ("region_rmsd_no_h12", P + "figures_ext__Y42_region_rmsd__02_noh12_223_positions.png",
     "MDS of the RMSD over 223 positions without 528-550: median cl7-cl4 0.82 A, the clusters overlap.",
     [("without residues 528-550 it is 0.82 A and the clusters merge", RS, "§2")]),
    ("displacement_profile", "figures/ext/Y10_displacement_profile.png",
     "cl7-cl4 C-alpha distance of the cluster means per residue (and of the complete-H12 medoids), within-cluster "
     "spread and share of chains observing each position; the large values sit at 528-550.",
     [("Positions 528-550 carry 99% of the squared difference between the cluster means", RS, "§2")]),
    ("h12_shift_by_cluster", P + "figures_ext__Y11_features_by_cluster__01_h12_shift_a_median_0_2_vs_12_3_r_1_00.png",
     "H12 shift from the cl7 mean position by cluster (violins): median cl7 0.2 A vs cl4 12.3 A.",
     [("H12 sits 12.3 A from its cl7 position in cl4 (cl7 0.2 A)", RS, "§2")]),
    ("h12_h3_angle_by_cluster", P + "figures_ext__Y11_features_by_cluster__05_h12_h3_angle_deg_median_81_7_vs_42_5_r.png",
     "H12-H3 axis angle by cluster: median cl7 81.7 deg vs cl4 42.5 deg.",
     [("the H12-H3 axis angle is 82 vs 43 deg", RS, "§2")]),
    ("k362_e542_by_cluster", P + "figures_ext__Y11_features_by_cluster__07_k362_e542_a_median_19_3_vs_13_6_r_1_00.png",
     "K362-E542 C-alpha distance by cluster: median cl7 19.3 A vs cl4 13.6 A.",
     [("K362-E542 C-alpha distance is 19.3 vs 13.6 A", RS, "§2")]),
    ("medoids_superposed_front", P + "figures_ext__Y30_pymol_medoids__01_front_view.png",
     "PyMOL front view of the complete-H12 medoids 5TM1|A|1 (cl7) and 7UJ7|C|1 (cl4) superposed on core v3, H12 "
     "opaque: the two H12 positions.",
     [("This is the agonist (H12 capping the pocket, coactivator groove free) vs antagonist (H12 in the groove) "
       "geometry.", RS, "§2"),
      ("H12 caps the ligand pocket in cl7 and lies in the coactivator groove in cl4.", RS, "Storyline")]),
    ("+medoids_superposed_legend", P + "figures_ext__Y30_pymol_medoids__legend.png",
     "Colour key of the graph with the same number: cl7 medoid 5TM1|A|1 blue, cl4 medoid 7UJ7|C|1 orange "
     "(H12 opaque).",
     [("H12 caps the ligand pocket in cl7 and lies in the coactivator groove in cl4.", RS, "Storyline")]),
    ("phipsi_pca", P + "figures_ext__Y12_phipsi_pca__01_phi_psi_pca_190_positions_silhouette_of.png",
     "PCA of backbone phi/psi over 190 positions (no superposition) coloured by S08 cluster: silhouette of cl7/cl4 "
     "0.66.",
     [("backbone dihedrals separate the clusters (silhouette 0.66)", RS, "§2")]),
    ("phipsi_pc1_loadings", P + "figures_ext__Y12_phipsi_pca__02_h10_11.png",
     "Dihedral PC1 loading per residue along the LBD: one bar, at 537, carries almost all of it.",
     [("PC1 is essentially one residue, 537, the H11-H12 hinge", RS, "§2")]),
    ("split_binary_hypotheses", P + "figures_ext__Y14_split_hypotheses__01_binary_hypotheses_violet_chain_grey.png",
     "Log odds ratio of being in cl7 (95% CI) per binary hypothesis at chain and entry/ligand level: peptide, Y537S, "
     "L536S, D538G, wild type, H12 complete, P21, tier A, basic amine.",
     [("414 of 445 cl7 chains have a co-crystallised peptide vs 4 of 329 cl4 chains (entry-level log odds ratio "
       "6.7); Y537S (log OR 4.2) goes with cl7, L536S (log OR -6.7; 0 cl7 chains) with cl4; ligands with a basic "
       "amine go with cl4 (ligand-level log OR -5.3)", RS, "§2")]),
    ("split_continuous_hypotheses", P + "figures_ext__Y14_split_hypotheses__02_continuous_hypotheses.png",
     "Rank-biserial r (positive = larger in cl7) of continuous variables at chain and entry level: heavy atoms, "
     "resolution, R-free, mean B, pocket size, release year, pH.",
     [("cl4 structures are, if anything, better (median resolution 1.88 vs 2.22 A)", RS, "§2")]),
    ("cl4_subclusters_features", P + "figures_ext__Y45_cl4_substructure__02_what_differs_between_the_cl4_sub.png",
     "Feature heatmap of the cl4 sub-clusters at the 1.0 A cut: L536S share 84% in cl4-1 vs 0% in cl4-2, with H12, "
     "ligand and space-group columns.",
     [("Inside cl4, the 1.0 A sub-partition separates L536S constructs from the others", RS, "§2")]),
    ("cl4_subclusters_profile", P + "figures_ext__Y45_cl4_substructure__03_where_cl4_1_and_cl4_2_differ_largest_at.png",
     "cl4-1 vs cl4-2 mean C-alpha distance per residue: largest at 534, 533 and 339.",
     [("they differ at the H11-H12 loop (534, 533, 339) rather than H12", RS, "§2")]),
    # ---- 3. outliers
    ("outlier_detector_agreement", P + "figures_ext__Y15_outliers__07_how_many_methods_agree.png",
     "Chains by the number of the five outlier detectors that flag them: 0: 431, 1: 244, 2: 91, 3: 13, 4: 9.",
     [("Five detectors flag 357 chains at least once; 22 chains (17 entries) are flagged by >= 3", RS, "§3")]),
    ("outlier_quality", "figures/ext/Y16_outlier_quality.png",
     "Resolution, R-free, mean B and missing C-alpha of consensus outliers vs other chains, effect sizes per outlier "
     "group, and the verdict per consensus outlier (16 conformation, 6 quality); original image, not cut by VIS1.",
     [("Consensus outliers are not worse in resolution, R-free, B or missing residues", RS, "§3"),
      ("6 of 22 have an individual quality reason, 16 look like real conformations", RS, "§3")]),
    ("metadata_vs_cluster", P + "figures_ext__Y17_metadata_effects__01_metadata_vs_structural_cluster.png",
     "Effect size (Cramer's V or epsilon^2) of 25 metadata variables vs the structural cluster at chain and entry "
     "level: space group, mutation, peptide and oligomer state highest.",
     [("crystal form, construct and peptide track the cluster strongly; precipitant class, PEG, glycerol, DMSO and "
       "bound additives do not", RS, "§3")]),
    ("h12_lattice_contacts", P + "figures_ext__Y44_h12_contacts__01_packing_at_h12_cl7_vs_cl4_p_6e_70_cl4.png",
     "Crystal-contact partner atoms within 4 A of H12 (peptide excluded) for cl7, cl4 with complete H12 and cl4 "
     "with partial H12.",
     [("median 1 lattice atoms touch H12 in cl7 (the peptide sits there) vs 10 in cl4; in cl4 more contacts go with "
       "a fully ordered H12 (17.5 vs 8.0)", RS, "§3")]),
    ("h12_contacts_vs_shift_cl4", P + "figures_ext__Y44_h12_contacts__03_within_cl4_contacts_vs_h12_position.png",
     "Within cl4: H12 crystal-contact atoms vs H12 shift from the cl7 mean (Spearman 0.55).",
     [("a larger H12 shift (Spearman 0.55)", RS, "§3")]),
    # ---- 4. B-factors and core
    ("helix_b_ranking", "figures/ext/Y21_helix_b.png",
     "Normalised C-alpha B per helix by B-group, cluster and tier (ranks), per-chain B on core v3 / old core / H12 "
     "by cluster, and the median B profile; original image, not cut by VIS1.",
     [("The helix ranking by normalised B is the same in every B-group, cluster and tier", RS, "§4"),
      ("H12 is relatively more mobile in cl4 (H12 / core-v3 B ratio 1.53 vs 1.33)", RS, "§4")]),
    ("core_v3_selection", P + "figures_ext__Y21_core_selection__01_k04_rule_helices_in_the_better_half_in.png",
     "Rank of median normalised B per helix in the low, mid and high B-groups: the five core-v3 helices stay in "
     "ranks 1-5, H12 at rank 8.",
     [("core v3 = the five most stable helices (338-362; 371-396; 420-439; 441-456; 496-531)", RS, "§4")]),
    ("core_v3_vs_old_core_b", P + "figures_ext__Y21_core_selection__02_per_chain_n_788_median_old_core_1_038.png",
     "Per-chain normalised B on core v3 vs the old core for 788 chains: median 0.928 vs 1.038.",
     [("Per chain, normalised B on core v3 is 0.93 vs 1.04 on the old core", RS, "§4")]),
    # ---- 5. ligand chemistry
    ("chem_zoo_runs", P + "figures_ext__Y08_chem_zoo_metrics__03_all_chemistry_runs.png",
     "Silhouette vs bootstrap stability of all 276 chemistry clustering runs, the degenerate declared-rule winner "
     "and the auto-selected run marked.",
     [("Chemistry zoo: 276 runs; auto-selected MACCS|kmeans_mds|k=6 (K = 6, silhouette 0.29, stability 0.97) after "
       "the declared rule's winner proved degenerate", RS, "§5")]),
    ("ligand_physchem_by_cluster", "figures/ext/Y09_physchem_by_cluster.png",
     "Ligand descriptors of cl7 (192) vs cl4 (104) ligands: MW, cLogP, TPSA, heavy atoms, rotatable bonds, H-bond "
     "donors/acceptors, rings, basic-amine share (1% vs 74%); original image, not cut by VIS1.",
     [("Ligands of cl4 chains are larger (34 vs 24 heavy atoms, MW 470 vs 344) and carry a basic amine far more "
       "often (74% vs 1%)", RS, "§5")]),
    ("chem_struct_permutation", P + "figures_ext__Y25_chem_struct_link__02_cramer_s_v_0_43_95_ci_0_40_0_55_mi_0_31.png",
     "Chi-square of chemistry cluster x structural majority under 10,000 ligand-level permutations vs the observed "
     "170 (permutation p 0.0001; Cramer's V 0.43).",
     [("permutation p 0.0001, Cramer's V 0.43", RS, "§5")]),
    ("nn_ligand_h12_agreement", P + "figures_ext__Y25_nn_agreement__01_do_similar_ligands_give_the_same_h12.png",
     "Share of ligands whose nearest chemical neighbour has the same structural majority, by ECFP4 Tanimoto bin: "
     "overall 0.89 vs permutation null 0.51.",
     [("89% of ligands share the H12 state of their nearest chemical neighbour (null 51%)", RS, "§5")]),
    ("potency_cl4_vs_cl7", P + "figures_ext__Y25_chem_struct_link__03_potency_cl4_cl7_1_13_raw_p_3e_05_0_83.png",
     "Median pChEMBL at ESR1 of cl7 (27) and cl4 (53) ligands by chemistry cluster: cl4-cl7 +1.13 raw, +0.83 within "
     "chemistry (p 0.02).",
     [("cl4 ligands are more potent at ESR1 (+1.13 pChEMBL; +0.83 within chemistry clusters, p 0.02)", RS, "§5")]),
    # ---- 6. redundancy
    ("hill_numbers", P + "figures_ext__Y26_hill_numbers__01_hill_numbers_many_clusters_2_effective.png",
     "Hill numbers q0 (clusters), q1 and q2 (effective conformations) for S08 cuts, other methods and subsets: q1 "
     "about 2 throughout.",
     [("Effective number of conformations (Hill q1): 2.19 at 1.5 A (q0 = 7 clusters); q1 stays about 2 for every "
       "method, core and subset", RS, "§6")]),
    ("set_cover", "figures/ext/Y26_set_cover.png",
     "Greedy set cover: representatives needed vs coverage radius for 90/95/99% of the 788 chains, and per subset "
     "and core at 1.0 A; original image, not cut by VIS1.",
     [("Two complete-H12 representatives place 90% of all 788 chains within 1.0 A of one of them; at 0.5 A it takes "
       "48", RS, "§6")]),
    ("cluster_discovery_curve", "figures/diversity_curve.png",
     "Share of S08 clusters seen vs chains added in random order (100 orderings, 788 chains) at cuts 1.0-2.5 A: at "
     "1.5 A, 90% at N=488.",
     [("488 random chains to see 90% of the 1.5 A clusters, because of rare singletons", RS, "§6")]),
    ("ncs_copies", P + "figures_ext__Y43_ncs_copies__01_pair_classes.png",
     "C-alpha RMSD distributions of NCS copies (median 0.54 A), the same ligand in other entries (0.65 A) and other "
     "pairs (1.38 A).",
     [("Copies in one crystal differ by a median 0.54 A, the same ligand in another entry by 0.65 A", RS, "§6")]),
    ("novelty_by_year", P + "figures_ext__Y27_novelty__01_each_new_structure_vs_everything.png",
     "Novelty of each chain (minimum RMSD to earlier complete-H12 chains) by release year, coloured by cluster, with "
     "the yearly median and the 1.5 A line.",
     [("since 2010 only 0.4% of new chains are > 1.5 A from every earlier one (median novelty 0.38 A)", RS, "§6")]),
    ("release_timeline", "figures/ext/Y27_timeline.png",
     "Chains released per year by cluster (cl7 and cl4 both present in 1998), with median resolution and R-free per "
     "release year.",
     [("both states have existed since 1998", RS, "§6")]),
    # ---- 7. ESR2
    ("esr2_peptide", P + "figures_ext__Y22_esr2_parallel__03_coloured_by_co_crystallised_peptide.png",
     "ESR2 MDS map (86 chains) coloured by co-crystallised peptide (63 with, 23 without): the main basin carries the "
     "peptide.",
     [("ESR2 (86 chains, 39 entries) shows the same split: main basin 93% peptide-bound, second basin 0%", RS,
       "§7")]),
    ("esr2_esr1_state_distances", P + "figures_ext__Y23_joint_embedding__02_receptor_vs_state_which_distance_is.png",
     "Median C-alpha RMSD between ESR1 clusters and ESR2 basins in the joint frame: cl7-ESR2 main 1.33 A, cl4-ESR2 "
     "second 2.31 A, cl7-cl4 3.49 A.",
     [("ESR2 agonist-state chains sit 1.33 A from ESR1 cl7, and the ESR2 second basin 2.31 A from ESR1 cl4, while "
       "the two ESR1 states are 3.49 A apart", RS, "§7")]),
    ("esr1_esr2_residue_shift", "figures/ext/Y23_receptor_residue_diff.png",
     "Per-residue distance between ESR1 and ESR2 mean C-alpha at the same H12 state, sequence differences flagged; "
     "draft selectivity positions 384 and 421 marked.",
     [("The draft pocket selectivity residues (384, 421) shift < 0.25 A in backbone", RS, "§7")]),
    ("esr1_vs_esr2_potency", P + "figures_ext__Y24_potency_esr1_vs_esr2__01_potency_at_the_two_receptors_spearman_0.png",
     "Median pChEMBL at ESR1 vs ESR2 for 1615 ChEMBL compounds with the 10-fold band; PDB ligands of cl7 and cl4 "
     "chains marked.",
     [("1615 compounds measured at both receptors, Spearman 0.50, 54% within 10-fold", RS, "§7")]),
    # ---- Supervisor requests 2026-10-28 (SV paragraphs in their order)
    ("sv01a_completeness_profile", "figures/sv/SV01a/SV01a_completeness_profile.png",
     "Fraction of chains with C-alpha observed at residues 520-550, cl7 vs cl4, chain- and entry-weighted; H12 "
     "window 538-547 shaded.",
     [("H12 538-547 is complete in 431 of 445 cl7 chains (0.969) but in 120 of 329 cl4 chains (0.365); 207 of the "
       "212 partial H12 lack only C-terminal residues", RS, "SV01a")]),
    ("sv01b_forest_completeness", "figures/sv/SV01b/SV01b_forest_completeness.png",
     "Odds ratios (95% CI) for complete H12 538-547 per predictor, blocks Q / X (with lattice terms) / B, within cl7 "
     "and cl4, chain and entry level, with the block gains; original image, not cut by VIS1.",
     [("H12 complete, chain level, Q+X vs B gain cl7 0.414 vs 0.225 and cl4 0.523 vs 0.079", RS, "SV01b"),
      ("per SD of symmetry-mate contact atoms the odds of a complete H12 rise 20.2-fold in cl4 (12.6-fold in cl7)",
       RS, "SV01b")]),
    ("sv01b_forest_mobility", "figures/sv/SV01b/SV01b_forest_mobility.png",
     "Coefficients (95% CI) for H12 mobility (C-alpha B z-score of 528-547 vs core v3) per predictor, blocks Q / X / "
     "B, within cl7 and cl4, with the block gains; original image, not cut by VIS1.",
     [("H12 mobility cl7 0.184 vs 0.076 and cl4 0.210 vs 0.143", RS, "SV01b"),
      ("H12 mobility falls by 0.268 z in cl4 (0.17 in cl7)", RS, "SV01b")]),
    ("sv02_helix_track", "figures/sv/SV02/SV02_helix_track.png",
     "pydssp helix fraction per residue for all chains, cl7 and cl4, with the 8 consensus helices, the historical "
     "H12 windows and the share of chains observing each position.",
     [("The data-derived H12 window is 538-544 (cl7 538-545, cl4 537-544)", RS, "SV02"),
      ("H10/11 ends at 530 in cl7 and 525 in cl4, H9 starts at 473 in cl7 and 466 in cl4", RS, "SV02")]),
    ("sv02_rama_537", "figures/sv/SV02/rama/SV02_rama_537.png",
     "Ramachandran plot of residue 537 by cluster (cl7 429, cl4 304, minor 11 chains): cl7 psi near +160 deg, cl4 "
     "near -50 deg.",
     [("Residue 537 alone separates the clusters: circular-mean psi +156.6 deg in cl7 vs -46.8 deg in cl4 "
       "(separation 157.2 deg; overlap 0 at chain and 0 at entry level)", RS, "SV02")]),
    ("sv03_h12_h3_angle", PS + "SV03/SV03_descriptors_by_cluster__angle_H12_H3_deg.png",
     "Angle H12-H3 by S08 cluster in the core v3 frame (cl7 n=431, cl4 n=307, minor n=10).",
     [("the largest are the H12-H3 angle (median 79.05 vs 46.63 deg)", RS, "SV03")]),
    ("sv03_h12_h11axis_distance", PS + "SV03/SV03_descriptors_by_cluster__dist_H12c_H11axis_A.png",
     "Distance from the H12 centroid to the H11 axis by S08 cluster (cl7 n=431, cl4 n=307, minor n=10).",
     [("the distance from the H12 centroid to the H11 axis (12.48 vs 22.26 A)", RS, "SV03")]),
    ("sv04_h12_map_partial_a", P + "figures_sv_SV04__h12_map_partial__01_a_h12_map_549_complete_h12_chains_at.png",
     "(a) H12 map: 549 complete-H12 chains at SV03 positions (filled), 199 partial-H12 chains at SV04 midpoint "
     "estimates (hollow); colour = cluster, marker = ligand chemotype.",
     [("Of 212 partial-H12 chains, 199 have a usable H12 midpoint, and every one of the 191 in cl7 or cl4 lands on "
       "its own S08 side", RS, "SV04")]),
    ("sv04_h12_map_partial_b", P + "figures_sv_SV04__h12_map_partial__02_b_same_chains_all_at_sv04_estimates_549.png",
     "(b) The same 549 complete and 199 partial-H12 chains, all placed at SV04 midpoint estimates.",
     [("panel (b) uses the SV04 estimator for all chains, because midpoint and centroid differ by a "
       "cluster-dependent offset in x (complete cl4 +1.51 A, cl7 -0.68 A)", RS, "SV04")]),
    ("sv05_core_v4_spread", P + "figures_sv_SV05__SV05_core_v4__01_a_128_of_132_candidate_positions.png",
     "(a) C-alpha spread after iterative fit per position: 128 of 132 candidate positions have spread <= 0.5 A "
     "(core v4).",
     [("iterative superposition keeps 128 positions with spread at most 0.5 A as core v4", RS, "SV05")]),
    ("sv05_core_v4_membership", P + "figures_sv_SV05__SV05_core_v4__02_b_core_positions_v3_123_v4_128_jaccard.png",
     "(b) Positions of core v3 (123), core v4 (128) and the lowest-spread 123 along the sequence (Jaccard v4-v3 "
     "0.64).",
     [("it adds H1 and H9, drops 25 v3 positions at helix ends, and has Jaccard 0.64 with v3", RS, "SV05")]),
    ("sv05_displacement_v3_v4", P + "figures_sv_SV05__SV05_core_v4__03_c_cluster_mean_distance_cl7_cl4_per.png",
     "(c) cl7-cl4 cluster-mean distance per position with core v3 and core v4: 528-550 hold 98.6% and 98.5% of the "
     "squared difference.",
     [("98.5 % of the squared cl7-cl4 difference stays in 528-550 (v3 98.6 %)", RS, "SV05")]),
    ("sv05_clusters_ari", P + "figures_sv_SV05__SV05_core_v4__04_d_clusters_per_cut_ari_of_v4_against.png",
     "(d) Clusters per average-linkage cut with core v3 (S08) and core v4, with ARI and chains moved (1.5 A: ARI "
     "0.990, 3 chains moved).",
     [("ARI with S08 0.990 at 1.5 A, 3 chains change cluster", RS, "SV05")]),
    ("sv06_outliers_map", "figures/sv/SV06/SV06_outliers_map.png",
     "(A) The 22 consensus outliers on the shared MDS map (2 cl7, 12 cl4, 8 minor); (B) their RMSD to the own vs the "
     "other complete-H12 medoid; original image, not cut by VIS1.",
     [("2 in cl7 against 8.0 expected (q 0.0023)", RS, "SV06"),
      ("None (0 of 14) of the main-cluster outliers is closer to the other cluster's medoid than to its own", RS,
       "SV06")]),
    ("sv07_spike_profile", PS + "SV07/SV07_spike_profile__profile.png",
     "(a) Y10 cl7-cl4 C-alpha distance profile (symlog) with the 1.0 A threshold and medoid curves; every spike "
     "labelled with its audit verdict.",
     [("19 positions outside 528-550 exceed 1.0 A in the Y10 profile", RS, "SV07"),
      ("14 are genuine loop differences (331-340, 464, 466-468)", RS, "SV07")]),
    ("sv07_observation_share", PS + "SV07/SV07_spike_profile__observation.png",
     "(b) Share of cl7 and cl4 chains observing each position (SV00 mask), spike positions as thin lines; lowest at "
     "a spike: 462 (cl7 12%, cl4 30%).",
     [("Y10 already averages only the chains that observe a position", RS, "SV07")]),
    ("sv11_pocket_site_share", "figures/sv/SV11/SV11_site_share.png",
     "Share of chains with a P2Rank pocket at the ligand pocket, the AF-2 groove and the dimer interface, cl7 vs "
     "cl4.",
     [("A pocket in the AF-2 groove is found in 0.53 of cl7 chains (234 of 445) but 0.02 of cl4 chains (5 of 329",
       RS, "SV11"),
      ("a dimer-interface pocket in 0.18 vs 0.31 of chains", RS, "SV11")]),
    ("sv13_forest_c1_complete", PS + "SV13/SV13_forest__C1_logit.png",
     "C1 rows (H12 complete and fraction observed; SERD vs SERM and antagonist vs SERM, all chains and within cl4): "
     "each row shows its coverage-gate reason in place of an estimate.",
     [("Every class comparison fails the coverage gate (0 of 22 comparison rows testable): SERD vs SERM (SERD 4 "
       "entries, all in cl4), antagonist vs SERM (0 antagonist entries)", RS, "SV13")]),
    ("sv13_forest_c1_mobility_c2a", PS + "SV13/SV13_forest__C1_mobility_C2a.png",
     "C1 mobility rows and C2a rows (partial agonist vs agonist within cl7: mobility, H12 geometry, psi 537): each "
     "row shows its coverage-gate reason.",
     [("partial agonist vs agonist (0 partial-agonist entries)", RS, "SV13")]),
    ("sv13_forest_c2b", PS + "SV13/SV13_forest__C2b_rho.png",
     "C2b rows (agonists within cl7, Spearman rho of each H12 outcome with the SV15 dpX): each row shows its "
     "coverage-gate reason.",
     [("agonist affinity vs H12 within cl7 (3 agonist entries with a tested SV15 dpX)", RS, "SV13")]),
    ("sv13_h12_panel_agonist", PS + "SV13/SV13_h12_class_panels__agonist.png",
     "Agonist class (55 chains, 27 entries): representative 5WGD|B|1 with EST, its 528-550 coloured by H12 B z, and "
     "the 528-550 traces of all class chains by cluster.",
     [("`5WGD|B|1`, agonist, H12 10 of 10 residues modelled", RS, "SV13"),
      ("complete H12 538-547 in agonist 53 of 55 chains", RS, "SV13")]),
    ("sv13_h12_panel_serm", PS + "SV13/SV13_h12_class_panels__SERM.png",
     "SERM class (53 chains, 20 entries): representative 2BJ4|A|1 with OHT, H12 not modelled, and the 528-550 traces "
     "of all class chains by cluster.",
     [("`2BJ4|A|1`, SERM, H12 not modelled", RS, "SV13"),
      ("SERM 19 of 53", RS, "SV13")]),
    ("sv13_h12_panel_serd", PS + "SV13/SV13_h12_class_panels__SERD.png",
     "SERD class (9 chains, 4 entries): representative 7MSA|C|1 with ZNM, H12 7 of 10 residues modelled, and the "
     "528-550 traces of the class chains.",
     [("`7MSA|C|1`, SERD, 7 of 10", RS, "SV13"),
      ("SERD 5 of 9", RS, "SV13")]),
    ("sv13_h12_class_overlay", PS + "SV13/SV13_h12_class_panels__overlay.png",
     "Overlay of the three class representatives (agonist 5WGD|B|1, SERM 2BJ4|A|1, SERD 7MSA|C|1) in the same "
     "view; matplotlib C-alpha traces, ligands as sticks.",
     [("one view for all panels", RS, "SV13"),
      ("drawn with matplotlib because PyMOL is not available", RS, "SV13")]),
    ("sv13_class_mds_map", "figures/sv/SV13/SV13_mds_map.png",
     "Shared MDS map (788 chains), fill = S08 cluster, marker = SV12 action class: agonist 55, SERM 53, SERD 9 "
     "chains, 671 unlabelled.",
     [("agonist entries are 25 of 27 in cl7 and SERM entries 19 of 20 in cl4, so this restates the clusters", RS,
       "SV13")]),
    ("sv14_pocket_vs_ligand", "figures/sv/SV14/SV14_pocket_vs_ligand.png",
     "Pocket vs ligand volume per chain (701 chains) as % of the E2 references; fill = cluster, marker = action "
     "class; equal-filling lines at 25/50/75%.",
     [("Ligand and pocket volume go together (Spearman rho 0.44 at ligand level, 0.47 at entry level), within cl7 "
       "(0.41 / 0.46) but not within cl4 (0.11 / 0.14)", RS, "SV14")]),
    ("sv14_filling", "figures/sv/SV14/SV14_filling.png",
     "Pocket filling % (entry medians) by cluster and by action class; original image with two graphs, not cut by "
     "VIS1.",
     [("Filling has a median of 33.8 % and does not differ between cl7 and cl4 (entry medians 33.0 vs 32.3 %, q "
       "0.14) or between agonist and SERM (q 0.94)", RS, "SV14")]),
    ("sv15_dpx_by_group", "figures/sv/SV15/SV15_dpx_by_group.png",
     "dpX (ligand pX minus E2 pX in the same binding assay) per ligand by SV10 group and by SV12 class; relative "
     "filled, absolute-only hollow.",
     [("dpX does not differ significantly between cl7 and cl4 ligands (median -2.17 vs -0.23, q 0.15", RS,
       "SV15")]),
    ("sv15_dpx_vs_filling", "figures/sv/SV15/SV15_dpx_vs_filling.png",
     "dpX vs SV14 pocket filling % for 14 ligands with a relative dpX (filled, labelled) and 48 absolute-only "
     "ligands (hollow).",
     [("does not follow pocket filling (rho -0.03)", RS, "SV15")]),
]

# ----------------------------------------------------------------------------------------------------------------
# Exclusion reasons for the rest of the pool, first match wins (pattern on the pool key, then on the source figure).
# {slug} is replaced by the visuals/ file name of that selected graph.
# ----------------------------------------------------------------------------------------------------------------
EXCL = [
    # superseded / September versions
    ("results/archive_v2_strict/*", "superseded: archived strict-filter version"),
    ("results/figures/*", "September (351-entry) core-task version; keep-all graphs preferred, no claim of its own"),
    (P + "results_figures__*", "September (351-entry) core-task version; keep-all graphs preferred, no claim of its own"),
    (P + "results_archive_v2_strict__*", "superseded: archived strict-filter version"),
    ("figures/provisional_*", "superseded: September (351-entry) provisional version; keep-all graphs preferred"),
    (P + "figures__provisional_*", "superseded: September (351-entry) provisional version; keep-all graphs preferred"),
    ("figures/qc_outliers.png", "September (351-entry) QC graph; no stated claim"),
    (P + "figures__qc_outliers__*", "September (351-entry) QC graph; no stated claim"),
    ("figures/0*_*.png", "September (351-entry) core-task figure; keep-all graphs preferred, no claim of its own"),
    ("figures/1*_*.png", "September (351-entry) core-task figure; keep-all graphs preferred, no claim of its own"),
    (P + "figures__0*", "September (351-entry) core-task figure; keep-all graphs preferred, no claim of its own"),
    (P + "figures__1*", "September (351-entry) core-task figure; keep-all graphs preferred, no claim of its own"),
    ("esr_analysis/resources/reference/reference_context.png", "no stated claim (reference structure render)"),
    ("*reference_render/reference_context.png", "no stated claim (reference structure render)"),
    # other keep-all figures in figures/
    ("figures/helix_bfactor_stability.png", "same claim as {helix_b_ranking}, which also covers cluster and tier"),
    ("figures/potency_by_cluster.png", "same claim as {potency_cl4_vs_cl7}, which states the numbers"),
    ("figures/esr2_pca_landscape.png", "same claim as {esr2_peptide} (ESR2 split), which shows the peptide"),
    ("figures/structural_sensitivity.png", "no stated claim (threshold sensitivity; the declared 1.5 A cut is used)"),
    (P + "figures__structural_sensitivity__*", "no stated claim (threshold sensitivity; the declared 1.5 A cut is used)"),
    ("figures/pocket_pca.png", "no stated claim (P2Rank descriptor PCA)"),
    ("figures/deposits_per_year.png", "no stated claim (deposits per year, ESR1 and ESR2)"),
    ("figures/displacement_heatmap.png", "no stated claim (ESR1-ESR2 pair displacement; {esr1_esr2_residue_shift} carries the claim)"),
    ("figures/displacement_structure.png", "no stated claim (structure render)"),
    ("figures/medoids.png", "no stated claim (render of cluster medoids; {medoids_superposed_front} carries the H12 claim)"),
    # PyMOL
    ("figures/ext/pymol_panels/medoids_*", "duplicate: raw PyMOL render of a Y30 medoid view (no key); {medoids_superposed_front} kept"),
    ("figures/ext/pymol_panels/bputty_*", "duplicate: raw PyMOL render of a Y30 B-putty panel; see the B-putty reason"),
    (P + "figures_ext__Y30_pymol_bputty__legend.png", "legend-only file whose graphs are not kept"),
    (P + "figures_ext__Y30_pymol_bputty__*", "same claim as {helix_b_ranking} (H12 / core-v3 B ratio 1.53 vs 1.33 over all chains); the putty shows two medoids"),
    (P + "figures_ext__Y30_pymol_medoids__02_*", "duplicate: the same superposition rotated 90 deg; {medoids_superposed_front} kept"),
    # ext figures, duplicates of a kept claim
    (P + "figures_ext__Y02_metadata_overview__06_*", "no stated claim shown: the title gives the peptide share, the graph shows chains per ASU"),
    (P + "figures_ext__Y03_embeddings__*", "same claim as {embedding_pca} (all five embeddings separate the two); the PCA panel carries the cluster key"),
    (P + "figures_ext__Y04_zoo_grid__*", "same claim as {zoo_pairwise_ari} (which families find the split); one run per panel"),
    ("figures/ext/Y06_alluvial.png", "same claim as {zoo_pairwise_ari}"),
    (P + "figures_ext__Y04_zoo_metrics__01_*", "same claim as {zoo_best_per_family}, which names the selected run with its numbers"),
    (P + "figures_ext__Y06_agreement__02_*", "duplicate: AMI repeats the ARI matrix of {zoo_pairwise_ari}"),
    (P + "figures_ext__Y06_agreement__03_*", "same claim as {zoo_pairwise_ari} (ARI of the selected run with S08)"),
    (P + "figures_ext__Y28_robustness_grid__02_*", "duplicate: ARI vs S08 restates the split agreement of {robustness_split_agreement}"),
    (P + "figures_ext__Y42_metric_check__03_*", "same claim as {region_rmsd_h12_only} and {region_rmsd_no_h12}, which give the stated distances"),
    ("figures/ext/Y07_dendrogram_metadata.png", "same claims as {dendrogram_1p5_cut} (1.5 A cut) and {metadata_vs_cluster} (what tracks the cluster)"),
    (P + "figures_ext__Y12_phipsi_pca__03_*", "same claim as {phipsi_pc1_loadings} (residue 537)"),
    (P + "figures_ext__Y13_feature_importance__*", "same claim as {split_binary_hypotheses} / {split_continuous_hypotheses}, which give the stated effects"),
    ("figures/ext/Y13_decision_tree.png", "same claim as {split_binary_hypotheses}"),
    (P + "figures_ext__Y14_split_hypotheses__03_*", "same claim as {split_binary_hypotheses} (strength of evidence restates the tests)"),
    (P + "figures_ext__Y15_outliers__06_*", "same content as {sv06_outliers_map} (consensus outliers on the shared MDS map)"),
    (P + "figures_ext__Y15_outliers__*", "diagnostic: single-detector map or overlap of flagged sets; {outlier_detector_agreement} carries the claim"),
    (P + "figures_ext__Y17_metadata_effects__02_*", "same claim as {outlier_quality} (outliers differ in experimental conditions)"),
    (P + "figures_ext__Y17_metadata_effects__03_*", "same claim as {metadata_vs_cluster} (crystal form tracks the cluster)"),
    (P + "figures_ext__Y18_*", "colour-by map: RESULTS_SUMMARY section 3 says these maps 'show the same thing "
     "visually'; that claim is kept with {metadata_vs_cluster}"),
    (P + "figures_ext__Y45_cl4_substructure__01_*", "same claim as {cl4_subclusters_features} (the 1.0 A cl4 sub-partition)"),
    (P + "figures_ext__Y08_chem_map__*", "same claim as {chem_zoo_runs} (selected chemistry clustering)"),
    (P + "figures_ext__Y08_chem_zoo_metrics__*", "same claim as {chem_zoo_runs}"),
    ("figures/ext/Y09_scaffolds.png", "same claim as {ligand_physchem_by_cluster}, which gives the numbers"),
    (P + "figures_ext__Y25_chem_struct_link__01_*", "same claim as {chem_struct_permutation}, which gives p and Cramer's V"),
    (P + "figures_ext__Y26_rarefaction__01_*", "same claim as {cluster_discovery_curve}, which states 90% at N=488"),
    (P + "figures_ext__Y26_rarefaction__02_*", "same claim as {set_cover}"),
    (P + "figures_ext__Y22_esr2_parallel__01_*", "same claim as {esr2_peptide} (the ESR2 basins)"),
    (P + "figures_ext__Y23_joint_embedding__01_*", "same claim as {esr2_peptide} and {esr2_esr1_state_distances} (ESR2 repeats the two states)"),
    (P + "figures_ext__Y44_h12_contacts__02_*", "same claim as {h12_lattice_contacts} (the peptide sits on H12 in cl7)"),
    # SV figures
    ("figures/sv/SV01a/SV01a_forest_completeness.png", "superseded: first pass without lattice terms; SV01b refit {sv01b_forest_completeness} kept"),
    ("figures/sv/SV02/SV02_sse_heatmap.png", "same claims as {sv02_helix_track} (helix ends per cluster)"),
    ("figures/sv/SV02/rama/*", "no stated claim (only residue 537 is cited; {sv02_rama_537})"),
    (PS + "SV03/*", "covered by the SV03 claim of {sv03_h12_h3_angle} and {sv03_h12_h11axis_distance} (the two largest of the 16 descriptors)"),
    ("figures/sv/SV03/SV03_h12_map.png", "superseded: the same map with the SV04 overlay, {sv04_h12_map_partial_a}"),
    (P + "figures_sv_SV05__SV05_core_v4__05_*", "same claim as {sv05_clusters_ari} (3 chains change cluster at 1.5 A)"),
    ("figures/sv/SV06/SV06_rmsd_heatmap_outlier_bar.png", "same claim as {sv06_outliers_map} (outliers per cluster)"),
    (PS + "SV13/SV13_h12_class_panels__antagonist.png", "empty panel: no chain carries the class (text only)"),
    (PS + "SV13/SV13_h12_class_panels__partial_agonist.png", "empty panel: no chain carries the class (text only)"),
    (PS + "SV13/SV13_h12_class_panels__other.png", "empty panel: no chain carries the class (text only)"),
    ("figures/sv/SV13/SV13_h12_map.png", "same claim as {sv13_class_mds_map} (classes follow the clusters)"),
]
DEFAULT_REASON = "no stated claim in RESULTS_SUMMARY.md or an SV verdict file (QC, diagnostic or context graph)"


def sha1(path: Path) -> str:
    return hashlib.sha1(path.read_bytes()).hexdigest()


def esc(text: str) -> str:
    return text.replace("|", "\\|")


def load_pool():
    pool = []
    with MANIFEST.open(newline="") as fh:
        for row in csv.DictReader(fh):
            sources = row["source_path(s)"].split(";")
            if row["method"] in ("single", "not_cut"):
                key = sources[0]
            else:
                key = row["panel_path"]
            pool.append({"key": key, "method": row["method"], "source": sources[0],
                         "all_sources": row["source_path(s)"], "title": row["panel_title"],
                         "sha1_manifest": row["sha1"]})
    return pool


def main() -> int:
    pool = load_pool()
    by_key = {p["key"]: p for p in pool}
    if len(by_key) != len(pool):
        sys.exit("duplicate pool keys in the manifest")

    # numbering
    plan, n = [], 0
    for slug, key, what, claims in SEL:
        if key not in by_key:
            sys.exit(f"not in the VIS1 pool: {key}")
        if slug.startswith("+"):
            nn, slug = f"{n:02d}b", slug[1:]
        else:
            n += 1
            nn = f"{n:02d}"
        plan.append({"nn": nn, "slug": slug, "file": f"{nn}_{slug}.png", "key": key, "what": what,
                     "claims": claims})
    fname = {p["slug"]: p["file"] for p in plan}
    keys_sel = {p["key"]: p for p in plan}

    # visuals/: refuse unknown files (never delete), then copy
    VIS.mkdir(exist_ok=True)
    planned = {p["file"] for p in plan} | {"INDEX.md"}
    stray = sorted(f.name for f in VIS.iterdir() if f.name not in planned)
    if stray:
        sys.exit(f"visuals/ holds files outside the plan, not touched: {stray}")
    for p in plan:
        src = ROOT / p["key"]
        dst = VIS / p["file"]
        shutil.copy2(src, dst)
        if sha1(src) != sha1(dst):
            sys.exit(f"copy mismatch {src}")

    # selection table over the whole pool
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for item in pool:
        k = item["key"]
        if k in keys_sel:
            rows.append({**item, "decision": "included", "visual_file": keys_sel[k]["file"],
                         "reason": "supports a stated claim (INDEX.md)"})
            continue
        reason = DEFAULT_REASON
        for pat, text in EXCL:
            if fnmatch.fnmatchcase(k, pat) or fnmatch.fnmatchcase(item["source"], pat):
                reason = text
                break
        for slug, f in fname.items():
            reason = reason.replace("{" + slug + "}", f)
        if "{" in reason:
            sys.exit(f"unresolved placeholder in reason for {k}: {reason}")
        rows.append({**item, "decision": "excluded", "visual_file": "", "reason": reason})
    with (OUT / "selection.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["key", "method", "source", "title", "decision", "visual_file", "reason"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # INDEX.md
    n_graphs = sum(1 for p in plan if not p["nn"].endswith("b"))
    n_leg = len(plan) - n_graphs
    lines = [
        "# Visuals index",
        "",
        "Valuable graphs of the ESR1 project, copied (not moved) from the VIS1 pool (`figures/panels/manifest.csv`: "
        "panels, `single` images, `not_cut` originals) by `scripts/sv/VIS2_build.py`, numbered in the story order of "
        "`docs/RESULTS_SUMMARY.md` (The data, sections 1-7, Storyline, then the SV paragraphs of 'Supervisor requests "
        "2026-10-28').",
        f"{n_graphs} graphs and {n_leg} legend (suffix b). Rule: one graph per claim stated in RESULTS_SUMMARY.md or "
        "an SV verdict file, keep-all (788-chain) over September (351-entry) versions, the declared 1.5 A threshold, "
        "and every main SV figure of SV09. Decisions for the whole pool: `results/sv/VIS2/selection.csv`.",
        "Claims are verbatim quotes (line breaks of the document read as spaces); `scripts/sv/VIS2_check.py` checks "
        "them, the files and their SHA-1.",
        "",
        "| NN | file | what it shows | claim supported | source path |",
        "|---|---|---|---|---|",
    ]
    for p in plan:
        item = by_key[p["key"]]
        claim_cell = " ; ".join(f'"{esc(q)}" ({doc}, {loc})' for q, doc, loc in p["claims"])
        src_cell = f"`{p['key']}`"
        if item["method"] in ("png_gutter", "sv_script"):
            src_cell += f" (panel of `{item['source']}`)"
        lines.append(f"| {p['nn']} | `{p['file']}` | {esc(p['what'])} | {claim_cell} | {src_cell} |")
    (VIS / "INDEX.md").write_text("\n".join(lines) + "\n")

    n_inc = sum(r["decision"] == "included" for r in rows)
    print(f"pool {len(pool)}; included {n_inc} ({n_graphs} graphs + {n_leg} legend); excluded {len(rows) - n_inc}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
