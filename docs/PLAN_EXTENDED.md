# Extended analysis plan (final run, 2026-10-25)

Professor's guidance at the last meeting: do not stick to the original description; try many clustering methods
(scikit-learn zoo), report **how compact** clusters are, auto-evaluate chemistry clusterings and pick the best by a
metric, check whether outliers come from **structure quality** or **experimental conditions**, colour plots by
R-free / B-factors / solvents / anything else, and find out **why there are two clusters**. He wants graphs.

Standing rules: tracker tasks Y01-Y41 in `tracker/esr1.tasks` (phase `ext`, exact-path checks); keep-all (all 788
matrix chains in every full run, subsets are sensitivity checks only); frozen outputs are read, never rewritten; new
outputs in `results/ext/`, `figures/ext/` (200 dpi PNG + SVG), scripts `scripts/ext_*.py`; seed 42; one palette
(`scripts/ext_common.py`: cl7 blue, cl4 orange, minor clusters aqua + marker). Ligand action labels (Q10/Q11) are
labelled by hand by a teammate; nothing here drafts them. The only action labels shown are the 7 evidence-supported
block-07 labels; everything else is "unlabelled".

The structural reference partition is S08 (average linkage, 1.5 A cut): two main clusters, cl7 (445 chains) and cl4
(329 chains), plus five minor clusters (14 chains).

## Pre-declared selection rule (structural and chemistry zoo)

Declared before any run: every run with >= 2 clusters and noise fraction <= 0.10 is scored by
`score = mean(rank(silhouette_precomputed), rank(stability_ARI))` (higher better for both). Davies-Bouldin and
Calinski-Harabasz are reported but not used in the score (both assume convex Euclidean clusters, which an RMSD space
need not have). Ties are broken by fewer clusters. The winner is reported together with the whole ranking table.
For chemistry, the same rule is applied over fingerprint x method x K, and the result is marked
`selection = auto-selected` (C08 still needs a human sign-off).

## Tasks, priority, expected figures

| id | priority | what | figures (figures/ext/) |
|---|---|---|---|
| Y01 | infra | style, palette, card registry, `pymol_safe.py` | - |
| Y02 | P1-E | mmCIF experimental metadata (pH, temperature, method, solution/additives, space group, cell, ASU, year, expression, mutations, peptide) | Y02_metadata_overview |
| Y03 | P1-F | shared embeddings: PCA, classical MDS, metric MDS, t-SNE, UMAP + trustworthiness | Y03_embeddings |
| Y04 | P1-A | clustering zoo (13 families, ~60 settings) + metrics + bootstrap stability + ranking | Y04_zoo_metrics, Y04_zoo_grid |
| Y05 | P1-A | compactness per cluster (mean/max intra, medoid radius, diameter, cohesion/separation) | Y05_compactness |
| Y06 | P1-A | method agreement ARI/AMI, contingency, alluvial | Y06_agreement, Y06_alluvial |
| Y07 | P1-A/F | dendrogram with metadata bars, ordered RMSD heatmap with sidebars | Y07_dendrogram_metadata, Y07_rmsd_heatmap_annotated |
| Y08 | P1-B | chemistry zoo: 6 fingerprints x 7 methods x K, auto-selection (C08) | Y08_chem_zoo_metrics, Y08_chem_map |
| Y09 | P1-B | Murcko scaffolds + physchem | Y09_physchem_by_cluster, Y09_scaffolds |
| Y10 | P1-C | per-residue displacement cl7 vs cl4 with helix bands | Y10_displacement_profile |
| Y11 | P1-C | H12 geometry, helix angles, charge clamp, contact maps, pocket/ligand size | Y11_features_by_cluster, Y11_contact_map_diff |
| Y12 | P1-C | phi/psi PCA | Y12_phipsi_pca |
| Y13 | P1-C | supervised prediction of the split, permutation importance, MI | Y13_feature_importance, Y13_decision_tree |
| Y14 | P1-C | hypothesis tests (peptide, mutation, ligand size, H12, quality, crystal form) | Y14_split_hypotheses |
| Y15 | P1-D | outliers: kNN, medoid distance, LOF, isolation forest, Mahalanobis | Y15_outliers |
| Y16 | P1-D | outlier vs quality (MWU/Fisher/chi2, effect sizes, BH) | Y16_outlier_quality |
| Y17 | P1-E | metadata vs cluster / outlier (Cramer's V, KW, MI, BH) | Y17_metadata_effects |
| Y18 | P1-F | colour-by small multiples (quality / experiment / biology) | Y18_colour_by_* (3) |
| Y19 | P1-F | interactive plotly embedding | docs/interactive_embedding.html |
| Y20 | P2-G | PCA of superposed C-alpha coordinates | Y20_coord_pca, Y20_coord_pca_loadings |
| Y21 | P2-H | B-factor by helix / group / cluster / tier, core selection (K03 with core v3) | Y21_helix_b, Y21_core_selection |
| Y22 | P2-I | ESR2 parallel analysis | Y22_esr2_parallel |
| Y23 | P2-I | joint ESR1+ESR2 matrix + embedding, per-residue receptor differences | Y23_joint_embedding, Y23_receptor_residue_diff |
| Y24 | P2-I | ESR1 vs ESR2 potency / selectivity (ChEMBL) | Y24_potency_esr1_vs_esr2 |
| Y25 | P2-J | ligand <-> conformation link on the 788 matrix (P04 redone) | Y25_chem_struct_link, Y25_nn_agreement |
| Y26 | P2-K | rarefaction, Hill numbers, set cover (explicit "90% of diversity"), representatives | Y26_rarefaction, Y26_hill_numbers, Y26_set_cover |
| Y27 | P2-L | time trends, novelty | Y27_timeline, Y27_novelty |
| Y28 | P3-M | robustness grid method x core x subset | Y28_robustness_grid |
| Y29 | Phase A | ESR2 reference core-fraction question | (table) |
| Y30 | figures | PyMOL renders via pymol_safe | Y30_pymol_medoids, Y30_pymol_bputty |
| Y31 | Phase A | upstream rebuild scripts + frozen-vs-rebuilt report | - |
| Y32 | Phase A | human-step preparation page | - |
| Y33-Y41 | Phase C | index + gallery, exploration log, numbers, METHODOLOGY, PIPELINE/REPO_MAP, RESULTS_SUMMARY, run_all + lock, evidence, final pack | - |

Extra ideas are added as further Y tasks when they come up (see docs/EXPLORATION.md).

## Phase A items and where they are handled

- F01-F04 PyMOL: F01/F02 are DONE; F03 needs verified selectivity positions (E02 <- E01, human); F04 needs F03.
  New PyMOL renders go through `scripts/pymol_safe.py` (Y30).
- P04 chemistry x structure on the new 788 x 788 matrix: Y25 (L02/L03 stay blocked on C08 sign-off).
- C08 choice of K: auto-selected in Y08 by the declared rule; C08 stays MANUAL until a human signs.
- S15 "90% of diversity": Y26 gives the explicit set-cover definition (fraction of chains within radius r of a chosen
  representative) next to the S15 cluster-count definition.
- ESR2 reference core fraction: Y29.
- K03 core B with the new core: Y21 (`core_b_v3.csv`).
- Missing upstream scripts: Y31.
- Human tasks: Y32 (`docs/HUMAN_STEPS.md`).
