# Methodology: ESR1 structural diversity and redundancy

This document describes the whole pipeline in the order it runs, from the PDB query to the extended analyses of the
final run. For every step: purpose, inputs (file and format), outputs, the code (script and function), the
parameters and why, and what is recorded for reproducibility. Tracker task IDs (Q.., K.., S.., R.., C.., E.., X..,
Y..) point to `tracker/esr1.tasks` (the spec) and `tracker/esr1.log` (what was done, with evidence paths).
Numbers in this file are filled in by `scripts/collect_numbers_ext.py` from `results/ext/numbers_ext.json`.

Conventions used everywhere:

- **Keep-all.** Every chain that passes the completeness screen ({{R.n_chains}} chains from {{R.n_entries}} entries)
  is in every full analysis. Quality is a label (tier A/B/C, H12 label, outlier flags), never a filter. Subsets
  (tier A only, one chain per entry, complete H12) are sensitivity runs next to the full run.
- **Frozen inputs are loaded, never recomputed** (`results/rmsd_entry.npy`, `results/tanimoto_*.npy`). New outputs of
  the final run go to `results/ext/` and `figures/ext/` (200 dpi PNG + SVG); scripts are `scripts/ext_*.py`.
- **Seed** {{R.set_random_seed}} for every random step. **Units:** distances in Angstrom (A), B-factors in A^2.
- **Numbering:** UniProt P03372 (ESR1). ESR2 residues carry ESR1 numbers through the E03 numbering map.
- **Settings:** `esr_analysis/resources/config/analysis_settings.json` (thresholds) and
  `esr_analysis/resources/config/alignment_core_v3.json` (alignment core).

## 1. Data acquisition (Q01-Q03, X13, W04, Y31)

- **Purpose:** collect every X-ray structure of the ESR1 (P03372) and ESR2 (Q92731) ligand-binding domain.
- **Code:** `esr_analysis/query_structures.py` (RCSB search API), `esr_analysis/download_structures.py`
  (`download_one`), `scripts/fetch_data.py` (one command to rebuild git-ignored inputs: mmCIFs, CCD files, data/lbd);
  `scripts/rebuild_fetch_raw.py` (final run: fetches any missing mmCIF and checks it against
  `data/raw/sha256sums.txt`).
- **Inputs/outputs:** `data/pdb_ids_esr1.txt`, `data/pdb_ids_esr2.txt` (one PDB ID per line) -> `data/raw/<ID>.cif`
  (mmCIF), `data/raw/download_manifest.csv`, `data/raw/sha256sums.txt`; ESR2 mmCIFs for metadata in
  `data/cache/esr2_cif/`; ChEMBL activities `data/chembl_esr1_activities.csv` (CHEMBL206, IC50/Ki/Kd, relation '=')
  and `results/ext/chembl_esr2_activities.csv` (CHEMBL242, downloaded in the final run, date in
  `results/ext/chembl_esr2_meta.json`).
- **Recorded:** sha256 per mmCIF, download manifest, ChEMBL query URL and date.

## 2. Chain extraction, mapping and completeness screen (Q04, K01, K06)

- **Purpose:** find every receptor chain in every entry, map its residues to canonical ESR1 numbering, and decide
  whether it has enough of the LBD to be measured.
- **Code:** `esr_analysis/extract_lbd.py` (`inventory_entry`), `esr_analysis/map_residues.py`; rebuilt from scratch by
  `scripts/rebuild_lbd_extract.py` (Y31). Attrition: `scripts/attrition_funnel.py`; master status table:
  `scripts/chain_status.py`.
- **Outputs:** `results/lbd_chains.csv` (one row per chain: identity, mapped LBD/core/H12 C-alpha counts, fractions,
  mismatches, eligibility and reason), `results/lbd_extract/{entries,receptor_chains,residue_mapping}.csv`
  (residue_mapping = one row per residue with ESR1 position, author numbering, C-alpha coordinates and B-factor),
  `results/chain_status.csv` (one row per ESR1 chain, {{R.n_chains_all_esr1}} rows, none dropped).
- **Parameters:** LBD = ESR1 {{R.set_esr1_lbd_start}}-{{R.set_esr1_lbd_end}}; a chain is measurable when it has at
  least {{R.set_min_chain_lbd_ca}} mapped LBD C-alpha atoms and {{R.set_min_core_ca_for_alignment}} core C-alpha atoms
  (old core). This is a completeness screen (a chain with nothing to superpose cannot be compared), not a quality
  filter: {{R.n_chains_out_of_matrix}} ESR1 chains fail it, each with a stated reason (no LBD in the construct, or
  below the minimum).
- **Recorded:** the reason for every chain that is not in the matrix (`chain_status.csv`, `chain_loss_reasons.csv`).

## 3. Quality labels (Q06, K05, S16, K03)

- **Purpose:** label, never remove, structures by quality.
- **Code:** `esr_analysis/flag_quality.py` (entry flags), `scripts/qc_tiers.py` (tiers), `scripts/h12_metric.py`
  (H12 completeness and position), `scripts/helix_bfactors.py` (normalised B per helix and B-group).
- **Outputs:** `results/entries.csv`, `results/qc_tiers.csv` (tier per metric and overall qc_label),
  `results/h12_metric.csv` (h12_label complete / partial / not_modelled), `results/helix_bfactor_by_chain.csv`.
- **Parameters:** tier A = the original flags (resolution <= {{R.set_qc_resolution_A}} A, R-free <=
  {{R.set_qc_r_free}}, mean B <= {{R.set_qc_mean_B_A2}} A^2); tier B = relaxed limits in
  `esr_analysis/resources/config/qc_tiers.json`; tier C = beyond. Matrix chains: {{R.tier_A_chains}} tier A,
  {{R.tier_B_chains}} tier B, {{R.tier_C_chains}} tier C. B-factors are normalised per chain (C-alpha B / the chain's
  median C-alpha B) so refinement differences cancel; B-group = tertile of the chain mean B.
- **Recorded:** a plain-language warning per entry; every threshold in config files.

## 4. Alignment core (K02, K04; Y21 refresh)

- **Purpose:** superpose on the part of the LBD that is most stable, so that real motions are not hidden by a bad fit.
- **Code:** `scripts/helix_table.py` (helix ranges from the mmCIF HELX_P records of reference 7NFB), 
  `scripts/helix_stability.py` (rule and core), `scripts/ext_bfactors.py` (Y21: the same ranking per cluster and tier,
  and per-chain core B on the new core, `results/ext/core_b_v3.csv`).
- **Rule:** helices whose median normalised B ranks in the better half in every B-group, never H12, >= 60 residues.
- **Result:** core v3 = {{R.core_v3_ranges}} ({{R.core_v3_n}} positions), stored in
  `alignment_core_v3.json`; the old core ({{R.old_core_ranges}}) is kept for comparison. Per chain, the median normalised
  B on core v3 is {{Y21.core_v3_median|.3f}} vs {{Y21.old_core_median|.3f}} on the old core.

## 5. C-alpha superposition (S03)

- **Purpose:** put every chain in one frame.
- **Code:** `scripts/align_selected_chains.py` -> `esr_analysis/align_structures.fit_to_reference` (Kabsch on the
  shared core C-alpha atoms).
- **Inputs/outputs:** residue_mapping coordinates -> `results/alignment/aligned_coordinates.npz` (788 x 246 x 3,
  NaN where unobserved), `alignment_transforms.json` (rotation + centroids per chain), `global_alignment.csv`,
  `results/alignment_core_rmsd.csv`; the same for the old core (`*_old_core*`).
- **Parameters:** reference chain {{Y23.reference}} (chosen by the block-08 reference policy, `esr_analysis/choose_reference.py`,
  frozen in `esr_analysis/resources/frozen/reference.json`); minimum
  {{R.set_min_core_ca_for_alignment}} shared core C-alpha atoms. Core RMSD > 2 A is labelled, never dropped.

## 6. RMSD matrix (S04, S05; Y42 check)

- **Purpose:** one distance per pair of chains.
- **Code:** `scripts/build_selected_rmsd.py` -> `esr_analysis/calculate_rmsd.compute_distances`.
- **Definition:** C-alpha RMSD over the LBD positions observed in both chains, in the S03 frame (no refit per pair).
  Pairs need >= {{R.set_min_pairwise_common_ca}} common positions.
- **Outputs:** `results/rmsd_entry.npy` ({{R.n_chains}} x {{R.n_chains}}, float64), `results/rmsd_entry_index.csv`
  (row, chain_key, entry_id).
- **Known property (Y42, `scripts/ext_metric_check.py`):** because the position set changes from pair to pair, the
  matrix is not a metric: {{Y42.full.triangle_violation_frac|.1%}} of random triples violate the triangle
  inequality, and chains without a modelled H12 are close to everyone ("hubs"). Average linkage is robust to this;
  medoid- and density-based methods and nearest-neighbour statistics are not, so the final run uses complete-H12
  chains as medoids, representatives and novelty references.
- **Recorded:** `results/alignment/distance_exclusions.json`, overlap counts (`overlap_entry.npy`).

## 7. Clustering (S08, S12, S14; Y04-Y07, Y28)

- **Reference partition (S08):** `scripts/cluster_selected.py` -> `esr_analysis/cluster_structures.structural_groups`
  (SciPy average linkage, distance criterion) at 1.0 / 1.5 / 2.0 / 2.5 A; primary cut
  {{R.set_primary_structural_threshold_A}} A gives {{R.nclusters_cut1p5}} clusters: cl7 ({{R.cl7_size}} chains),
  cl4 ({{R.cl4_size}}) and five minor clusters ({{R.minor_size}} chains). Outputs `results/entry_clusters.csv`,
  `results/linkage_entry.npy`; medoids `results/medoids.csv` (`scripts/medoids.py`).
- **Clustering zoo (Y04, `scripts/ext_zoo.py`):** {{R.zoo_n_runs}} runs in {{R.zoo_n_families}} families -
  hierarchical average / complete / single on the matrix, Ward on a 10-D classical-MDS embedding, k-means, PAM
  k-medoids (implemented in the script: BUILD + SWAP), spectral (affinity exp(-d^2 / 2 s^2), s = 0.5 / 1 / 2 A),
  DBSCAN, HDBSCAN, OPTICS, affinity propagation, mean shift, BIRCH, Gaussian mixture. Metrics: silhouette on the
  precomputed RMSD, Davies-Bouldin and Calinski-Harabasz on the embedding, bootstrap stability (50 subsamples of 80%,
  ARI), noise fraction. **Selection rule declared before running:** eligible = >= 2 clusters and noise <= 10%; score =
  mean percentile rank of silhouette and stability. Winner: {{Y04.best_run}} (silhouette {{Y04.best_sil|.2f}},
  stability {{Y04.best_stab|.3f}}, ARI vs S08 {{Y04.best_ari_s08|.2f}}). Outputs `results/ext/zoo_runs.csv`,
  `zoo_labels.csv`, `zoo_ranking.csv`.
- **Compactness (Y05):** per cluster mean / max intra-cluster RMSD, medoid radius, diameter, cohesion / separation
  (`scripts/ext_compactness.py`, `results/ext/zoo_compactness.csv`).
- **Agreement (Y06):** ARI / AMI between the best run of every family, contingency and alluvial
  (`scripts/ext_agreement.py`). **Dendrogram and ordered heatmap (Y07):** `scripts/ext_dendrogram.py`.
- **Robustness (Y28):** 6 methods x 2 cores (+ 2 negative-control position sets) x 5 subsets
  (`scripts/ext_robustness.py`, `results/ext/robustness_grid.csv`).

## 8. Embeddings (S11-like, Y03, Y20)

- **Code:** `scripts/ext_embeddings.py`: PCA of RMSD rows, classical MDS, metric SMACOF MDS, t-SNE, UMAP (umap-learn
  0.5), all on the same matrix; trustworthiness (k = 10) and Spearman between RMSD and map distance. **Rule:** the
  map with the highest Spearman (global faithfulness) is the shared map for every colour-by figure: metric MDS.
- **Coordinate PCA (Y20, `scripts/ext_coord_pca.py`):** PCA of the superposed C-alpha coordinates (positions observed
  in >= 95% of chains, mean-imputed); explained-variance ratios of modes 1-5: {{Y20.ev}}.

## 9. Why two clusters: features, tests and models (Y10-Y14, Y42, Y45)

- **Per-residue displacement (Y10):** distance between cl7 and cl4 mean C-alpha positions, spread, observation rate
  (`scripts/ext_displacement_profile.py`).
- **Structural features (Y11):** H12 centroid shift, H12-H3 / H12-H11 axis angles (first principal axis of each
  helix), K362-E542 C-alpha distance, L540 to pocket centroid, contact maps (C-alpha < 8 A)
  (`scripts/ext_structural_features.py`).
- **Dihedrals (Y12):** phi/psi from `data/clean/*.pdb` (Bio.PDB PPBuilder), sin/cos encoding, PCA
  (`scripts/ext_phipsi.py`) - an alignment-free check.
- **Supervised check (Y13):** decision tree (depth 3), random forest (300 trees), logistic regression; GroupKFold
  (5 folds, groups = entries, so copies of one crystal never sit in train and test); permutation importance and
  mutual information (`scripts/ext_supervised.py`).
- **Hypothesis tests (Y14):** Fisher exact (log odds ratio, Woolf CI), Mann-Whitney (rank-biserial, entry-bootstrap
  CI), logistic likelihood-ratio tests, Mantel-Haenszel; chain level plus entry / ligand level; Benjamini-Hochberg
  (`scripts/ext_split_hypotheses.py`).
- **Region-restricted matrices (Y42):** fixed positions, no 528-550, H12 only, core only (new files
  `results/ext/rmsd_*.npy`).

## 10. Outliers, quality and experimental metadata (Y02, Y15-Y17, Y44)

- **Metadata (Y02, `scripts/ext_metadata.py`):** gemmi 0.7 parse of each local mmCIF (`_exptl`, `_exptl_crystal_grow`,
  `_symmetry`, `_cell`, `_pdbx_audit_revision_history`, `_entity_src_gen`, `_entity`, `_pdbx_struct_assembly`,
  `_refine`, `_reflns`, `_diffrn`); RCSB Data API only for gaps (cached in `data/cache/rcsb_entry/`); pH parsed from
  the free-text growth conditions when the pH item is empty (marked in `ph_source`); missing = not_reported.
  Construct mutations come from the residue-level mapping (`lbd_chains.csv` mismatches with cause engineered
  mutation / conflict); peptide = any non-receptor polymer of <= {{R.set_peptide_max_length}} residues.
- **Outliers (Y15):** kNN distance and distance to the nearer complete-H12 main medoid (robust z > 3.5), LOF (> 1.5),
  isolation forest, robust Mahalanobis (MCD, p < 0.001); consensus = >= 3 methods (`scripts/ext_outliers.py`).
- **Outlier vs quality (Y16):** Mann-Whitney / Fisher with BH against all chains and against the same main cluster;
  a per-outlier verdict (`scripts/ext_outlier_quality.py`).
- **Metadata tests (Y17):** chi-square + Cramer's V, Kruskal-Wallis + epsilon^2, mutual information, BH per target and
  level (`scripts/ext_metadata_tests.py`).
- **Crystal contacts at H12 (Y44):** gemmi NeighborSearch with symmetry, partner atoms within 4 A of H12 atoms
  (`scripts/ext_h12_contacts.py`).

## 11. Diversity, entropy, saturation and redundancy (R01-R04, S15; Y26, Y27, Y43)

- **Earlier:** Shannon entropy of cluster sizes (`scripts/cluster_entropy.py`), discovery curve over 100 random
  orderings and n_90 (`scripts/diversity_curve.py`, `scripts/saturation.py`): n_90 = {{R.n90_cut1p5}} chains at
  1.5 A; representative set by clusters (`scripts/representative_set.py`).
- **Final run (Y26, `scripts/ext_redundancy.py`):** Hill numbers q0 / q1 / q2 (effective number of conformations);
  greedy set cover at radius r = 0.3-2.0 A for 90 / 95 / 99 / 100% of chains, representatives restricted to
  complete-H12 chains (hubs excluded; the unrestricted result is reported too); rarefaction of cluster discovery and
  of coverage; farthest-first representative sets; redundancy per cluster. This makes "90% of the conformational
  diversity" explicit as coverage at a radius, next to the S15 cluster-count definition.
- **Time (Y27):** novelty = minimum RMSD to complete-H12 chains released in earlier years (`scripts/ext_time.py`).
- **Copies (Y43):** within-entry vs between-entry RMSD (`scripts/ext_ncs_copies.py`).

## 12. Chemistry (Q08, Q16, Q18, C01-C11; Y08, Y09, Y25)

- **Ligand table and ligand of interest:** CCD SMILES (`esr_analysis/build_ligand_table.py`; rebuild
  `scripts/rebuild_ligands.py`); the ligand in each chain's pocket = a heavy atom within {{R.set_contact_cutoff_A}} A of
  a pocket residue (`scripts/ligand_of_interest.py`, `results/chain_ligand_of_interest.csv`).
- **Fingerprints and clustering (Y08, `scripts/ext_chem_zoo.py`):** RDKit 2026.09 generators - ECFP4, ECFP6, FCFP4
  (Morgan radius 2/3, 2048 bits), MACCS, hashed atom pair and topological torsion (2048); Tanimoto distance; Butina,
  PAM, spectral, average linkage, k-means on MDS, HDBSCAN, DBSCAN; {{R.chem_n_runs}} runs; same declared rule as the
  structural zoo. The rule's winner was degenerate ({{Y08.best_run}}, largest cluster {{Y08.best_largest|.0%}}); a
  non-triviality constraint (largest cluster <= 60%) was added after seeing it and is reported as post hoc:
  auto-selected {{Y08.balanced_run}}, K = {{Y08.balanced_k}} (C08 MANUAL until signed).
- **Descriptors (Y09):** MW, cLogP, TPSA, HBD, HBA, rotatable bonds, heavy atoms, rings, Fsp3, basic amine (SMARTS),
  Murcko scaffolds (`scripts/ext_ligand_descriptors.py`). Unit of test = ligand.
- **Chemistry x structure (Y25, `scripts/ext_chem_struct_link.py`):** ligand-level chi-square with 10,000
  label permutations, Cramer's V with bootstrap CI, MI, purity, nearest-neighbour agreement, Mantel-style test of
  Tanimoto vs RMSD, potency OLS with chemistry-cluster fixed effects.
- **ChEMBL (X13-X16, Y24):** per-compound median pChEMBL per receptor; selectivity = pESR1 - pESR2
  (`scripts/ext_selectivity.py`).

## 13. ESR2 (X01-X05, E03-E10; Y22, Y23, Y29)

- **Matrix:** `scripts/esr2_rmsd.py` (86 chains, ESR2 reference 3OLL|A|1, core v3), clusters `scripts/esr2_clusters.py`.
- **Final run:** the zoo, outliers and metadata on ESR2 (`scripts/ext_esr2.py`); ESR2 chains fitted on the ESR1
  reference and a joint 874 x 874 matrix (`scripts/ext_joint_esr1_esr2.py`, `results/ext/rmsd_joint.npy`; its ESR1
  block equals S04 to {{Y23.esr1_block_maxdiff|.0e}} A); the reference core-fraction question
  (`scripts/ext_esr2_core_fraction.py`, `docs/esr2_core_fraction.md`).

## 14. Figures, documents and reproducibility

- **Figures:** every `ext_*` script saves 200 dpi PNG + SVG through `ext_common.save_fig` and registers a card
  (what it shows, a takeaway computed from its result file, source table, script) in `results/ext/figure_cards.csv`;
  `scripts/ext_build_index.py` builds `docs/figures/INDEX.md` and `docs/gallery.html`. One palette
  (`scripts/ext_common.py`, validated with the dataviz palette checker): cl7 blue, cl4 orange, minor clusters aqua.
- **PyMOL:** only through `scripts/pymol_safe.py` (separate process, faulthandler, 150 s timeout, PNG check,
  matplotlib fallback) - `scripts/ext_pymol_renders.py`.
- **Numbers:** `scripts/collect_numbers_ext.py` (called by `scripts/collect_numbers.py`) writes
  `results/ext/numbers_ext.json` and renders the documents in `docs/templates/`.
- **Upstream rebuild (Y31):** `scripts/rebuild_all.py` regenerates the upstream tables into `results/rebuild/` and
  compares them with the frozen files: {{R.rebuild_identical}} identical, {{R.rebuild_equivalent}} equivalent
  (float noise or line endings), {{R.rebuild_differs}} differ (UMAP / t-SNE are stochastic; one documented variant).
- **One command:** `bash run_all.sh` (see `docs/PIPELINE.md`); exact package versions in `requirements-lock.txt`.
- **Tracker:** every unit is a task with exact-path checks; `python tracker/tl.py status` shows
  {{R.tasks_DONE}} DONE, {{R.tasks_MANUAL}} MANUAL, {{R.tasks_TODO}} TODO of {{R.tasks_total}} at the time this file was
  generated.

## Limitations

- C-alpha only: side-chain differences (e.g. the ESR1/ESR2 pocket residues L384M, M421I) are invisible to the RMSD.
- The pairwise-common-position RMSD is not a metric (section 6); results that depend on nearest neighbours use
  complete-H12 references.
- Chains of one entry are not independent; tests are repeated at the entry (or ligand) level and permutation tests
  shuffle whole ligands.
- Ligand action labels are pending (Q10/Q11, teammate); the agonist / antagonist reading of the two clusters rests
  on peptide, construct, ligand chemistry and H12 geometry, not on curated labels.
- ChEMBL potencies mix assay types; potency results are descriptive.
