# Exploration log (final run, from 2026-10-25)

Running log of everything tried in the extended analyses, including what showed nothing. One entry per attempt:
what, result (with the file it is in), and whether it led anywhere. Numbers are copied from the named result files
and the tracker log; the documents in docs/ read them through scripts/collect_numbers_ext.py.

## Data and metadata

- **mmCIF metadata (Y02).** Parsed 351 local mmCIFs with gemmi; RCSB API only for gaps (cached in
  data/cache/rcsb_entry/). pH is missing from `_exptl_crystal_grow.pH` in 109 entries; 65 of them state the pH in the
  free-text crystallisation details and were parsed (marked `ph_source`), 44 remain not_reported. Almost all crystals
  are PEG-grown (323/351), so "solvent" has little variance to explain anything. File: results/ext/entry_metadata.csv.
- **Construct mutations are the striking variable.** The residue-level mapping (results/lbd_chains.csv `mismatches`)
  shows 443 chains carry Y537S and 243 carry L536S (often with L372S or the C381S/C417S/C530S triple). Cross-tab with
  the S08 1.5 A clusters: Y537S chains are 405/445 of cl7, L536S chains are 237/329 of cl4. A co-crystallised
  peptide is present in 414/445 cl7 chains (378 of them NCOA2/SRC-2/GRIP1 or NCOA1/SRC-1) and in 4/329 cl4 chains
  (none a coactivator). -> led to Y14
  hypothesis tests. File: results/ext/chain_metadata.csv.

## Embeddings (Y03)

- Five maps of the same matrix (PCA of rows, classical MDS, metric MDS, t-SNE, UMAP). All separate cl7 from cl4.
  PCA of RMSD rows puts 90% of variance on PC1. t-SNE / UMAP have the best local trustworthiness (0.94 / 0.91) but
  distort global distances (Spearman 0.68 / 0.58); metric MDS has the best global agreement (0.92) and became the
  shared map. File: results/ext/embedding_quality.csv.

## Clustering zoo (Y04)

- 81 runs in 14 families. The declared rule (silhouette + bootstrap stability ranks) selects spectral clustering
  (sigma 2 A, K = 2); the top runs are all K = 2-3 and agree with S08 (ARI 0.93-0.98).
- **Did not work:** DBSCAN (eps >= 0.75 A gives one cluster; smaller eps gives many tiny clusters + noise), HDBSCAN
  and OPTICS (most chains noise or one blob), affinity propagation (10 clusters, negative silhouette), single
  linkage (chaining). Density-based methods fail because a few bridging chains connect the two clusters.
- **PAM surprise:** k-medoids (K = 2) picks two chains with no modelled H12 (5T92|A|1, 2QE4|B|1) as medoids and gives
  a partition that mixes cl7 and cl4 (silhouette 0.21). Its total cost (328) is lower than with the S08 medoids
  (359), so this is not a bug: it exposed a property of the matrix (next entry).

## Is the RMSD a metric? (Y42)

- S04 RMSD uses the positions observed in both chains of each pair. 5.4% of random triples violate the triangle
  inequality (max excess 6.5 A). Chains with few observed positions (no H12) have a mean RMSD of ~0.5 A to everyone:
  they are hubs. Average linkage is robust to this; medoid- and density-based methods are not.
- Same superposition, other position sets: on the 181 positions observed in >= 99% of chains, without 528-550, or
  on the core only, cl7 and cl4 overlap (median cl7-cl4 0.44-0.82 A, one cluster at 1.5 A). On H12 alone (538-547)
  the median cl7-cl4 distance is 14.0 A. **The two clusters are two H12 positions.** File: results/ext/metric_check.csv.

## Compactness and agreement (Y05, Y06, Y07)

- cl7 is tighter (mean pairwise RMSD 0.70 A, medoid radius 0.44 A) than cl4 (0.91 / 0.50 A); both far below the
  3.5 A between them. cl4 has visible sub-blocks in the ordered heatmap (Y07).
- The S12 medoid of cl4 (2QE4|B|1) has no modelled H12 - a hub chain (Y42). For pictures and profiles a
  "complete-H12 medoid" is used instead (cl4: 7UJ7|C|1).
- Top-7 families agree with S08 at ARI 0.89-0.98; density, medoid, single-linkage and affinity propagation do not
  (ARI <= 0.16). Alluvial: ribbons stay together through spectral / k-means / average / GMM.

## Chemistry (Y08, Y09)

- 276 chemistry runs. The declared rule picked a degenerate partition (MACCS average linkage K=2: 305 + 1 ligands).
  Lesson: silhouette + stability alone reward "one outlier vs the rest". A post-hoc constraint (largest cluster
  <= 60%) gives MACCS k-means on MDS, K=6 (silhouette 0.29, stability 0.97) -> C08 auto-selected, human sign-off.
- Silhouettes in chemistry space stay low (<= 0.29 for non-degenerate runs): chemical series overlap.
- Unit = ligand: ligands found only in cl4 chains are bigger (median 34 vs 24 heavy atoms), more flexible and carry
  a basic amine in 74% vs 1% (Fisher q 7e-44). TPSA shows nothing (q 0.26). 7 of the 12 most common Murcko scaffolds
  occur in one structural cluster only.

## Why two clusters (Y10-Y14)

- Per-residue profile: positions 528-550 carry 99% of the squared cl7-cl4 mean-position difference; outside them the
  largest shifts are the H2-H3 loop (~3.6 A at 333) and the H8-H9 loop (~2.2 A at 467), both poorly observed loops.
- Geometry: no overlap at all between clusters in H12 position (shift 0.2 vs 12.3 A from the cl7 mean), H12-H3 axis
  angle (82 vs 43 deg), Lys362-Glu542 C-alpha distance (19.3 vs 13.6 A). This is the textbook agonist (H12 capping
  the pocket, AF-2 groove free for a coactivator) vs antagonist (H12 lying in the AF-2 groove) geometry.
- phi/psi PCA (no superposition at all): 81% of dPC1 is one residue, 537 (psi +157 deg in cl7, -47 deg in cl4), the
  H11-H12 hinge - and the site of the Y537S construct mutation that stabilises the agonist hinge.
- Supervised: quality alone predicts the cluster with balanced accuracy 0.80, ligand chemistry 0.89, construct
  (peptide + mutation) 0.96, crystal form 0.97. The tree splits first on peptide.
- Hypothesis tests: peptide (entry log OR 6.7), Y537S (4.2), L536S (-6.7; no cl7 entry has L536S), P21 crystal form
  (6.1). cl4 structures have *better* quality (resolution 1.87 vs 2.22 A). Quality adds McFadden R2 0.05 once
  peptide and mutation are known. D538G and wild-type status show nothing.
- Interpretation (hypothesis, consistent with all tests): the two clusters are two H12 states, and which state a
  crystal shows is set by the experiment design - agonist-state crystals are grown with a coactivator peptide and
  often the Y537S mutation (P21 form); antagonist/SERM crystals use L536S constructs, no peptide, and ligands with a
  basic side chain. Ligand action labels (Q10, teammate) will allow the direct agonist/antagonist test.

## Outliers, quality and metadata (Y15-Y17, Y44)

- Five detectors agree poorly (Jaccard 0.02-0.31). Robust Mahalanobis (MCD, chi-square p < 0.001) flags 339 of 788
  chains: the clusters are not Gaussian, so this detector alone is not useful; it is kept but only counts in the
  >= 3-method consensus (22 chains, 17 entries).
- Answer to "quality or conformation?": consensus outliers do not differ from the rest in resolution, R-free, B or
  missing residues (q >= 0.24); 6 of 22 have an individual quality reason (resolution / R-free beyond the clean 95th
  percentile, or no H12), 16 look like real conformations. They differ in conditions: fewer peptides, Y537S, P21.
- Metadata vs cluster: crystal form (space group V 0.73), oligomer state, construct and peptide dominate.
  Precipitant class, PEG, glycerol, DMSO and bound additives show nothing at entry level - "solvent" does not explain
  the clusters (nearly every crystal is PEG-grown, 323/351).
- cl4 structures have **better** resolution (median 1.87 vs 2.22 A) and lower B than cl7.
- Crystal contacts at H12 (gemmi with symmetry): cl7 H12 is touched by the coactivator peptide, hardly by the lattice
  (median 1 atom); cl4 H12 touches lattice/ASU neighbours (median 10 atoms), and within cl4 more contacts go with a
  complete H12 (17.5 vs 8 atoms) and a larger H12 shift (Spearman 0.55). Packing helps order H12 in the
  antagonist-like state.

## Coordinate PCA, B-factors (Y20, Y21)

- Coordinate PCA: mode 1 = 86.5% of variance and is the H12 switch; modes 2 and 3 are also H12 motions.
- Helix stability ranking by normalised B is the same in every B-group, cluster and tier (H8 > H4/5 > H10/11 ~ H3 >
  H7 > ... > H12 > H6). Core v3 per-chain normalised B is lower than the old core (0.93 vs 1.04). H12/core-v3 B ratio:
  1.33 (cl7) vs 1.53 (cl4).

## ESR2 (Y22-Y24, Y29)

- ESR2 (86 chains) shows the same split: main basin 93% peptide-bound, second basin 0% peptide.
- Joint ESR1+ESR2 space (ESR2 fitted on the ESR1 reference): same-state cross-receptor distance 1.3 A (agonist) and
  2.3 A (second basin) vs 3.5 A between ESR1 states: the H12 state matters more than the receptor.
- The two draft pocket selectivity residues (L384M, M421I) move < 0.2 A in backbone: selectivity there is a
  side-chain effect. Backbone differences between receptors sit in the H8-H9 region (458-500), where the sequence
  differs most (E10 block).
- ChEMBL: 1615 compounds with potency at both receptors; Spearman 0.50; 54% within 10-fold; PDB ligands of cl4 chains
  lean ESR1-selective (median +1.16 vs +0.08 log units, p 0.037).
- ESR2 core fraction 1.0 is unreachable because ESR1 position 466 has no ESR2 counterpart (deletion). Changing the
  ESR2 reference (3 candidates or the ESR1 frame) gives identical clusters (ARI 1.00).

## Ligand and conformation (Y25)

- P04 redone on 788 chains with the auto chemistry clusters: ligand-level permutation p 1e-4, Cramer's V 0.43
  (CI 0.40-0.55). 89% of ligands share the H12 state of their nearest chemical neighbour (null 51%).
- Potency: cl4 ligands are more potent at ESR1 (+1.13 pChEMBL raw; +0.83 within chemistry clusters, p 0.02).

## Redundancy (Y26, Y27, Y43)

- Hill numbers: whatever the method, core or subset, q1 = 1.7-3.3 effective conformations (S08 1.5 A: q0 7, q1 2.19).
- Set cover with complete-H12 representatives: 2 structures place 90% of all chains within 1.0 A; 48 needed at 0.5 A.
  Allowing any chain as representative gives 1 at 1.0 A - the hub effect again (shown, not used).
- Time: since 2010, 3 of 673 new chains are > 1.5 A from every earlier complete-H12 chain; median novelty 0.38 A.
- Copies in one crystal: median 0.54 A, 98% in the same cluster; the same ligand in another entry 0.65 A (but 15% of
  same-ligand pairs sit in different clusters).

## Robustness (Y28) and cl4 substructure (Y45)

- 6 methods x 2 cores x 5 subsets: split agreement 0.98-1.00 everywhere. Negative controls without H12 drop to
  0.53-0.97: k-means / GMM still partly recover the split from the loops (secondary signature).
- Inside cl4 the 1.0 A sub-partition separates L536S constructs (cl4-1: 84% L536S, P1) from non-L536S ones (cl4-2:
  0% L536S, P6522); they differ at the H11-H12 loop 533-534 (where residue 536 sits) and the H2-H3 loop, not in H12.
  A silhouette-chosen average-linkage split of cl4 was degenerate (327 + 2) - recorded, not used.

## Things that showed nothing

- TPSA of ligands (q 0.26), D538G and wild-type status (n.s.), precipitant/additive class, X-ray source, diffraction
  temperature, completeness, Wilson B vs cluster; release year vs novelty (flat).
