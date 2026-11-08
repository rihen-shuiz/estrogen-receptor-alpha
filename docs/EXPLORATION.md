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

## Suggested, not done (SV11, 2026-11-08)

- Ligand-pocket site set from 4 A heavy-atom contacts over all residues of the SV00 structures (SV11 used the Q18
  pocket_residues as the run file asks: 4.5 A over 19 fixed positions), to check that the pocket labels do not
  depend on that list.
- In cl7, AF-2 groove pockets occur only in peptide-bound chains (234/414 vs 0/31; results/sv/SV11/site_by_chain.csv):
  compare the empty-groove geometry of peptide-bound and peptide-free cl7 chains.
- A small-cluster-robust test (bias-reduced sandwich or entry bootstrap) for the AF-2 groove best rank / probability
  contrast, which rests on 5 cl4 chains in 4 entries (results/sv/SV11/site_tests.csv).

## Suggested, not done (SV02, 2026-11-08)

- Re-run the SV02 secondary structure with mkdssp >= 4 where it can be installed: pydssp's 'H' merges 3-, 4- and
  5-turn helices, so the helix ends (e.g. the proposed H12 end at 544, results/sv/SV02/h12_window.md) may move when
  3_10 and pi turns are separated.
- cl4 starts H9 at 466 and cl7 at 473 (results/sv/SV02/helix_boundaries.csv), but cl7 observes 462-469 in < 50 % of
  chains: check whether the H8-H9 loop/H9 start differs between clusters on chains that observe it in both.

## Suggested, not done (SV01b, 2026-11-08)

- Run the SV01b step-8 medoid H12 clash test (scripts/sv/SV01b_clash.py, written but gated) even though the
  re-issued verdict is "artifact": the lattice terms are counted on modelled atoms of 528-550, so for the
  completeness outcome they are partly a consequence of the outcome (cl4 complete chains carry 12.3 of their 24.3
  mean symmetry-contact atoms on 546-550, which partial chains rarely model; results/sv/SV01b/lattice_contacts_by_group.csv).
  The clash test places the medoid H12 independently of what was modelled.
- A lattice term that does not depend on the chain's own modelled atoms (e.g. neighbour atoms within 4 A of the
  superposed medoid H12, or contacts of 528-545 only), refitted as a sensitivity of the SV01b completeness verdict.
- The mobility models without the lattice terms, to show whether the mobility "artifact" verdict rests on them
  (lattice_sym has the smallest q of the X predictors at chain level in cl7 and cl4; results/sv/SV01b/odds_ratios.csv).

## Suggested, not done (SV03, 2026-11-08)

- H12 axis descriptors on one fixed residue set for every chain (or stratified by the number of C-alpha used): on
  ideal helices the [DEFAULT] window-centroid axis is off by 9.6 / 4.5 / 0.6 deg with 6 / 7 / 8 C-alpha
  (results/sv/SV03/axis_method_check.csv), and cl4 H12 axes rest mostly on 7 C-alpha, cl7 on 8
  (results/sv/SV03/numbers_geometry.json h12_axis_n_ca_by_group), so within-cluster H12 angle differences of a few
  degrees may partly follow how many H12 residues are helical.
- A small-cluster-robust test (entry bootstrap or bias-reduced sandwich) for the 4 SV03 [DEFAULT] hits whose smaller
  group has <= 5 entries (D538G in cl7: 2 entries; peptide-bound chains in cl4: 3 entries; other mutation and H9
  length in cl7: 5 entries; results/sv/SV03/sensitivity_check.csv).
- An interval-censored length comparison that also uses the censored helices: censoring differs by cluster (H3 0.094
  cl7 vs 0.447 cl4, H10/11 0.079 vs 0.413, H12 0.140 vs 0.589; results/sv/SV03/helix_length_summary.csv), so the
  uncensored subsets are selected differently in the two clusters.
- One model per within-cluster hit with mutation, ligand chemistry and peptide together (descriptor rank ~ Y537S /
  L536S + amine + acid + heavy atoms + peptide, entry-clustered): the SV03 rank tests take one covariate at a time,
  and mutation and ligand series co-vary by entry.

## Suggested, not done (SV04, 2026-11-08)

- Midpoint clustering without H12 (the 7 other consensus helices only), to see whether they alone separate cl7 from
  cl4: in SV04 the cl7-cl4 midpoint distance is 12.6 A for H12 and <= 0.30 A for every other helix
  (results/sv/SV04/distance_summary.csv), so the ARI with H12 mostly restates the H12 split, and all 9 K = 2
  mismatches are chains without a usable H12 (results/sv/SV04/numbers_cluster.json k2_chain_mismatch).
- H12 midpoints on the SV02 data-derived window 538-544 as a sensitivity: with 538-547 only 75 chains have every
  residue helical (results/sv/SV04/numbers_midpoints.json complete_helices), so most H12 midpoints extrapolate over
  546-547, which are observed but not helical in most chains; the truncation curve on those 75 fully helical H12s
  may be optimistic for real disordered ends.

## Suggested, not done (SV05, 2026-11-08)

- Repeat any result that uses the S08 1.0 A sub-clusters of cl4 with the core-v4 labels: at 1.0 A the cl4 sub-clusters
  regroup (S08 262/48/11 chains vs v4 211/114; ARI 0.817, 66 chains change; results/sv/SV05/ari.csv,
  cluster_sizes.csv), while at 1.5-2.5 A only 2-3 chains change, so the 1.0 A sub-structure of cl4 depends on the core.
- A look at the three chains whose 1.5 A cluster depends on the core (2JF9|A|1, 2JFA|B|1, 3DT3|B|1;
  results/sv/SV05/SUMMARY.md step 5): which positions carry their v4-v3 RMSD change.

## Suggested, not done (SV07, 2026-11-08)

- Within-cluster look at the two genuine spike regions of the Y10 profile (results/sv/SV07/covariate_tests.csv): at
  335-337 the cl4 chains with a basic-amine ligand sit 3.3-3.8 A along the cl7->cl4 displacement and the cl4 chains
  with 'other' ligands 0.4-0.9 A from the cl7 mean (335: 3.75 vs 0.44 A), and at 467 the 22 peptide-free cl7 chains sit near
  cl4 (467: 1.86 vs 2.19 A). A within-cl4 (chemotype) and within-cl7 (peptide) comparison with matched crystal forms
  would show whether these loops follow the ligand / peptide rather than the cluster.
- The re-picked medoid pair (results/sv/SV07/restricted_profile.csv) exceeds 1.0 A at 463, 465, 470 and 471, where the
  Y10 medoids do not observe the residue; a medoid-free loop comparison (e.g. per-entry distributions) would cover them.

## Suggested, not done (SV01c, 2026-11-08)

- A lattice term that is not counted on modelled atoms: the step-8 medoid H12 clash counts (5TM1|A|1 and 7UJ7|C|1
  positions) computed for every chain, complete ones included, as a block-X predictor. On 538-544 the two SV01b lattice
  terms are dropped by joint separation in cl7 because the not-modelled chains have no atoms on 528-550 to contact
  (results/sv/SV01c/odds_ratios.csv, lattice_contacts_by_group.csv), so the current terms cannot say whether the
  lattice holds H12 in place.
- Completeness on the cluster-specific SV02 windows (cl7 538-545, cl4 537-544; results/sv/SV02/helix_boundaries.csv):
  SV01c used the global consensus 538-544 for both clusters.

## Suggested, not done (SV12, 2026-11-08)

- Validate route (c) against the teammate labels once data/sv/action_class_teammate.csv exists: route (c) labels 24
  entry x ligand rows (SERD 19, antagonist 4, agonist 1; 22 of them have no ChEMBL mechanism record), but only 2 rows
  overlap with route (b), so the [DEFAULT] gate (>= 80 % of >= 10 rows) cannot pass on ChEMBL mechanism records alone
  (results/sv/SV12/route_c_validation.csv, numbers.json).
- Cross-check the InChIKey -> ChEMBL match with the ChEMBL/DrugBank links in the RCSB chemcomp records
  (rcsb_chem_comp_related, cached in data/sv/rcsb_chemcomp/); 146 of 306 CCD codes have no ChEMBL molecule by full or
  connectivity InChIKey (results/sv/SV12/ligand_ids.csv).

## Suggested, not done (SV14, 2026-11-08)

- Pocket volume with a larger grid box or a buriedness-bounded region: at >= 10 rays the pocket reaches the [DEFAULT]
  box (ligand-pocket residues + 6 A) in 203 of 788 chains (612 at >= 8 rays, 4 at >= 12;
  results/sv/SV14/numbers.json), so those volumes are bounded by the box.
- A pocket definition that separates the ligand pocket from the buried region past E353/R394: in 9 of 17 wild-type
  EST reference chains the >= 10-ray pocket joins it (405-495 A^3 farther than 4 A from E2), in 8 it does not
  (<= 21 A^3), which makes the E2 pocket reference bimodal (results/sv/SV14/e2_reference_pocket_lining.csv).
- Apo pocket baseline from more apo structures: only 3 entries pass the SV14 apo check (2B23, 4Q13, 9VTU); 7Y8F and
  7YMK carry a PEG in the pocket (results/sv/SV14/apo_blocking_residues.csv).

## Suggested, not done (SV15, 2026-11-08)

- Bridge relative affinities through a second shared reference (e.g. raloxifene or 4-hydroxytamoxifen measured in
  the same assay as the ligand and, elsewhere, together with E2): 48 primary-set ligands have only an absolute
  binding pChEMBL because none of their binding assays also measured E2, against 14 with a same-assay E2 pair
  (results/sv/SV15/numbers.json, affinity_per_ligand.csv).
- Check which activities the ChEMBL potential_duplicate flag removes: it drops 142 of 551 activities with pChEMBL
  (results/sv/SV15/numbers.json cleaning_cascade); where the flagged copy is the one measured next to E2, keeping the
  original instead could add same-assay pairs.

## Suggested, not done (SV13, 2026-11-08)

- Class representative among chains with a complete H12: the [DEFAULT] (lowest median RMSD in results/rmsd_entry.npy)
  picks 2BJ4|A|1 for SERM, which has no modelled H12 (the matrix compares chains on their common C-alpha positions,
  so chains without H12 are compared on fewer positions); runner-up 2JFA|A|1 (results/sv/SV13/representatives.csv).
- Agonist vs SERM is the only class pair with >= 5 entries per side, but it is not one of the SV13 comparisons and it
  almost coincides with cl7 vs cl4 (agonist 25 of 27 entries in cl7, SERM 19 of 20 in cl4;
  results/sv/SV13/class_by_cluster.csv), which SV03 already tests.
- Re-run SV13 Part B with the teammate labels: SERD needs one more entry (4 of 5) for C1 SERD vs SERM, and C2b needs
  >= 5 agonist entries with a same-assay dpX (3 now: DES 2, J3Z 1; results/sv/SV13/gate.csv).

## Suggested, not done (SV09, 2026-11-08)

- docs/RESULTS_SUMMARY.md is re-rendered from docs/templates/RESULTS_SUMMARY.md by scripts/collect_numbers_ext.py, and
  a re-render drops the appended section "Supervisor requests 2026-10-28". Adding a hook in that script that re-runs
  scripts/sv/SV09_docs.py (or moving the section into the template) would keep it; until then, re-run
  scripts/sv/SV09_docs.py and scripts/sv/SV09_number_check.py after any re-render.

## Suggested, not done (VIS1, 2026-11-08)

- Split the multi images VIS1 left uncut (results/sv/VIS1/cut_log.csv, status not_cut) with method (b) once a
  faithful SVG renderer is available (libcairo for cairosvg or rsvg-convert, plus the DejaVu Sans font, since the
  SVGs carry live text), or re-save those figures with a little more space between axes: in each of them a title or
  tick label of one graph reaches past the edge of its neighbour, so no rectangle separates them.
- scripts/sv/common/sv_fig.py crops each panel to the axes' tight bounding box, which takes in fragments of a close
  neighbour or cuts a label at the panel edge (SV01a/SV01b forests, SV04, SV05, SV06 outliers map, SV14 filling;
  results/sv/VIS1/SUMMARY.md); hiding the other axes before saving each panel would avoid both.

## Suggested, not done (VIS2, 2026-11-08)

- figures/ext/Y02_metadata_overview.png, panel 6: the title says "Peptide in 60% of entries" but the bars show ESR1
  chains per asymmetric unit, so the data claim "210 carry a co-crystallised peptide" has no graph in visuals/
  (results/sv/VIS2/selection.csv); a peptide bar in that figure would give it one.
- RESULTS_SUMMARY section 6 cites figures/ext/Y26_rarefaction.png for "488 random chains to see 90% of the 1.5 A
  clusters", but the number is printed only in figures/diversity_curve.png (used for that claim in visuals/); citing
  that file, or labelling N=488 in Y26, would align the citation with the graph.
- figures/sv/SV13/SV13_h12_class_panels.png draws three empty panels (antagonist, partial agonist, other: 0 chains);
  when SV13 is re-run with the teammate labels, drawing only classes that have chains would make the figure denser.
