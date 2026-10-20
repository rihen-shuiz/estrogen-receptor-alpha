# Evidence for tasks logged DONE after 2026-10-18T23:56+05:00

Repo: local clone, branch tracker-tasks (Windows; Python `.venv/Scripts/python`). Generated from the checks in `scripts/check_evidence.py`; every value below was produced by running those checks. Recount commands use the helpers defined at the top of that script (`csv`, `npy`, `txt`, `up`, `nonwhite`, `nclust`, `medoid`, `core_positions`, `CUTS`). Floats in a claim match when they agree to the decimals written in the claim. Output first rows are cut at 90 characters per line.

Tasks in scope: 34 (latest state DONE, DONE line after the cutoff).

| id | plain name | log time | output path(s) | rows or shape | tl.py check | recount |
|---|---|---|---|---|---|---|
| S06 | Note that distances use only C-alpha atoms | 2026-10-20T22:07+05:00 | docs/limitations.md | 3 lines | PASS | MATCH |
| S08 | Group all chains by structural similarity at four distance cuts | 2026-10-20T22:07+05:00 | results/entry_clusters.csv, results/linkage_entry.npy | 788 rows; shape (787, 4) float64 | PASS | MATCH |
| S09 | Table of cluster counts and sizes per cut | 2026-10-20T22:14+05:00 | results/cluster_stats.csv | 4 rows | PASS | MATCH |
| S10 | Tree drawing (dendrogram) of the clustering | 2026-10-20T22:15+05:00 | figures/dendrogram.png | 1272x725 px | PASS | MATCH |
| S12 | Most central structure of each cluster | 2026-10-20T22:15+05:00 | results/medoids.csv | 7 rows | PASS | MATCH |
| S13 | Does keeping only good-quality structures change the clusters? | 2026-10-20T22:21+05:00 | results/qc_clustering_delta.csv | 12 rows | PASS | MATCH |
| S14 | Plot of cluster count against distance cut | 2026-10-20T22:23+05:00 | figures/structural_sensitivity.png | 1480x654 px | PASS | MATCH |
| S07 | Compare C-alpha distances with all-atom distances on 50 entries | 2026-10-20T22:14+05:00 | results/rmsd_allatom_subset.csv, results/rmsd_allatom_subset_summary.json | 1225 rows; 17 top-level keys | PASS | MATCH |
| E05 | Superimpose each matched ESR2 chain on its ESR1 partner | 2026-10-20T22:16+05:00 | results/pair_superposition.csv | 671 rows | PASS | MATCH |
| R01 | How evenly the chains spread over clusters (entropy) | 2026-10-20T22:16+05:00 | results/entropy.csv | 4 rows | PASS | MATCH |
| R02 | How many clusters appear as structures are added at random | 2026-10-20T22:17+05:00 | results/diversity_curve.csv | 3152 rows | PASS | MATCH |
| E06 | Per-residue shift between ESR1 and ESR2 for the first 5 pairs | 2026-10-20T22:18+05:00 | results/first_five_displacements.csv | 1125 rows | PASS | MISMATCH |
| E07 | Per-residue shift for all matched pairs (position x pair table) | 2026-10-20T22:18+05:00 | results/displacement_matrix.csv, results/displacement_pairs.csv | 246 rows; 135 rows | PASS | MATCH |
| R04 | Number of structures needed to see 90% of clusters | 2026-10-20T22:18+05:00 | results/saturation.csv | 4 rows | PASS | MATCH |
| R03 | Figure of the cumulative diversity curve with a 90% line | 2026-10-20T22:18+05:00 | figures/diversity_curve.png | 1200x720 px | PASS | MATCH |
| Q14 | Are rare structural groups worse in resolution or R-free? | 2026-10-20T22:20+05:00 | results/qc_outliers.csv, docs/qc_outliers.md | 788 rows; 39 lines | PASS | MATCH |
| Q15 | Do high B-factor structures gather in particular clusters? | 2026-10-20T22:20+05:00 | results/qc_bfactor_clusters.csv, results/qc_bfactor_clusters_tests.json | 7 rows; 7 top-level keys | PASS | MATCH |
| K06 | One status row for every ESR1 chain, none dropped | 2026-10-20T22:21+05:00 | results/chain_status.csv | 901 rows | PASS | MATCH |
| K09 | Does the choice of alignment core change the clusters? | 2026-10-20T22:23+05:00 | results/core_sensitivity.csv, results/rmsd_old_core.npy | 2 rows; shape (788, 788) float64 | PASS | MATCH |
| E08 | Heatmap of the ESR1/ESR2 per-residue shifts | 2026-10-20T22:24+05:00 | figures/displacement_heatmap.png | 1472x899 px | PASS | MATCH |
| E10 | Positions that shift by more than 1.5 A in most pairs | 2026-10-20T22:24+05:00 | results/consistent_positions.csv | 246 rows | PASS | MATCH |
| R06 | Entropy with and without structures failing R-free | 2026-10-20T22:25+05:00 | results/entropy_qc_compare.csv | 12 rows | PASS | MATCH |
| K07 | One-chain-per-entry view next to the all-chain view | 2026-10-20T22:25+05:00 | results/entry_level_clusters.csv, results/entry_vs_chain_clusters.csv | 351 rows; 4 rows | PASS | MATCH |
| S15 | Smallest set of representative structures | 2026-10-20T22:26+05:00 | results/representative_set.csv | 7 rows | PASS | MATCH |
| K08 | Plain note for the professor on what data is kept | 2026-10-20T22:26+05:00 | docs/data_retention.md | 70 lines | PASS | MATCH |
| X02 | Distance matrix between all ESR2 structures | 2026-10-20T22:28+05:00 | results/rmsd_esr2.npy, results/rmsd_esr2_index.txt, results/esr2_reference.json | shape (86, 86) float64; 86 lines; 14 top-level keys | PASS | MATCH |
| X15 | Plot of ChEMBL potency per structural cluster | 2026-10-20T22:29+05:00 | figures/potency_by_cluster.png, results/potency_by_cluster.csv | 1125x675 px; 82 rows | PASS | MATCH |
| X16 | Statistical tests: potency against cluster | 2026-10-20T22:29+05:00 | results/potency_tests.csv | 3 rows | PASS | MATCH |
| X03 | Group the ESR2 structures the same way as ESR1 | 2026-10-20T22:31+05:00 | results/esr2_clusters.csv, results/esr2_linkage.npy | 86 rows; shape (85, 4) float64 | PASS | MATCH |
| X04 | Two-axis map (PCA) of the ESR2 structures | 2026-10-20T22:32+05:00 | figures/esr2_pca_landscape.png, results/esr2_pca_coords.csv | 1050x780 px; 86 rows | PASS | MATCH |
| X05 | ESR2 cluster sizes, with small clusters marked not interpretable | 2026-10-20T22:32+05:00 | results/esr2_basin_counts.csv | 5 rows | PASS | MATCH |
| W10 | Manuscript section: dataset and quality control | 2026-10-20T22:33+05:00 | docs/manuscript/02_dataset_qc.md | 71 lines | PASS | MATCH |
| W09 | Worked example notebook that rebuilds the main figure | 2026-10-20T22:34+05:00 | notebooks/worked_example.ipynb | 9 cells | PASS | MATCH |
| W16 | Script that collects the headline numbers for checking the manuscript | 2026-10-20T22:35+05:00 | results/numbers.json, docs/numbers.md | 1 top-level keys; 279 lines | PASS | MATCH |

Not in scope (latest state is not DONE; no evidence card):

- C08: MANUAL (2026-10-14T18:19+05:00): prepared docs/chemistry_k_choice.md (silhouettes 0.10-0.13, K=6 provisional). Human must decide whether crystallization additives (MG, CL, GOL...) stay in the chemistry set, then choose K and sign off
- X07: BLOCK (2026-10-14T18:34+05:00): needs JDK 11+ (only /usr/lib/jvm java-8-openjdk present, java -version 1.8.0_504) and prank on PATH (command -v prank fails). System JDK install needs root; not attempted. Human to install JDK 11+ and P2Rank, or approve 
- F01: PART (2026-10-20T22:30+05:00): chk pml passes, chk png fails: PyMOL not installed on this machine (no binary, no module). medoids.pml loads 7 S12 medoids with saved S03 transforms (verified vs aligned_coordinates.npz, max 2e-6 A) and colours H12 by au
- F02: PART (2026-10-20T22:30+05:00): chk pml passes, chk png fails: PyMOL not installed. displacement_structure.pml writes E07 median displacement (238 positions) into b of 5TM1\|A\|1 (largest-cluster medoid), blue_white_red 0-3 A, grey = no value. To finis

## Evidence cards

### S06: Note that distances use only C-alpha atoms
- **What was done:** Rewrote the C-alpha-only limitation so it states the current run. The paragraph says side-chain changes are not measured.
- **Code:** none (text task). Run: `none (text edited by hand)`
- **Input:** results/alignment_core_rmsd.csv (CSV, 788 rows); results/rmsd_entry.npy (NumPy, shape (788, 788) float64)
- **Output:** docs/limitations.md (md, 3 lines, sha256 b67738fe): # Limitation: C-alpha only / Every distance in this analysis is a C-alpha RMSD. The core fit (788 eligible chains from …
- **Log line:** `2026-10-20T22:07+05:00 \| S06 \| DONE \| who=opus \| chk=1/1 \| ev=docs/limitations.md,results/alignment_core_rmsd.csv,results/rmsd_entry.npy,esr_analysis/resources/config/alignment_core_v3.json \| note=updated for keep-all-data rerun: 788 chains / 351 entries, 123-residue core v3`
- **Recount:** claim "788 chains / 351 entries, 123-residue core v3"; command `(npy('results/rmsd_entry.npy').shape[0], csv('results/alignment_core_rmsd.csv').entry_id.nunique(), json.loads(txt('esr_analysis/resources/config/alignment_core_v3.json'))['residues_total'], '788 x 788' in txt('docs/limitations.md'))` → `(788, 351, 123, True)` (expected (788, 351, 123, True)): **MATCH**
- **tl.py check:** PASS  chk grep docs/limitations.md (?i)all-atom|side.chain  [docs\limitations.md]
- **Provenance:** tracker/prov/S06.json exists; 4 ev files, no drift
- **Caveats:** Text task: no producing script. The paragraph itself is unchanged except the first sentence.

### S08: Group all chains by structural similarity at four distance cuts
- **What was done:** Grouped all 788 chains by their pairwise distances at four cut-offs. Each chain also carries its quality label for reading.
- **Code:** scripts/cluster_selected.py:main L24-43; esr_analysis/cluster_structures.py:structural_groups L16-43. Run: `PYTHONPATH=. .venv/Scripts/python scripts/cluster_selected.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/rmsd_entry_index.csv (CSV, 788 rows); results/qc_tiers.csv (CSV, 526 rows)
- **Output:** results/entry_clusters.csv (CSV, 788 rows, sha256 6ce7d962): entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5,qc_label / 1A52,1A52\|A\|1,1,1,1,1,B / 1A52,1A52\|B\|1,2,1,1,1,B / 1ERE,1ERE\|A\|1,15,7,4,3,B ‖ results/linkage_entry.npy (NumPy, shape (787, 4) float64, sha256 3b6f4785): shape (787, 4), dtype float64
- **Log line:** `2026-10-20T22:07+05:00 \| S08 \| DONE \| who=opus \| chk=1/1 \| ev=results/entry_clusters.csv,results/linkage_entry.npy,scripts/cluster_selected.py \| note=788 chains/351 entries, average linkage on 788x788 core-v3 matrix; clusters 16/7/4/3 at 1.0/1.5/2.0/2.5 A; qc_label copied from K05 (A 544, B 240, C 4), not used for filtering`
- **Recount:** claim "clusters 16/7/4/3 at 1.0/1.5/2.0/2.5 A; qc_label A 544, B 240, C 4"; command `[csv('results/entry_clusters.csv')[c].nunique() for c in CUTS] + [int(fcluster(linkage(squareform(npy('results/rmsd_entry.npy'), checks=False), 'average'), 1.5, 'distance').max())] + [int((csv('results/entry_clusters.csv').qc_label == k).sum()) for k in 'ABC']` → `[16, 7, 4, 3, 7, 544, 240, 4]` (expected [16, 7, 4, 3, 7, 544, 240, 4]): **MATCH**
- **tl.py check:** PASS  chk csv results/entry_clusters.csv rows>=700 cols=entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5,qc_label nonempty=chain_key  [788 rows (need 700)]
- **Provenance:** tracker/prov/S08.json exists; 3 ev files, no drift
- **Caveats:** The 5th value re-clusters the saved matrix independently with scipy (7 at 1.5 A). qc_label is the entry-level K05 label copied per chain.

### S09: Table of cluster counts and sizes per cut
- **What was done:** Counted clusters and how much of the data sits in the three biggest at each cut-off. Also gave typical distances inside and between clusters.
- **Code:** scripts/cluster_stats.py:main L31-45; scripts/cluster_stats.py:cut_stats L19-28. Run: `PYTHONPATH=. .venv/Scripts/python scripts/cluster_stats.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/cluster_stats.csv (CSV, 4 rows, sha256 9e631cff): cut,n_chains,n_entries,n_clusters,top3_fraction,median_within_A,median_between_A / 1.0,788,351,16,0.9543147208121827,0.6964918321498813,3.417287216953521 / 1.5,788,351,7,0.9885786802030457,0.7288879654958179,3.479935270533569 / 2.0,788,351,4,0.9974619289340102,0.7340370464525812,3.4868299378698397
- **Log line:** `2026-10-20T22:14+05:00 \| S09 \| DONE \| who=opus \| chk=1/1 \| ev=results/cluster_stats.csv,scripts/cluster_stats.py \| note=788 chains/351 entries; clusters 16/7/4/3; top3 chain share 0.954/0.989/0.997/1.000; median within 0.70-0.74 A, between 3.42-3.49 A`
- **Recount:** claim "clusters 16/7/4/3; top3 chain share 0.954/0.989/0.997/1.000"; command `[csv('results/entry_clusters.csv')[c].nunique() for c in CUTS] + [round(csv('results/entry_clusters.csv')[c].value_counts().head(3).sum() / len(csv('results/entry_clusters.csv')), 3) for c in CUTS]` → `[16, 7, 4, 3, 0.954, 0.989, 0.997, 1]` (expected [16, 7, 4, 3, 0.954, 0.989, 0.997, 1]): **MATCH**
- **tl.py check:** PASS  chk csv results/cluster_stats.csv rows>=4 cols=cut,n_clusters,top3_fraction  [4 rows (need 4)]
- **Provenance:** tracker/prov/S09.json exists; 2 ev files, no drift
- **Caveats:** Spec asks for the share of entries in the 3 largest clusters; the file gives the share of chains (rows are chains). Entry-level view is K07.

### S10: Tree drawing (dendrogram) of the clustering
- **What was done:** Drew the clustering tree for all 788 chains. A dashed line marks the main 1.5 A cut-off.
- **Code:** scripts/figure_dendrogram.py:main L34-62. Run: `PYTHONPATH=. .venv/Scripts/python scripts/figure_dendrogram.py`
- **Input:** results/linkage_entry.npy (NumPy, shape (787, 4) float64); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** figures/dendrogram.png (PNG, 1272x725 px, sha256 b975cfe3): binary/structured (see shape)
- **Log line:** `2026-10-20T22:15+05:00 \| S10 \| DONE \| who=opus \| chk=1/1 \| ev=figures/dendrogram.png,scripts/figure_dendrogram.py \| note=788 leaves (351 entries), cut 1.5 A = 7 clusters; static PNG, table twin results/entry_clusters.csv`
- **Recount:** claim "788 leaves, cut 1.5 A = 7 clusters"; command `(npy('results/linkage_entry.npy').shape[0] + 1, int(fcluster(npy('results/linkage_entry.npy'), 1.5, 'distance').max()), nonwhite('figures/dendrogram.png') > 0)` → `(788, 7, True)` (expected (788, 7, True)): **MATCH**
- **tl.py check:** PASS  chk file figures/dendrogram.png  [figures\dendrogram.png]
- **Provenance:** tracker/prov/S10.json exists; 2 ev files, no drift
- **Caveats:** Static PNG; leaves unlabeled, table twin results/entry_clusters.csv.

### S12: Most central structure of each cluster
- **What was done:** Picked, for each cluster, the chain that is closest on average to the others in it. These serve as cluster representatives.
- **Code:** scripts/medoids.py:main L20-48. Run: `PYTHONPATH=. .venv/Scripts/python scripts/medoids.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/medoids.csv (CSV, 7 rows, sha256 385f02bd): cluster,entry_id,chain_key,cluster_size,cluster_entries,mean_rmsd_to_mates_A / 1,2JF9,2JF9\|A\|1,5,3,0.571782326822496 / 2,1R5K,1R5K\|B\|1,3,1,0.5166378987207454 / 3,9B25,9B25\|B\|1,3,2,0.8086270340434355
- **Log line:** `2026-10-20T22:15+05:00 \| S12 \| DONE \| who=opus \| chk=1/1 \| ev=results/medoids.csv,scripts/medoids.py \| note=7 medoids at 1.5 A over 788 chains; sizes 5/3/3/329/1/2/445 chains; 9B25 chains A and B fall in different clusters (3 and 6)`
- **Recount:** claim "sizes 5/3/3/329/1/2/445 chains; largest-cluster medoid 5TM1|A|1"; command `csv('results/medoids.csv').cluster_size.tolist() + [medoid(7)]` → `[5, 3, 3, 329, 1, 2, 445, '5TM1\|A\|1']` (expected [5, 3, 3, 329, 1, 2, 445, '5TM1\|A\|1']): **MATCH**
- **tl.py check:** PASS  chk csv results/medoids.csv rows>=1 cols=cluster,entry_id  [7 rows (need 1)]
- **Provenance:** tracker/prov/S12.json exists; 2 ev files, no drift
- **Caveats:** Spec says 'the entry'; the file picks a chain (an entry can have chains in two clusters: 9B25 A and B are in clusters 6 and 3).

### S13: Does keeping only good-quality structures change the clusters?
- **What was done:** Re-grouped only the better-quality chains and compared with the full grouping. Nothing was removed from the files.
- **Code:** scripts/qc_clustering_delta.py:main L34-58. Run: `PYTHONPATH=. .venv/Scripts/python scripts/qc_clustering_delta.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/qc_clustering_delta.csv (CSV, 12 rows, sha256 7dd247a8): cut,subset,n_chains,n_entries,n_clusters,top3_fraction,clusters_all,clusters_without_membe… / 1.0,all,788,351,16,0.9543147208121827,16,0,1.0 / 1.5,all,788,351,7,0.9885786802030457,7,0,1.0 / 2.0,all,788,351,4,0.9974619289340102,4,0,1.0
- **Log line:** `2026-10-20T22:21+05:00 \| S13 \| DONE \| who=opus \| chk=1/1 \| ev=results/qc_clustering_delta.csv,scripts/qc_clustering_delta.py \| note=slices of S04 788 matrix, S08 method; tier A+B (784 chains) identical to all at every cut (ARI 1.0); tier A only (544 chains/247 entries): clusters 14/4/3/2 vs 16/7/4/3, ARI 0.91/0.98/0.99/0.99, 2/1/0/0 all-chain clusters with no tier-A member; conclusions hold`
- **Recount:** claim "tier A only clusters 14/4/3/2; tier A+B ARI 1.0 at every cut"; command `[nclust(csv('results/entry_clusters.csv').qc_label.eq('A').to_numpy(), t) for t in (1.0, 1.5, 2.0, 2.5)] + csv('results/qc_clustering_delta.csv').query("subset == 'tier_A+B'").ari_vs_all.round(3).tolist()` → `[14, 4, 3, 2, 1, 1, 1, 1]` (expected [14, 4, 3, 2, 1, 1, 1, 1]): **MATCH**
- **tl.py check:** PASS  chk csv results/qc_clustering_delta.csv rows>=8 cols=cut,subset,n_chains,n_clusters,top3_fraction,clusters_without_members,ari_vs_all  [12 rows (need 8)]
- **Provenance:** tracker/prov/S13.json exists; 2 ev files, no drift
- **Caveats:** Tier-A cluster counts re-derived by re-clustering the matrix slice independently.

### S14: Plot of cluster count against distance cut
- **What was done:** Plotted how the number of clusters falls as the cut-off grows. The four reported cut-offs are marked.
- **Code:** scripts/figure_sensitivity.py:main L44-96; scripts/figure_sensitivity.py:stats_at L38-41. Run: `PYTHONPATH=. .venv/Scripts/python scripts/figure_sensitivity.py`
- **Input:** results/linkage_entry.npy (NumPy, shape (787, 4) float64); results/cluster_stats.csv (CSV, 4 rows)
- **Output:** figures/structural_sensitivity.png (PNG, 1480x654 px, sha256 3717501e): binary/structured (see shape)
- **Log line:** `2026-10-20T22:23+05:00 \| S14 \| DONE \| who=opus \| chk=1/1 \| ev=figures/structural_sensitivity.png,scripts/figure_sensitivity.py \| note=788 chains, sweep 1.0-2.5 A every 0.05 A from S08 linkage; reported cuts match cluster_stats.csv (16/7/4/3; top-3 share 95.4-100%)`
- **Recount:** claim "reported cuts 16/7/4/3"; command `[int(fcluster(npy('results/linkage_entry.npy'), t, 'distance').max()) for t in (1.0, 1.5, 2.0, 2.5)] + [nonwhite('figures/structural_sensitivity.png') > 0]` → `[16, 7, 4, 3, True]` (expected [16, 7, 4, 3, True]): **MATCH**
- **tl.py check:** PASS  chk file figures/structural_sensitivity.png  [figures\structural_sensitivity.png]
- **Provenance:** tracker/prov/S14.json exists; 2 ev files, no drift
- **Caveats:** Static PNG; numbers in the figure checked against cluster_stats.csv by the producing script.

### S07: Compare C-alpha distances with all-atom distances on 50 entries
- **What was done:** Measured distances with all heavy atoms on 50 entries and compared them with the C-alpha distances. This tests whether C-alpha only changes the ordering.
- **Code:** scripts/allatom_rmsd_subset.py:main L111-175; scripts/allatom_rmsd_subset.py:pair_rmsd L99-108. Run: `PYTHONPATH=. .venv/Scripts/python scripts/allatom_rmsd_subset.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/alignment/alignment_transforms.json (JSON, 788 top-level keys); results/lbd_extract/residue_mapping.csv (CSV, 222348 rows)
- **Output:** results/rmsd_allatom_subset.csv (CSV, 1225 rows, sha256 41f7dcf6): chain_a,chain_b,entry_a,entry_b,ca_rmsd,allatom_rmsd,ca_rmsd_same_residues,n_heavy_atoms,n… / 1A52\|B\|1,1R5K\|A\|1,1A52,1R5K,6.357535094375079,6.055490516080119,6.371076037532446,1847,234… / 1A52\|B\|1,1XP1\|A\|1,1A52,1XP1,7.190611196998688,6.769347519273019,7.2056367715578515,1865,23… / 1A52\|B\|1,1XP9\|A\|1,1A52,1XP9,7.153808558974437,6.739734728573631,7.168759407740879,1865,237… ‖ results/rmsd_allatom_subset_summary.json (JSON, 17 top-level keys, sha256 00190436): { /   "task": "S07", /   "n_entries": 50,
- **Log line:** `2026-10-20T22:14+05:00 \| S07 \| DONE \| who=opus \| chk=1/1 \| ev=results/rmsd_allatom_subset.csv,results/rmsd_allatom_subset_summary.json,scripts/allatom_rmsd_subset.py \| note=50 entries (historical chain, stratified by S08 cut_1.5, seed 20261020), 1225 pairs, heavy atoms after saved S03 core fit, no refit; Spearman CA vs all-atom 0.985, per-chain median 0.952 (min 0.823 2JF9\|A\|1); ARI same K 0.92-1.0; 50 mmCIFs from RCSB match sha256sums`
- **Recount:** claim "1225 pairs over 50 entries; Spearman CA vs all-atom 0.985"; command `(lambda d: (len(d), round(spearmanr(d.ca_rmsd, d.allatom_rmsd)[0], 3), pd.concat([d.entry_a, d.entry_b]).nunique()))(csv('results/rmsd_allatom_subset.csv'))` → `(1225, 0.985, 50)` (expected (1225, 0.985, 50)): **MATCH**
- **tl.py check:** PASS  chk csv results/**/rmsd_allatom*.csv rows>=50  [1225 rows (need 50)]
- **Provenance:** tracker/prov/S07.json exists; 3 ev files, no drift
- **Caveats:** Uses the saved S03 superposition for both measures (no all-atom refit). Needs mmCIFs in data/raw (git-ignored). Subset = historical chain per entry, stratified by cluster.

### E05: Superimpose each matched ESR2 chain on its ESR1 partner
- **What was done:** Laid each ESR2 chain over its ESR1 partner that holds the same ligand. The fit uses the stable core helices.
- **Code:** scripts/pair_superposition.py:main L40-86; scripts/pair_superposition.py:shared_rmsd L33-37. Run: `PYTHONPATH=. .venv/Scripts/python scripts/pair_superposition.py`
- **Input:** results/matched_pairs.csv (CSV, 135 rows); results/lbd_chains.csv (CSV, 987 rows); results/lbd_extract/residue_mapping.csv (CSV, 222348 rows)
- **Output:** results/pair_superposition.csv (CSV, 671 rows, sha256 c3ffa2fe): esr1_entry,esr2_entry,ligand,esr1_chain,esr2_chain,primary_pair,n_common_core_ca,core_rmsd… / 1A52,3OLL,EST,1A52\|A\|1,3OLL\|A\|1,False,123,0.6458521973880651,231,5.666900167438007,7,28.15… / 1A52,3OLL,EST,1A52\|A\|1,3OLL\|B\|1,False,123,0.6439060650708466,230,5.671379421522448,7,28.14… / 1A52,3OLL,EST,1A52\|B\|1,3OLL\|A\|1,True,123,0.6349131557583904,231,5.561594306519177,7,27.903…
- **Log line:** `2026-10-20T22:16+05:00 \| E05 \| DONE \| who=opus \| chk=1/1 \| ev=results/pair_superposition.csv,scripts/pair_superposition.py \| note=135 entry pairs -> 671 chain pairs (all eligible chains, 135 flagged primary_pair), all fitted on core v3 (115-123 shared CA); core RMSD median 0.524 A (max 1.876, 1QKT\|A\|1 vs 5TOA); LBD median 1.434 A; H12 median 1.957 A`
- **Recount:** claim "671 chain pairs, 135 primary; core median 0.524 (max 1.876), LBD median 1.434, H12 median 1.957 A"; command `(lambda d: (len(d), int(d.primary_pair.sum()), d.core_rmsd_A.median(), d.core_rmsd_A.max(), d.lbd_rmsd_A.median(), d.h12_rmsd_A.median()))(csv('results/pair_superposition.csv'))` → `(671, 135, 0.523739, 1.87637, 1.43361, 1.95702)` (expected (671, 135, 0.524, 1.876, 1.434, 1.957)): **MATCH**
- **tl.py check:** PASS  chk csv results/pair_superposition.csv rows>=5  [671 rows (need 5)]
- **Provenance:** tracker/prov/E05.json exists; 2 ev files, no drift
- **Caveats:** All chain combinations kept; primary_pair marks one pair per entry pair.

### R01: How evenly the chains spread over clusters (entropy)
- **What was done:** Measured how evenly chains spread across clusters (entropy). Low values mean many near-identical structures.
- **Code:** scripts/cluster_entropy.py:main L25-43; scripts/cluster_entropy.py:entropy_bits L19-22. Run: `PYTHONPATH=. .venv/Scripts/python scripts/cluster_entropy.py`
- **Input:** results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/entropy.csv (CSV, 4 rows, sha256 a7dd4a77): cut,n_chains,n_entries,n_clusters,largest_cluster_chains,entropy_bits,entropy_bits_entry_w… / 1.0,788,351,16,442,1.595844356895745,1.542466102520555,4.0,0.39896108922393625 / 1.5,788,351,7,445,1.1332722559505624,1.0765573203063588,2.807354922057604,0.40367972251971… / 2.0,788,351,4,445,1.0581080713180215,1.0150736821025923,2.0,0.5290540356590108
- **Log line:** `2026-10-20T22:16+05:00 \| R01 \| DONE \| who=opus \| chk=1/1 \| ev=results/entropy.csv,scripts/cluster_entropy.py \| note=788 chains; entropy 1.60/1.13/1.06/1.03 bits at 1.0/1.5/2.0/2.5 A (entry-weighted 1.54/1.08/1.02/1.00); normalised 0.40-0.65: highly redundant`
- **Recount:** claim "entropy 1.60/1.13/1.06/1.03 bits"; command `(lambda d: [float(-(p * np.log2(p)).sum()) for p in (d[c].value_counts(normalize=True) for c in CUTS)])(csv('results/entry_clusters.csv'))` → `[1.59584, 1.13327, 1.05811, 1.03465]` (expected [1.6, 1.13, 1.06, 1.03]): **MATCH**
- **tl.py check:** PASS  chk csv results/entropy.csv rows>=4 cols=cut,entropy_bits  [4 rows (need 4)]
- **Provenance:** tracker/prov/R01.json exists; 2 ev files, no drift
- **Caveats:** none found

### R02: How many clusters appear as structures are added at random
- **What was done:** Added structures in 100 random orders and counted clusters seen. This gives a mean curve with its spread.
- **Code:** scripts/diversity_curve.py:main L26-42; scripts/provisional_common.py:clusters_seen L63-75. Run: `PYTHONPATH=. .venv/Scripts/python scripts/diversity_curve.py`
- **Input:** results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/diversity_curve.csv (CSV, 3152 rows, sha256 59cfc5d2): cut,n,mean_clusters,sd_clusters,n_clusters_total,n_90,seed,n_orderings / 1.0,1,1.0,0.0,16,520,42,100 / 1.0,2,1.54,0.5009082659620332,16,520,42,100 / 1.0,3,1.91,0.6046119049072346,16,520,42,100
- **Log line:** `2026-10-20T22:17+05:00 \| R02 \| DONE \| who=opus \| chk=1/1 \| ev=results/diversity_curve.csv,scripts/diversity_curve.py \| note=788 chains x 4 cuts, 100 orderings seed 42; n_90 = 520/488/315/147 chains at 1.0/1.5/2.0/2.5 A (late saturation driven by singleton/small clusters)`
- **Recount:** claim "3152 rows; n_90 = 520/488/315/147"; command `(lambda d: (len(d), [int(g.n[g.mean_clusters >= 0.9 * g.n_clusters_total].min()) for _, g in d.groupby('cut')]))(csv('results/diversity_curve.csv'))` → `(3152, [520, 488, 315, 147])` (expected (3152, [520, 488, 315, 147])): **MATCH**
- **tl.py check:** PASS  chk csv results/diversity_curve.csv rows>=300 cols=n,mean_clusters,sd_clusters  [3152 rows (need 300)]
- **Provenance:** tracker/prov/R02.json exists; 2 ev files, no drift
- **Caveats:** none found

### E06: Per-residue shift between ESR1 and ESR2 for the first 5 pairs
- **What was done:** For the first five ESR1/ESR2 pairs, measured how far each residue moves after the fit. One row per pair and residue.
- **Code:** scripts/receptor_displacements.py:main L49-90; scripts/receptor_displacements.py:displacements L33-46. Run: `PYTHONPATH=. .venv/Scripts/python scripts/receptor_displacements.py --first 5`
- **Input:** results/pair_superposition.csv (CSV, 671 rows); results/lbd_extract/residue_mapping.csv (CSV, 222348 rows); results/numbering_map_esr1_esr2.csv (CSV, 520 rows)
- **Output:** results/first_five_displacements.csv (CSV, 1125 rows, sha256 6a6abaca): pair,esr1_entry,esr2_entry,ligand,esr1_chain,esr2_chain,esr1_pos,esr2_pos,esr1_aa,esr2_aa,… / 1,1A52,3OLL,EST,1A52\|B\|1,3OLL\|A\|1,309,262,S,A,False,False,False,4.960943387628074 / 1,1A52,3OLL,EST,1A52\|B\|1,3OLL\|A\|1,310,263,L,L,True,False,False,2.3917722009654714 / 1,1A52,3OLL,EST,1A52\|B\|1,3OLL\|A\|1,311,264,T,S,False,False,False,1.8571620399721274
- **Log line:** `2026-10-20T22:18+05:00 \| E06 \| DONE \| who=opus \| chk=1/1 \| ev=results/first_five_displacements.csv,scripts/receptor_displacements.py \| note=first 5 primary pairs (all 1A52\|B\|1 vs EST-bound ESR2), 1125 position rows; median 0.53-0.73 A; >30 A only at H12 538-544 where 1A52 H12 is displaced (S16 h12_rmsd 28.4 A); refit core RMSD equals E05`
- **Recount:** claim "1125 rows, 5 pairs; per-pair median 0.53-0.73 A; >30 A only at 538-544"; command `(lambda d: (len(d), d.pair.nunique(), d.groupby('pair').displacement_A.median().min(), d.groupby('pair').displacement_A.median().max(), int(d.esr1_pos[d.displacement_A > 30].min()), int(d.esr1_pos[d.displacement_A > 30].max())))(csv('results/first_five_displacements.csv'))` → `(1125, 5, 0.529464, 0.734414, 540, 544)` (expected (1125, 5, 0.53, 0.73, 538, 544)): **MISMATCH**
- **tl.py check:** PASS  chk csv results/first_five_displacements.csv rows>=200  [1125 rows (need 200)]
- **Provenance:** tracker/prov/E06.json exists; 2 ev files, no drift
- **Caveats:** All five first pairs use the same ESR1 chain (1A52|B|1, estradiol), so they are not a varied sample. Related NOTE: `2026-10-20T23:18+05:00 \| - \| NOTE \| who=opus \| note=correction to E06: displacements above 30 A occur at ESR1 positions 540-544 (recount of results/first_five_displacements.csv), not 538-544 as the DONE note says` Positions with values > 30 A are exactly 540, 541, 543, 544.

### E07: Per-residue shift for all matched pairs (position x pair table)
- **What was done:** Did the same for all 135 matched entry pairs. The result is a residue-by-pair table.
- **Code:** scripts/receptor_displacements.py:main L49-90. Run: `PYTHONPATH=. .venv/Scripts/python scripts/receptor_displacements.py --all`
- **Input:** results/pair_superposition.csv (CSV, 671 rows); results/lbd_extract/residue_mapping.csv (CSV, 222348 rows)
- **Output:** results/displacement_matrix.csv (CSV, 246 rows, sha256 2b939e9f): esr1_pos,1A52\|B\|1~3OLL\|A\|1,1A52\|B\|1~3OLS\|A\|1,1A52\|B\|1~4J24\|A\|1,1A52\|B\|1~4J26\|A\|1,1A52\|B\|1~… / 305,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,… / 306,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,… / 307,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,5.521700695938096,,,,,,,,,,,5.322009504468052,,,,… ‖ results/displacement_pairs.csv (CSV, 135 rows, sha256 12a44cb0): column,esr1_entry,esr2_entry,ligand,esr1_chain,esr2_chain,core_rmsd_A / 1A52\|B\|1~3OLL\|A\|1,1A52,3OLL,EST,1A52\|B\|1,3OLL\|A\|1,0.6349131557583904 / 1A52\|B\|1~3OLS\|A\|1,1A52,3OLS,EST,1A52\|B\|1,3OLS\|A\|1,0.4957956597922236 / 1A52\|B\|1~4J24\|A\|1,1A52,4J24,EST,1A52\|B\|1,4J24\|A\|1,0.4868279154174855
- **Log line:** `2026-10-20T22:18+05:00 \| E07 \| DONE \| who=opus \| chk=1/1 \| ev=results/displacement_matrix.csv,results/displacement_pairs.csv,scripts/receptor_displacements.py \| note=246 ESR1 LBD positions x 135 primary pairs; blank = unobserved in either chain; column key in displacement_pairs.csv`
- **Recount:** claim "246 positions x 135 pairs"; command `(lambda m: (m.shape[0], m.shape[1] - 1))(csv('results/displacement_matrix.csv'))` → `(246, 135)` (expected (246, 135)): **MATCH**
- **tl.py check:** PASS  chk csv results/displacement_matrix.csv rows>=200  [246 rows (need 200)]
- **Provenance:** tracker/prov/E07.json exists; 3 ev files, no drift
- **Caveats:** Columns are the 135 primary pairs only, not all 671 chain pairs of E05. Related NOTE: `2026-10-20T22:31+05:00 \| - \| NOTE \| who=opus \| note=correction to F02 PART note: displacement_structure.pml carries 219 positions with a value (positions observed in 5TM1\|A\|1 with a non-blank E07 median), not 238`

### R04: Number of structures needed to see 90% of clusters
- **What was done:** Read off how many structures are needed to see 90% of clusters. One number per cut-off.
- **Code:** scripts/saturation.py:main L15-30. Run: `PYTHONPATH=. .venv/Scripts/python scripts/saturation.py`
- **Input:** results/diversity_curve.csv (CSV, 3152 rows)
- **Output:** results/saturation.csv (CSV, 4 rows, sha256 3086465d): cut,n_90,n_clusters_total,threshold_90pct,mean_at_n90,n_chains,fraction_of_chains,seed,n_o… / 1.0,520,16,14.4,14.41,788,0.6598984771573604,42,100 / 1.5,488,7,6.3,6.3,788,0.6192893401015228,42,100 / 2.0,315,4,3.6,3.61,788,0.399746192893401,42,100
- **Log line:** `2026-10-20T22:18+05:00 \| R04 \| DONE \| who=opus \| chk=1/1 \| ev=results/saturation.csv,scripts/saturation.py \| note=n_90 = 520/488/315/147 of 788 chains (66/62/40/19%) at 1.0/1.5/2.0/2.5 A; agrees with R02`
- **Recount:** claim "n_90 = 520/488/315/147 (66/62/40/19% of chains)"; command `(lambda d: (d.n_90.tolist(), d.fraction_of_chains.tolist()))(csv('results/saturation.csv'))` → `([520, 488, 315, 147], [0.659898, 0.619289, 0.399746, 0.186548])` (expected ([520, 488, 315, 147], [0.66, 0.62, 0.4, 0.19])): **MATCH**
- **tl.py check:** PASS  chk csv results/saturation.csv rows>=4 cols=cut,n_90  [4 rows (need 4)]
- **Provenance:** tracker/prov/R04.json exists; 2 ev files, no drift
- **Caveats:** none found

### R03: Figure of the cumulative diversity curve with a 90% line
- **What was done:** Drew the diversity curves with a spread band and a 90% line. Saturation points are marked.
- **Code:** scripts/figure_diversity_curve.py:main L23-56. Run: `PYTHONPATH=. .venv/Scripts/python scripts/figure_diversity_curve.py`
- **Input:** results/diversity_curve.csv (CSV, 3152 rows); results/saturation.csv (CSV, 4 rows)
- **Output:** figures/diversity_curve.png (PNG, 1200x720 px, sha256 fb06e29f): binary/structured (see shape)
- **Log line:** `2026-10-20T22:18+05:00 \| R03 \| DONE \| who=opus \| chk=1/1 \| ev=figures/diversity_curve.png,scripts/figure_diversity_curve.py \| note=4 cuts, mean +/- 1 SD band, dashed 90% line, n_90 markers; table twins results/diversity_curve.csv, results/saturation.csv`
- **Recount:** claim "4 cuts drawn, non-blank figure"; command `(len(csv('results/saturation.csv')), nonwhite('figures/diversity_curve.png') > 0)` → `(4, True)` (expected (4, True)): **MATCH**
- **tl.py check:** PASS  chk file figures/diversity_curve.png  [figures\diversity_curve.png]
- **Provenance:** tracker/prov/R03.json exists; 2 ev files, no drift
- **Caveats:** Static PNG.

### Q14: Are rare structural groups worse in resolution or R-free?
- **What was done:** Compared resolution and R-free of entries in rare clusters with the rest. Gave a test line, a table and a plot.
- **Code:** scripts/qc_outliers.py:main L34-64. Run: `PYTHONPATH=. .venv/Scripts/python scripts/qc_outliers.py`
- **Input:** results/entry_clusters.csv (CSV, 788 rows); results/qc_tiers.csv (CSV, 526 rows)
- **Output:** results/qc_outliers.csv (CSV, 788 rows, sha256 2872721c): entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5,resolution_A,r_free,mean_b_A2,qc_label,… / 1A52,1A52\|A\|1,1,1,1,1,2.8,0.274,25.43,B,2,True,5,True,5,True,5,True / 1A52,1A52\|B\|1,2,1,1,1,2.8,0.274,25.43,B,3,True,5,True,5,True,5,True / 1ERE,1ERE\|A\|1,15,7,4,3,3.1,0.251,57.48,B,442,False,445,False,445,False,447,False ‖ docs/qc_outliers.md (md, 39 lines, sha256 f31f2425): # QC question 1: do conformational singletons have worse resolution? / Generated by `scripts/qc_outliers.py` from `results/entry_clusters.csv` (S08, 788 chains) … / A cluster is rare when it holds fewer than 10 chains. An entry is rare when any of its cha…
- **Log line:** `2026-10-20T22:20+05:00 \| Q14 \| DONE \| who=opus \| chk=2/2 \| ev=results/qc_outliers.csv,docs/qc_outliers.md,results/qc_outliers_stats.csv,figures/qc_outliers.png,scripts/qc_outliers.py \| note=per-entry Mann-Whitney, rare = cluster < 10 chains; at 1.5 A 7 rare vs 344 common entries: resolution median 2.40 vs 2.10 A (p 0.21), R-free 0.236 vs 0.238 (p 0.75); no cut reaches p<0.05 (lowest 0.061 resolution at 2.5 A, n=3)`
- **Recount:** claim "1.5 A: 7 rare vs 344 common entries; resolution 2.40 vs 2.10 A (p 0.21); R-free p 0.75"; command `(lambda e: (int(e.rare.sum()), int((~e.rare).sum()), e.resolution_A[e.rare].median(), e.resolution_A[~e.rare].median(), mannwhitneyu(e.resolution_A[e.rare], e.resolution_A[~e.rare]).pvalue, mannwhitneyu(e.r_free[e.rare], e.r_free[~e.rare]).pvalue))(csv('results/qc_outliers.csv').groupby('entry_id').agg(rare=('rare_1.5', 'any'), resolution_A=('resolution_A', 'first'), r_free=('r_free', 'first')))` → `(7, 344, 2.4, 2.1, 0.206738, 0.747674)` (expected (7, 344, 2.4, 2.1, 0.21, 0.75)): **MATCH**
- **tl.py check:** PASS  chk csv results/qc_outliers.csv rows>=1  [788 rows (need 1)]; PASS  chk file docs/qc_outliers.md  [docs\qc_outliers.md]
- **Provenance:** tracker/prov/Q14.json exists; 5 ev files, no drift
- **Caveats:** 'Rare' = cluster with fewer than 10 chains: a threshold chosen in the script, not given by the spec.

### Q15: Do high B-factor structures gather in particular clusters?
- **What was done:** Cross-tabulated high, mid and low B-factor chains against clusters. Tested whether high-B chains gather together.
- **Code:** scripts/qc_bfactor_clusters.py:main L28-75. Run: `PYTHONPATH=. .venv/Scripts/python scripts/qc_bfactor_clusters.py`
- **Input:** results/entry_clusters.csv (CSV, 788 rows); results/helix_bfactor_by_chain.csv (CSV, 788 rows); results/qc_tiers.csv (CSV, 526 rows)
- **Output:** results/qc_bfactor_clusters.csv (CSV, 7 rows, sha256 80a0e7ad): cluster,n_low,n_mid,n_high,n_chains,n_entries,median_mean_b_chain,high_share,high_enrichme… / 1,1,0,4,5,3,58.396276150627614,0.8,2.396958174904943,False,,,,1.5 / 2,0,0,3,3,1,46.95953389830509,1.0,2.996197718631179,False,,,,1.5 / 3,3,0,0,3,2,30.9508085106383,0.0,0.0,False,,,,1.5 ‖ results/qc_bfactor_clusters_tests.json (JSON, 7 top-level keys, sha256 dfb76439): { /   "task": "Q15", /   "cut": "1.5",
- **Log line:** `2026-10-20T22:20+05:00 \| Q15 \| DONE \| who=opus \| chk=1/1 \| ev=results/qc_bfactor_clusters.csv,results/qc_bfactor_clusters_tests.json,scripts/qc_bfactor_clusters.py \| note=high-B concentrates, not scatters: at 1.5 A high-B tertile share 42% in cluster 7 vs 21% in cluster 4 (chain chi2 p 5.6e-12); entry level tier_b B 66/219 vs 16/124 (chi2 p 7.4e-4); small clusters 1-2 mostly high-B (7/8 chains), listed untested`
- **Recount:** claim "high-B share 42% (cluster 7) vs 21% (cluster 4), chain chi2 p 5.6e-12; clusters 1-2: 7 of 8 chains high-B"; command `(lambda d: ((d[d['cut_1.5'] == 7].b_group == 'high').mean(), (d[d['cut_1.5'] == 4].b_group == 'high').mean(), f"{chi2_contingency(pd.crosstab(d[d['cut_1.5'].isin([4, 7])]['cut_1.5'], d.b_group))[1]:.1e}", int((d[d['cut_1.5'].isin([1, 2])].b_group == 'high').sum()), int(d['cut_1.5'].isin([1, 2]).sum())))(csv('results/entry_clusters.csv').merge(csv('results/helix_bfactor_by_chain.csv')[['chain_key', 'b_group']], on='chain_key'))` → `(0.420225, 0.209726, '5.6e-12', 7, 8)` (expected (0.42, 0.21, '5.6e-12', 7, 8)): **MATCH**
- **tl.py check:** PASS  chk csv results/qc_bfactor_clusters.csv rows>=1  [7 rows (need 1)]
- **Provenance:** tracker/prov/Q15.json exists; 3 ev files, no drift
- **Caveats:** Chain-level counts repeat entries; the entry-level test (p 7.4e-4) is the safer reading. Only clusters with >= 20 chains are tested.

### K06: One status row for every ESR1 chain, none dropped
- **What was done:** Built one row for each of the 901 ESR1 chains with all their labels. Chains outside the matrix keep their row and a reason.
- **Code:** scripts/chain_status.py:main L23-63. Run: `PYTHONPATH=. .venv/Scripts/python scripts/chain_status.py`
- **Input:** results/chain_loss_reasons.csv (CSV, 901 rows); results/qc_tiers.csv (CSV, 526 rows); results/helix_bfactor_by_chain.csv (CSV, 788 rows); results/h12_metric.csv (CSV, 788 rows)
- **Output:** results/chain_status.csv (CSV, 901 rows, sha256 3c8aceac): chain_key,entry_id,chain_id,auth_asym_id,mapped_lbd_ca,mapped_core_ca,completeness_status,… / 11SD\|C\|1,11SD,C,X,0,0,incomplete,"no declared LBD positions in this construct: no LBD over… / 11SD\|D\|1,11SD,D,Y,0,0,incomplete,"no declared LBD positions in this construct: no LBD over… / 11SE\|C\|1,11SE,C,X,0,0,incomplete,"no declared LBD positions in this construct: no LBD over…
- **Log line:** `2026-10-20T22:21+05:00 \| K06 \| DONE \| who=opus \| chk=1/1 \| ev=results/chain_status.csv,scripts/chain_status.py \| note=901 ESR1 chains / 444 entries, none dropped; in_matrix 788, out 113 all with completeness-screen reason (no LBD overlap 109, below LBD/core minimum 4); qc_label A 628 / B 261 / C 4 / unknown 8; old pipeline: 351 kept, 437 one-chain-per-entry, 113 completeness`
- **Recount:** claim "901 chains / 444 entries; in matrix 788, out 113; qc A 628 / B 261 / C 4 / unknown 8"; command `(lambda d: (len(d), d.entry_id.nunique(), int(d.in_matrix.sum()), int((~d.in_matrix).sum()), [int((d.qc_label == k).sum()) for k in ['A', 'B', 'C', 'unknown']]))(csv('results/chain_status.csv'))` → `(901, 444, 788, 113, [628, 261, 4, 8])` (expected (901, 444, 788, 113, [628, 261, 4, 8])): **MATCH**
- **tl.py check:** PASS  chk csv results/chain_status.csv rows>=788 cols=chain_key,entry_id,completeness_status,qc_label,h12_b_group,in_matrix,warning nonempty=chain_key,entry_id,in_matrix  [901 rows (need 788)]
- **Provenance:** tracker/prov/K06.json exists; 2 ev files, no drift
- **Caveats:** none found

### K09: Does the choice of alignment core change the clusters?
- **What was done:** Rebuilt the distance matrix using the old alignment core as a new file. Compared its clusters with the main result.
- **Code:** scripts/core_sensitivity.py:main L35-67. Run: `PYTHONPATH=. .venv/Scripts/python scripts/core_sensitivity.py`
- **Input:** results/alignment/aligned_coordinates_old_core.npz (NumPy npz, chain_keys (788,), positions (246,), coords (788, 246, 3)); results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** results/core_sensitivity.csv (CSV, 2 rows, sha256 7dc3c5a1): core,core_residues,n_chains,spearman_pairwise,median_abs_diff_A,max_abs_diff_A,n_clusters_… / alignment_core_v3 (primary),123,788,1.0,0.0,0.0,16,1.0,7,1.0,4,1.0,3,1.0,1.0 / old_core,84,788,0.987749902642129,0.04410863340130783,2.217611800209484,20,0.7748116661675… ‖ results/rmsd_old_core.npy (NumPy, shape (788, 788) float64, sha256 608405e4): shape (788, 788), dtype float64
- **Log line:** `2026-10-20T22:23+05:00 \| K09 \| DONE \| who=opus \| chk=2/2 \| ev=results/core_sensitivity.csv,results/rmsd_old_core.npy,scripts/core_sensitivity.py \| note=old-core 788x788 matrix (new file; equals archived strict matrix on the 351 shared chains, max diff 0.0); vs primary core v3: Spearman pairwise 0.988, median \|diff\| 0.044 A (max 2.22); clusters at 1.0/1.5/2.0/2.5 A 20/7/6/4 vs 16/7/4/3, ARI 0.77/0.96/0.97/0.98. Core choice barely changes the primary 1.5 A clusters; finer cut 1.0 A is core-sensitive`
- **Recount:** claim "old-core matrix 788x788; Spearman 0.988; median |diff| 0.044 A (max 2.22); old-core clusters 20/7/6/4"; command `(lambda a, b: (a.shape, spearmanr(up(a), up(b))[0], float(np.median(np.abs(up(a) - up(b)))), float(np.abs(up(a) - up(b)).max()), [int(fcluster(linkage(squareform(a, checks=False), 'average'), t, 'distance').max()) for t in (1.0, 1.5, 2.0, 2.5)]))(npy('results/rmsd_old_core.npy'), npy('results/rmsd_entry.npy'))` → `((788, 788), 0.98775, 0.0441086, 2.21761, [20, 7, 6, 4])` (expected ((788, 788), 0.988, 0.044, 2.22, [20, 7, 6, 4])): **MATCH**
- **tl.py check:** PASS  chk csv results/core_sensitivity.csv rows>=1 cols=core,n_chains,n_clusters_1_5,ari_vs_primary,spearman_pairwise  [2 rows (need 1)]; PASS  chk npy results/rmsd_old_core.npy square>=700  [shape (788, 788) (need square >= 700)]
- **Provenance:** tracker/prov/K09.json exists; 3 ev files, no drift
- **Caveats:** Spec allows this new matrix explicitly (the primary matrix is untouched).

### E08: Heatmap of the ESR1/ESR2 per-residue shifts
- **What was done:** Drew the residue-by-pair shift table as a heatmap. Grey marks residues not seen in a pair.
- **Code:** scripts/figure_displacement_heatmap.py:main L56-75; scripts/figure_displacement_heatmap.py:draw L31-53. Run: `PYTHONPATH=. .venv/Scripts/python scripts/figure_displacement_heatmap.py`
- **Input:** results/displacement_matrix.csv (CSV, 246 rows); results/displacement_pairs.csv (CSV, 135 rows)
- **Output:** figures/displacement_heatmap.png (PNG, 1472x899 px, sha256 1e01a612): binary/structured (see shape)
- **Log line:** `2026-10-20T22:24+05:00 \| E08 \| DONE \| who=opus \| chk=1/1 \| ev=figures/displacement_heatmap.png,scripts/figure_displacement_heatmap.py \| note=246 positions x 135 pairs, sequential blue clipped at 5 A, grey = unobserved; core v3 and H12 strip; table twin results/displacement_matrix.csv`
- **Recount:** claim "246 positions drawn, non-blank figure"; command `(csv('results/displacement_matrix.csv').shape[0], nonwhite('figures/displacement_heatmap.png') > 0)` → `(246, True)` (expected (246, True)): **MATCH**
- **tl.py check:** PASS  chk file figures/displacement_heatmap.png  [figures\displacement_heatmap.png]
- **Provenance:** tracker/prov/E08.json exists; 2 ev files, no drift
- **Caveats:** Static PNG; colour clipped at 5 A.

### E10: Positions that shift by more than 1.5 A in most pairs
- **What was done:** Listed residues that shift by more than 1.5 A in most pairs. Every residue is kept with its counts.
- **Code:** scripts/consistent_positions.py:main L23-48. Run: `PYTHONPATH=. .venv/Scripts/python scripts/consistent_positions.py`
- **Input:** results/displacement_matrix.csv (CSV, 246 rows); results/numbering_map_esr1_esr2.csv (CSV, 520 rows)
- **Output:** results/consistent_positions.csv (CSV, 246 rows, sha256 c40cb7da): position,n_pairs_observed,n_pairs_above,fraction_pairs,median_displacement_A,threshold_A,e… / 305,0,0,,,1.5,258.0,S,L,False,False,False,False / 306,0,0,,,1.5,259.0,L,L,True,False,False,False / 307,3,3,1.0,5.322009504468052,1.5,260.0,A,L,False,False,False,False
- **Log line:** `2026-10-20T22:24+05:00 \| E10 \| DONE \| who=opus \| chk=1/1 \| ev=results/consistent_positions.csv,scripts/consistent_positions.py \| note=all 246 positions written; 48 consistent (>1.5 A in >50% of >=20 observed pairs): none in core v3 (fit region), 7 in H12 (538-547), 29 with a different ESR1/ESR2 residue; main blocks N-term 308-311, 326-334 (next to the ESR1 331-332 gap in the ESR2 map), 365-369, loop 458-491, 535-547`
- **Recount:** claim "48 consistent positions; 0 in core, 7 in H12, 29 with a different residue"; command `(lambda m, nm: (lambda hit: (int(hit.sum()), int((hit & m.index.isin(core_positions())).sum()), int((hit & m.index.isin(range(538, 548))).sum()), int((hit & np.array([nm.get(p) == False for p in m.index])).sum())))(((m > 1.5).sum(axis=1) / m.notna().sum(axis=1) > 0.5) & (m.notna().sum(axis=1) >= 20)))(csv('results/displacement_matrix.csv', index_col='esr1_pos'), dict(zip(csv('results/numbering_map_esr1_esr2.csv').esr1_pos, csv('results/numbering_map_esr1_esr2.csv').same_residue)))` → `(48, 0, 7, 29)` (expected (48, 0, 7, 29)): **MATCH**
- **tl.py check:** PASS  chk csv results/**/consistent_positions*.csv rows>=1 cols=position,fraction_pairs  [246 rows (need 1)]
- **Provenance:** tracker/prov/E10.json exists; 2 ev files, no drift
- **Caveats:** 'Most pairs' = more than 50% of pairs, with at least 20 observed pairs: thresholds chosen in the script.

### R06: Entropy with and without structures failing R-free
- **What was done:** Compared entropy and the 90% number with and without structures that fail R-free. Re-grouped the passing set from the saved matrix.
- **Code:** scripts/entropy_qc_compare.py:main L37-69; scripts/entropy_qc_compare.py:n_90 L30-34. Run: `PYTHONPATH=. .venv/Scripts/python scripts/entropy_qc_compare.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows); results/qc_tiers.csv (CSV, 526 rows)
- **Output:** results/entropy_qc_compare.csv (CSV, 12 rows, sha256 3664938d): cut,subset,n_chains,n_entries,n_chains_left_out,n_clusters,entropy_bits,n_90,n_90_fraction… / 1.0,all,788,351,0,16,1.595844356895745,520,0.6598984771573604,0.0,0.0 / 1.5,all,788,351,0,7,1.1332722559505624,488,0.6192893401015228,0.0,0.0 / 2.0,all,788,351,0,4,1.0581080713180215,315,0.399746192893401,0.0,0.0
- **Log line:** `2026-10-20T22:25+05:00 \| R06 \| DONE \| who=opus \| chk=1/1 \| ev=results/entropy_qc_compare.csv,scripts/entropy_qc_compare.py \| note=R-free<=0.30 set (767 chains/343 entries, 21 chains out) reclustered from S04 slice: entropy within 0.06 bits of all at every cut, clusters 17/6/4/3 vs 16/7/4/3; tier_A row added (entropy delta <= 0.07 bits). n_90 shifts (e.g. 1.5 A 0.62 -> 0.47 of chains) because n_90 depends on which singleton clusters remain; QC does not change the redundancy conclusion`
- **Recount:** claim "R-free set 767 chains / 343 entries, 21 out; clusters 17/6/4/3"; command `(lambda d, m: (int(m.sum()), d.entry_id[m].nunique(), int((~m).sum()), [nclust(m, t) for t in (1.0, 1.5, 2.0, 2.5)]))(csv('results/entry_clusters.csv'), csv('results/entry_clusters.csv').entry_id.map(csv('results/qc_tiers.csv').set_index('entry_id').r_free).le(0.30).to_numpy())` → `(767, 343, 21, [17, 6, 4, 3])` (expected (767, 343, 21, [17, 6, 4, 3])): **MATCH**
- **tl.py check:** PASS  chk csv results/entropy_qc_compare.csv rows>=2  [12 rows (need 2)]
- **Provenance:** tracker/prov/R06.json exists; 2 ev files, no drift
- **Caveats:** A tier_A row was added beyond the spec. No tier A+B row (the S13 file has it). n_90 is unstable when clusters are few.

### K07: One-chain-per-entry view next to the all-chain view
- **What was done:** Took one chain per entry by a stated rule and re-grouped them. Compared that with the all-chain grouping.
- **Code:** scripts/entry_level_clusters.py:main L30-65. Run: `PYTHONPATH=. .venv/Scripts/python scripts/entry_level_clusters.py`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows); results/chain_status.csv (CSV, 901 rows); results/helix_bfactor_by_chain.csv (CSV, 788 rows)
- **Output:** results/entry_level_clusters.csv (CSV, 351 rows, sha256 a4a6ccf0): entry_id,chain_key,mapped_lbd_ca,mean_b_chain,n_chains_in_entry,cut_1.0,chain_level_cut_1.… / 1A52,1A52\|B\|1,240,54.611375,2,1,2,1,1,1,1,1,1 / 1ERE,1ERE\|D\|1,235,56.257021276595744,6,16,15,4,7,3,4,3,3 / 1ERR,1ERR\|A\|1,225,53.946044444444446,2,6,11,2,4,2,2,2,2 ‖ results/entry_vs_chain_clusters.csv (CSV, 4 rows, sha256 fdc032c2): cut,n_chains,n_entries,n_clusters_chain,n_clusters_entry,top3_chain,top3_entry,ari / 1.0,788,351,16,18,0.9543147208121827,0.8433048433048433,0.8407022832014079 / 1.5,788,351,7,6,0.9885786802030457,0.9914529914529915,0.9736325239143379 / 2.0,788,351,4,3,0.9974619289340102,1.0,0.9815305853014417
- **Log line:** `2026-10-20T22:25+05:00 \| K07 \| DONE \| who=opus \| chk=2/2 \| ev=results/entry_level_clusters.csv,results/entry_vs_chain_clusters.csv,scripts/entry_level_clusters.py \| note=351 entries, one chain each (most mapped LBD CA, then lowest mean CA B, then chain_key), S04 slice; entry vs chain clusters 18/6/3/3 vs 16/7/4/3, top-3 share entry 0.84/0.99/1.00/1.00 vs chain 0.95/0.99/1.00/1.00, ARI 0.84/0.97/0.98/0.99; both views kept`
- **Recount:** claim "351 entries; clusters 18/6/3/3; ARI 0.84/0.97/0.98/0.99"; command `(lambda d: (len(d), bool(d.entry_id.is_unique), [d[c].nunique() for c in CUTS], [ari(d['chain_level_' + c], d[c]) for c in CUTS]))(csv('results/entry_level_clusters.csv'))` → `(351, True, [18, 6, 3, 3], [0.840702, 0.973633, 0.981531, 0.988673])` (expected (351, True, [18, 6, 3, 3], [0.84, 0.97, 0.98, 0.99])): **MATCH**
- **tl.py check:** PASS  chk csv results/entry_level_clusters.csv rows>=350 cols=entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5 nonempty=entry_id,chain_key uniq=entry_id  [351 rows (need 350)]; PASS  chk csv results/entry_vs_chain_clusters.csv rows>=4 cols=cut,n_clusters_chain,n_clusters_entry,top3_chain,top3_entry,ari  [4 rows (need 4)]
- **Provenance:** tracker/prov/K07.json exists; 3 ev files, no drift
- **Caveats:** none found

### S15: Smallest set of representative structures
- **What was done:** Chose one representative per cluster until 90% of clusters are covered. Each pick has a one-line reason.
- **Code:** scripts/representative_set.py:main L25-64. Run: `PYTHONPATH=. .venv/Scripts/python scripts/representative_set.py`
- **Input:** results/medoids.csv (CSV, 7 rows); results/entry_clusters.csv (CSV, 788 rows); results/rmsd_entry.npy (NumPy, shape (788, 788) float64)
- **Output:** results/representative_set.csv (CSV, 7 rows, sha256 da8b9231): rank,cluster,entry_id,chain_key,qc_label,qc_pass,cluster_size_chains,cluster_entries,mean_… / 1,7,5WGD,5WGD\|B\|1,A,True,445,221,0.45123417387434506,False,True,"cluster 7 at 1.5 A (445 c… / 2,4,2QE4,2QE4\|B\|1,A,True,329,128,0.5039134183956449,True,True,"cluster 4 at 1.5 A (329 cha… / 3,1,2JF9,2JF9\|A\|1,A,True,5,3,0.571782326822496,True,True,"cluster 1 at 1.5 A (5 chains, ra…
- **Log line:** `2026-10-20T22:26+05:00 \| S15 \| DONE \| who=opus \| chk=1/1 \| ev=results/representative_set.csv,scripts/representative_set.py \| note=diversity = clusters (as R02/R04): 7 clusters at 1.5 A, 90% needs 7 -> 7 representatives, 6/7 QC tier A; cluster 2 (1R5K, 3 chains) has no tier-A chain, medoid 1R5K\|B\|1 (tier B) kept and flagged; cluster 7 medoid 5TM1\|A\|1 is tier B -> 5WGD\|B\|1 most central tier-A. By chain mass the top 2 already cover 98.2% (cum_chain_fraction column). Definition choice worth a human look`
- **Recount:** claim "7 representatives, 6 QC tier A; cluster 2 keeps 1R5K|B|1; top 2 clusters cover 98.2% of chains"; command `(lambda r, s: (len(r), int(r.qc_pass.sum()), r.chain_key[r.cluster == 2].item(), s.head(2).sum() / s.sum()))(csv('results/representative_set.csv'), csv('results/entry_clusters.csv')['cut_1.5'].value_counts())` → `(7, 6, '1R5K\|B\|1', 0.982234)` (expected (7, 6, '1R5K\|B\|1', 0.982)): **MATCH**
- **tl.py check:** PASS  chk csv results/**/representative*.csv rows>=1 cols=entry_id,justification  [7 rows (need 1)]
- **Provenance:** tracker/prov/S15.json exists; 2 ev files, no drift
- **Caveats:** Definition choice: '90% of diversity' read as 90% of clusters, so all 7 are needed; by chain share 2 would do. Needs a human decision.

### K08: Plain note for the professor on what data is kept
- **What was done:** Wrote a plain note for the professor on what used to be removed and what is kept now. It answers notes 1 and 2.
- **Code:** none (text task). Run: `none (text written by hand)`
- **Input:** results/chain_status.csv (CSV, 901 rows); results/attrition_funnel.csv (CSV, 16 rows)
- **Output:** docs/data_retention.md (md, 70 lines, sha256 1c580842): # Data retention: what was removed before, what is kept now / This note answers the professor's notes of 2026-10-18. Every number comes from a results f… / ## In short
- **Log line:** `2026-10-20T22:26+05:00 \| K08 \| DONE \| who=opus \| chk=3/3 \| ev=docs/data_retention.md,results/attrition_funnel.csv,results/chain_status.csv,results/qc_clustering_delta.csv,results/entropy_qc_compare.csv,results/entry_vs_chain_clusters.csv,results/helix_stability_by_bgroup.csv,results/h12_core_ratio_by_bgroup.csv,results/core_sensitivity.csv,results/h12_metric.csv,results/cluster_stats.csv,results/archive_v2_strict/cluster_stats.csv \| note=professor note: old strict vs keep-all, former filters (K01), QC labels (K05/K06), strict subsets (S13/R06), entry vs chain (K07), answers to notes 1 and 2; every number from an ev file`
- **Recount:** claim "244 chains with a warning, 113 outside the matrix; figures quoted in the note"; command `(lambda s, t: (int((s.in_matrix & (s.warning != 'none')).sum()), int((~s.in_matrix).sum()), all(x in t for x in ['244', '113', '1.255', '0.988'])))(csv('results/chain_status.csv'), txt('docs/data_retention.md'))` → `(244, 113, True)` (expected (244, 113, True)): **MATCH**
- **tl.py check:** PASS  chk lines docs/data_retention.md >=25  [44 lines (need 25)]; PASS  chk grep docs/data_retention.md (?i)label  [docs\data_retention.md]; PASS  chk grep docs/data_retention.md (?i)stable  [docs\data_retention.md]
- **Provenance:** tracker/prov/K08.json exists; 12 ev files, no drift
- **Caveats:** Text task; only these figures are recounted here.

### X02: Distance matrix between all ESR2 structures
- **What was done:** Aligned all 86 ESR2 chains to one ESR2 reference and measured every pairwise distance. Saved the matrix.
- **Code:** scripts/esr2_rmsd.py:main L35-83. Run: `PYTHONPATH=. .venv/Scripts/python scripts/esr2_rmsd.py`
- **Input:** results/lbd_chains.csv (CSV, 987 rows); results/esr2_entries.csv (CSV, 39 rows); results/lbd_extract/residue_mapping.csv (CSV, 222348 rows)
- **Output:** results/rmsd_esr2.npy (NumPy, shape (86, 86) float64, sha256 ec76852c): shape (86, 86), dtype float64 ‖ results/rmsd_esr2_index.txt (txt, 86 lines, sha256 8b1a9608): 1L2J\|A\|1 / 1L2J\|B\|1 / 1NDE\|A\|1 ‖ results/esr2_reference.json (JSON, 14 top-level keys, sha256 388cffad): { /   "chain_key": "3OLL\|A\|1", /   "pdb_id": "3OLL",
- **Log line:** `2026-10-20T22:28+05:00 \| X02 \| DONE \| who=opus \| chk=1/1 \| ev=results/rmsd_esr2.npy,results/rmsd_esr2_index.txt,results/esr2_alignment.csv,results/esr2_reference.json,scripts/esr2_rmsd.py \| note=86 eligible ESR2 chains (39 entries), all in matrix, ESR1 numbering, core v3 fit; reference 3OLL\|A\|1 (1.5 A). DEVIATION: reference_min_core_fraction 1.0 not applied (unreachable for ESR2, best 0.976), recorded in esr2_reference.json - human may want to confirm. Core RMSD median 0.572 A; pairwise median 0.947 A (max 6.64); min overlap 178`
- **Recount:** claim "86 chains / 39 entries; reference 3OLL|A|1; pairwise median 0.947 A (max 6.64); core median 0.572 A"; command `(lambda m, ix, al: (m.shape[0], bool(np.allclose(m, m.T)), len({k.split('|')[0] for k in ix}), float(np.median(up(m))), float(up(m).max()), al.core_rmsd_A.median(), json.loads(txt('results/esr2_reference.json'))['chain_key']))(npy('results/rmsd_esr2.npy'), txt('results/rmsd_esr2_index.txt').split(), csv('results/esr2_alignment.csv'))` → `(86, True, 39, 0.947077, 6.64165, 0.571783, '3OLL\|A\|1')` (expected (86, True, 39, 0.947, 6.64, 0.572, '3OLL\|A\|1')): **MATCH**
- **tl.py check:** PASS  chk npy results/rmsd_esr2.npy square>=30  [shape (86, 86) (need square >= 30)]
- **Provenance:** tracker/prov/X02.json exists; 5 ev files, no drift
- **Caveats:** Stated deviation: the reference core-fraction minimum (1.0) was not applied because no ESR2 chain reaches it (best 0.976).

### X15: Plot of ChEMBL potency per structural cluster
- **What was done:** Plotted ChEMBL potency of PDB ligands grouped by structural cluster. Potency is on a log scale.
- **Code:** scripts/potency_by_cluster.py:ligand_table L35-61; scripts/potency_by_cluster.py:figure L85-112. Run: `PYTHONPATH=. .venv/Scripts/python scripts/potency_by_cluster.py`
- **Input:** results/chain_ligand_of_interest.csv (CSV, 810 rows); results/chembl_pdb_matches.csv (CSV, 166 rows); data/chembl_esr1_activities.csv (CSV, 6714 rows); results/entry_clusters.csv (CSV, 788 rows)
- **Output:** figures/potency_by_cluster.png (PNG, 1125x675 px, sha256 9285151d): binary/structured (see shape) ‖ results/potency_by_cluster.csv (CSV, 82 rows, sha256 0022131c): ligand_id,cluster,n_chains,n_entries,share_in_cluster,median_h12_rmsd_A,pchembl,n_molecule… / OGJ,1,2,1,1.0,28.622209323384325,9.52,3,14.0,1.0,0.30199517204020193 / GW5,2,5,2,0.6,9.840633706095382,6.38,1,8.0,1.0,416.8693834703355 / 29S,4,6,2,1.0,13.512473706453406,7.05,9,28.0,0.9091,89.12509381337459
- **Log line:** `2026-10-20T22:29+05:00 \| X15 \| DONE \| who=opus \| chk=1/1 \| ev=figures/potency_by_cluster.png,results/potency_by_cluster.csv,scripts/potency_by_cluster.py \| note=82 PDB ligands with ChEMBL potency (median pChEMBL, IC50/Ki/Kd '='), unit = ligand assigned to the cluster holding most of its chains at 1.5 A: cluster 4 n=53, 7 n=27, 1 and 2 n=1; log nM axis, medians marked`
- **Recount:** claim "82 ligands; cluster 4 n=53, 7 n=27, 1 and 2 n=1"; command `(lambda d: (len(d), d.cluster.value_counts().to_dict(), nonwhite('figures/potency_by_cluster.png') > 0))(csv('results/potency_by_cluster.csv'))` → `(82, {4: 53, 7: 27, 1: 1, 2: 1}, True)` (expected (82, {4: 53, 7: 27, 1: 1, 2: 1}, True)): **MATCH**
- **tl.py check:** PASS  chk file figures/**/potency_by_cluster*.png  [figures\potency_by_cluster.png]
- **Provenance:** tracker/prov/X15.json exists; 3 ev files, no drift
- **Caveats:** Unit = one PDB ligand placed in the cluster holding most of its chains; assays are mixed (IC50/Ki/Kd).

### X16: Statistical tests: potency against cluster
- **What was done:** Tested whether potency differs between clusters and tracks the H12 shift. One table with n and p.
- **Code:** scripts/potency_by_cluster.py:tests L64-82. Run: `PYTHONPATH=. .venv/Scripts/python scripts/potency_by_cluster.py`
- **Input:** results/potency_by_cluster.csv (CSV, 82 rows)
- **Output:** results/potency_tests.csv (CSV, 3 rows, sha256 cadd948e): test,groups,n,statistic,p_value,median_a,median_b,cut,unit / kruskal_wallis,"4,7",80,16.115631667335084,5.9590145795987964e-05,,,1.5,"PDB ligand (media… / mann_whitney,4 vs 7,80,1110.0,6.088870798976637e-05,8.305,7.11,1.5,"PDB ligand (median pCh… / spearman_pchembl_vs_h12_rmsd,all,79,0.45412992312222167,2.6252530874850246e-05,,,1.5,"PDB …
- **Log line:** `2026-10-20T22:29+05:00 \| X16 \| DONE \| who=opus \| chk=1/1 \| ev=results/potency_tests.csv,results/potency_by_cluster.csv,scripts/potency_by_cluster.py \| note=per-ligand tests: Kruskal clusters 4,7 n=80 p 6.0e-5; Mann-Whitney 4 vs 7 median pChEMBL 8.31 vs 7.11 p 6.1e-5; Spearman pChEMBL vs median H12 RMSD rho 0.45 p 2.6e-5 (n=79). Descriptive: mixed assays, ligand class likely confounds (H12-displaced antagonists more potent)`
- **Recount:** claim "Kruskal p 6.0e-5; median pChEMBL 8.31 vs 7.11; Spearman vs H12 RMSD 0.45 (n=79)"; command `(lambda d: (f'{kruskal(d.pchembl[d.cluster == 4], d.pchembl[d.cluster == 7]).pvalue:.1e}', d.pchembl[d.cluster == 4].median(), d.pchembl[d.cluster == 7].median(), spearmanr(d.dropna(subset=['median_h12_rmsd_A']).pchembl, d.dropna(subset=['median_h12_rmsd_A']).median_h12_rmsd_A)[0], int(d.median_h12_rmsd_A.notna().sum())))(csv('results/potency_by_cluster.csv'))` → `('6.0e-05', 8.305, 7.11, 0.45413, 79)` (expected ('6.0e-05', 8.31, 7.11, 0.45, 79)): **MATCH**
- **tl.py check:** PASS  chk csv results/**/potency_tests*.csv rows>=1 cols=test,p_value,n  [3 rows (need 1)]
- **Provenance:** tracker/prov/X16.json exists; 3 ev files, no drift
- **Caveats:** Descriptive: ligand class probably confounds the link (no ligand-class labels yet).

### X03: Group the ESR2 structures the same way as ESR1
- **What was done:** Grouped the ESR2 chains the same way as ESR1. Quality flags are copied for reading only.
- **Code:** scripts/esr2_clusters.py:main L22-36. Run: `PYTHONPATH=. .venv/Scripts/python scripts/esr2_clusters.py`
- **Input:** results/rmsd_esr2.npy (NumPy, shape (86, 86) float64); results/esr2_entries.csv (CSV, 39 rows)
- **Output:** results/esr2_clusters.csv (CSV, 86 rows, sha256 f5fa7b4c): entry_id,chain_key,cut_1.0,cut_1.5,cut_2.0,cut_2.5,qc_pass / 1L2J,1L2J\|A\|1,6,5,4,3,fail / 1L2J,1L2J\|B\|1,6,5,4,3,fail / 1NDE,1NDE\|A\|1,3,3,3,3,fail ‖ results/esr2_linkage.npy (NumPy, shape (85, 4) float64, sha256 e8c0d327): shape (85, 4), dtype float64
- **Log line:** `2026-10-20T22:31+05:00 \| X03 \| DONE \| who=opus \| chk=1/1 \| ev=results/esr2_clusters.csv,results/esr2_linkage.npy,scripts/esr2_clusters.py \| note=86 ESR2 chains/39 entries, S08 method on X02 matrix; clusters 7/5/4/3 at 1.0/1.5/2.0/2.5 A; sizes at 1.5 A 68/12/2/2/2; qc_pass copied, not used`
- **Recount:** claim "86 chains / 39 entries; clusters 7/5/4/3; sizes at 1.5 A 68/12/2/2/2"; command `(lambda d: (len(d), d.entry_id.nunique(), [d[c].nunique() for c in CUTS], d['cut_1.5'].value_counts().tolist(), int(fcluster(linkage(squareform(npy('results/rmsd_esr2.npy'), checks=False), 'average'), 1.5, 'distance').max())))(csv('results/esr2_clusters.csv'))` → `(86, 39, [7, 5, 4, 3], [68, 12, 2, 2, 2], 5)` (expected (86, 39, [7, 5, 4, 3], [68, 12, 2, 2, 2], 5)): **MATCH**
- **tl.py check:** PASS  chk csv results/esr2_clusters.csv rows>=30  [86 rows (need 30)]
- **Provenance:** tracker/prov/X03.json exists; 3 ev files, no drift
- **Caveats:** none found

### X04: Two-axis map (PCA) of the ESR2 structures
- **What was done:** Projected the ESR2 distance matrix onto two axes (PCA) and coloured points by cluster.
- **Code:** scripts/esr2_pca.py:main L30-65. Run: `PYTHONPATH=. .venv/Scripts/python scripts/esr2_pca.py`
- **Input:** results/rmsd_esr2.npy (NumPy, shape (86, 86) float64); results/esr2_clusters.csv (CSV, 86 rows)
- **Output:** figures/esr2_pca_landscape.png (PNG, 1050x780 px, sha256 f9e794db): binary/structured (see shape) ‖ results/esr2_pca_coords.csv (CSV, 86 rows, sha256 6875920a): chain_key,entry_id,pc1,pc2,cluster / 1L2J\|A\|1,1L2J,10.840103014718354,4.261051359706125,5 / 1L2J\|B\|1,1L2J,8.557893916057033,4.551535987240276,5 / 1NDE\|A\|1,1NDE,18.90940785700647,2.547370231903846,3
- **Log line:** `2026-10-20T22:32+05:00 \| X04 \| DONE \| who=opus \| chk=1/1 \| ev=figures/esr2_pca_landscape.png,results/esr2_pca_coords.csv,results/esr2_pca_variance.json,scripts/esr2_pca.py \| note=86 ESR2 chains, block-12 PCA of RMSD rows; PC1 81.9%, PC2 13.0%; coloured by X03 cluster (ESR2 ligand classes not labelled in repo); main basin 68 chains, second 12; 1 entry (2 chains) far out at PC1 ~40`
- **Recount:** claim "PC1 81.9%, PC2 13.0%"; command `tuple(PCA(n_components=2).fit(npy('results/rmsd_esr2.npy')).explained_variance_ratio_.tolist())` → `(0.818628, 0.129961)` (expected (0.819, 0.13)): **MATCH**
- **tl.py check:** PASS  chk file figures/esr2_pca_landscape.png  [figures\esr2_pca_landscape.png]
- **Provenance:** tracker/prov/X04.json exists; 4 ev files, no drift
- **Caveats:** Spec says 'same method as S11', which colours by ligand class; ESR2 has no ligand-class labels, so the figure colours by structural cluster.

### X05: ESR2 cluster sizes, with small clusters marked not interpretable
- **What was done:** Counted entries per ESR2 cluster. Clusters with fewer than three entries are marked not interpretable.
- **Code:** scripts/esr2_basin_counts.py:main L19-29. Run: `PYTHONPATH=. .venv/Scripts/python scripts/esr2_basin_counts.py`
- **Input:** results/esr2_clusters.csv (CSV, 86 rows)
- **Output:** results/esr2_basin_counts.csv (CSV, 5 rows, sha256 6d7ae3ce): basin,n_chains,n_entries,entries,interpretable,caveat,cut / 2,68,29,1U3Q;1U3R;1U3S;1U9E;1X76;1X78;1X7B;1X7J;1YY4;1YYE;1ZAF;2NV7;2YJD;2YLY;3OLL;3OLS;3O… / 5,12,7,1L2J;1QKM;2GIU;2I0G;2JJ3;2QTU;2Z4B,True,,1.5 / 1,2,1,2FSZ,False,fewer than 3 entries: not interpretable,1.5
- **Log line:** `2026-10-20T22:32+05:00 \| X05 \| DONE \| who=opus \| chk=1/1 \| ev=results/esr2_basin_counts.csv,scripts/esr2_basin_counts.py \| note=5 ESR2 basins at 1.5 A: 2 interpretable (29 and 7 entries), 3 not interpretable (2FSZ 1 entry; 1NDE+2YLY 2 entries; 5TOA 1 entry), kept and marked`
- **Recount:** claim "5 basins with 29/7/2/1/1 entries; 2 interpretable"; command `(lambda n: (sorted(n.tolist(), reverse=True), int((n >= 3).sum())))(csv('results/esr2_clusters.csv').groupby('cut_1.5').entry_id.nunique())` → `([29, 7, 2, 1, 1], 2)` (expected ([29, 7, 2, 1, 1], 2)): **MATCH**
- **tl.py check:** PASS  chk csv results/esr2_basin_counts.csv rows>=1 cols=basin,n_entries,interpretable  [5 rows (need 1)]
- **Provenance:** tracker/prov/X05.json exists; 2 ev files, no drift
- **Caveats:** none found

### W10: Manuscript section: dataset and quality control
- **What was done:** Wrote the dataset and quality-control section of the manuscript. Numbers are taken from result files.
- **Code:** none (text task). Run: `none (text written by hand)`
- **Input:** results/attrition_funnel.csv (CSV, 16 rows); results/entries.csv (CSV, 526 rows); results/lbd_chains.csv (CSV, 987 rows)
- **Output:** docs/manuscript/02_dataset_qc.md (md, 71 lines, sha256 2bc9be74): # Dataset and quality control / ## Methods / ### Structure collection
- **Log line:** `2026-10-20T22:33+05:00 \| W10 \| DONE \| who=opus \| chk=1/1 \| ev=docs/manuscript/02_dataset_qc.md,data/pdb_ids_esr1.txt,data/pdb_ids_esr2.txt,data/raw/download_manifest.csv,results/entries.csv,results/attrition_funnel.csv,results/lbd_chains.csv,results/chain_status.csv,results/qc_tiers.csv,results/ligands.csv,results/ligands_ranked.csv,results/chain_ligand_of_interest.csv,results/near_duplicates.csv,results/representative_set.csv,results/entropy.csv,results/qc_outliers_stats.csv,results/qc_bfactor_clusters.csv,results/qc_bfactor_clusters_tests.json,results/qc_clustering_delta.csv,results/entropy_qc_compare.csv \| note=methods + results for Q01-Q18 incl. QC questions (Q14, Q15) and QC sensitivity; keep-all framing (788 chains, labels not removal); every number from an ev file`
- **Recount:** claim "1,666 chain-models; 247 selected entries pass all three QC flags; >= 25 lines"; command `(lambda t, f, e, c: ('1,666' in t, int(f.n_before.iloc[0]), int(e[e.entry_id.isin(set(c[(c.identity == 'esr1') & c.eligible_historical].pdb_id))].qc_pass.eq('pass').sum()), '247' in t, len([l for l in t.splitlines() if l.strip()]) >= 25))(txt('docs/manuscript/02_dataset_qc.md'), csv('results/attrition_funnel.csv'), csv('results/entries.csv'), csv('results/lbd_chains.csv'))` → `(True, 1666, 247, True, True)` (expected (True, 1666, 247, True, True)): **MATCH**
- **tl.py check:** PASS  chk lines docs/manuscript/02_dataset_qc.md >=25  [38 lines (need 25)]
- **Provenance:** tracker/prov/W10.json exists; 20 ev files, no drift
- **Caveats:** Text task; only these figures are recounted here.

### W09: Worked example notebook that rebuilds the main figure
- **What was done:** Made a notebook that loads saved results and redraws the main figure. It checks the 90% numbers against the saved table.
- **Code:** none (text task). Run: `python -m jupyter nbconvert --to notebook --execute notebooks/worked_example.ipynb --output worked_example.ipynb`
- **Input:** results/rmsd_entry.npy (NumPy, shape (788, 788) float64); results/entry_clusters.csv (CSV, 788 rows); results/diversity_curve.csv (CSV, 3152 rows); results/saturation.csv (CSV, 4 rows)
- **Output:** notebooks/worked_example.ipynb (notebook, 9 cells, sha256 eea839c9): binary/structured (see shape)
- **Log line:** `2026-10-20T22:34+05:00 \| W09 \| DONE \| who=opus \| chk=1/1 \| ev=notebooks/worked_example.ipynb,results/rmsd_entry.npy,results/entry_clusters.csv,results/diversity_curve.csv,results/saturation.csv \| note=notebook loads saved outputs only, asserts n_90 equals results/saturation.csv (520/488/315/147), redraws the diversity figure; executed with nbconvert from repo root (outputs saved in the notebook). README not edited (outside out); run command is in the notebook's first cell`
- **Recount:** claim "executed without error; n_90 values from saturation.csv appear in its outputs"; command `(lambda nb: (all(c.get('execution_count') for c in nb['cells'] if c['cell_type'] == 'code'), any(o.get('output_type') == 'error' for c in nb['cells'] for o in c.get('outputs', [])), all(str(x) in json.dumps(nb) for x in csv('results/saturation.csv').n_90)))(json.loads(txt('notebooks/worked_example.ipynb')))` → `(True, False, True)` (expected (True, False, True)): **MATCH**
- **tl.py check:** PASS  chk file notebooks/**/*.ipynb  [notebooks\worked_example.ipynb]
- **Provenance:** tracker/prov/W09.json exists; 5 ev files, no drift
- **Caveats:** The script that generated the notebook was a scratch file, not committed. README was not updated, although the spec says 'runs from README'.

### W16: Script that collects the headline numbers for checking the manuscript
- **What was done:** Extended the numbers script with the headline numbers of the current run. It writes them to numbers.json.
- **Code:** scripts/collect_numbers.py:keep_all L291-340; scripts/collect_numbers.py:write L343-360. Run: `.venv/Scripts/python scripts/collect_numbers.py`
- **Input:** results/*.csv (glob, 59 files); results/*.json (glob, 5 files)
- **Output:** results/numbers.json (JSON, 1 top-level keys, sha256 66ef1ffc): { /   "numbers": [ /     { ‖ docs/numbers.md (md, 279 lines, sha256 50cb6b38): query_esr1_entries: 487 (data/pdb_ids_esr1.txt) / query_esr2_entries: 39 (data/pdb_ids_esr2.txt) / entries_master_table: 526 (results/entries.csv)
- **Log line:** `2026-10-20T22:35+05:00 \| W16 \| DONE \| who=opus \| chk=1/1 \| ev=results/numbers.json,docs/numbers.md,scripts/collect_numbers.py \| note=collect_numbers.py extended with keep_all(): 64 new headline numbers (chains/matrix/tiers, clusters, top-3, entropy, n_90 per cut, QC-subset and entry-level ARIs, core sensitivity, Q14/Q15/S07/X16 p-values and n, ESR2 X02/X05); 279 total; the earlier 215 values unchanged (diffed). Needs Bio, run with .venv python`
- **Recount:** claim "279 numbers, 64 new keep-all numbers"; command `(lambda n: (len(n), sum(x['name'].startswith(('keep_', 'q14_', 'q15_', 's07_', 'x16_', 'x02_', 'x05_')) for x in n)))(json.loads(txt('results/numbers.json'))['numbers'])` → `(279, 64)` (expected (279, 64)): **MATCH**
- **tl.py check:** PASS  chk file results/numbers.json  [results\numbers.json]
- **Provenance:** tracker/prov/W16.json exists; 3 ev files, no drift
- **Caveats:** Needs Biopython (imports ligand_of_interest), so it runs only with the .venv Python.

## Invariants

- **FAIL** archive hashes equal MANIFEST (raw bytes): 6/19 byte-identical; differ: rmsd_entry_index.csv, entry_clusters.csv, cluster_stats.csv, medoids.csv, qc_clustering_delta.csv, h12_metric.csv, alignment_core_rmsd.csv, limitations.md, alignment/alignment_transforms.json, alignment/global_alignment.csv, alignment/chain_keys.txt, alignment/rmsd_entry_index.txt, alignment/distance_exclusions.json. Command: `scripts/check_evidence.py:inv_archive_raw`
- **PASS** archive hashes equal MANIFEST (after CRLF->LF): 19/19 match after CRLF->LF; differ: none. Command: `scripts/check_evidence.py:inv_archive_lf`
- **PASS** rmsd_entry.npy 788x788, no NaN, symmetric; index 788 unique: shape (788, 788), NaN 0, symmetric True, index rows 788, unique 788. Command: `scripts/check_evidence.py:inv_matrix`
- **PASS** chain-level tables have >= 788 rows: results/entry_clusters.csv 788; results/rmsd_entry_index.csv 788; results/alignment_core_rmsd.csv 788; results/h12_metric.csv 788; results/helix_bfactor_by_chain.csv 788; results/chain_status.csv 901; results/qc_outliers.csv 788; results/alignment/global_alignment.csv 788; results/chain_ligand_of_interest.csv (unique chain_key) 788. Command: `scripts/check_evidence.py:inv_chain_tables`
- **PASS** clustering files differ from archive; entry_clusters 788 rows: identical to archive: none; entry_clusters rows 788. Command: `scripts/check_evidence.py:inv_clustering_changed`
- **FAIL** S13 and R06 show all, tier A and tier A+B; main file uses all chains: S13 subsets ['all', 'tier_A', 'tier_A+B']; R06 subsets ['all', 'rfree<=0.30', 'tier_A']; main file all 788 chains True. Command: `scripts/check_evidence.py:inv_subsets`
- **FAIL** no archived file lost rows: rmsd_entry.npy 351->788; rmsd_entry_index.csv 351->788; entry_clusters.csv 351->788; linkage_entry.npy 350->787; cluster_stats.csv 4->4; medoids.csv 9->7; qc_clustering_delta.csv 4->12; h12_metric.csv 351->788; alignment_core_rmsd.csv 788->788; alignment/global_alignment.csv 788->788; alignment/chain_keys.txt 788->788; alignment/rmsd_entry_index.txt 351->788; alignment/overlap_entry.npy 351->788; lost rows: medoids.csv. Command: `scripts/check_evidence.py:inv_no_row_loss`
- **PASS** H12 present: 788 rows, labels add up to 788: rows 788; {'complete': 553, 'partial': 212, 'not_modelled': 23}. Command: `scripts/check_evidence.py:inv_h12`
- **PASS** every ESR2 entry appears in the ESR2 outputs: 39 ESR2 entries; missing: none. Command: `scripts/check_evidence.py:inv_esr2`
- **PASS** tl.py validate 0 errors; DONE tasks without prov: 125 tasks in 1 file(s), 140 log lines, 0 errors, 59 warnings; DONE tasks without prov file (23): C04, E03, E04, Q01, Q02, Q03, Q04, Q05, Q06, Q07, Q08, Q09, Q12, Q13, S01, W02, W03, W04, W07, W14, X01, X13, X14. Command: `scripts/check_evidence.py:inv_validate`

Notes on the invariants:
- Raw-byte archive hashes differ only because git on this Windows clone checks text files out with CRLF line endings (`git config core.autocrlf` = true, set in C:/Program Files/Git/etc/gitconfig). After CRLF->LF all 19 archived files equal MANIFEST.csv.
- medoids.csv has 7 rows instead of 9: one row per cluster, and the 788-chain run has 7 clusters at 1.5 A against 9 in the strict run. No chain is dropped; the row unit is the cluster.
- R06 compares all chains, the R-free <= 0.30 set and tier A; it has no tier A+B row (its task spec asks only for the R-free set). Tier A+B is in the S13 file.
- DONE tasks without a prov file were logged before tracker/prov/ existed; none of them is in this job.
